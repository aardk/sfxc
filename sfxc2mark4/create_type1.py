#!/usr/bin/env python
import sys
import struct
import argparse
import json
import time
import os
import numpy as np
from collections import namedtuple
from datetime import datetime, timedelta
from urlparse import urlparse
from sfxcdata import SFXCData
from experiment import experiment
from stationmap import create_one_letter_mapping
from ovex import create_rootfile, create_global_ovex, get_nbits
from vex import Vex 

def create_root_id(t):
  # rootid format changes after when the 'old' format with epoch in 1979 
  # rolls over, this is on Mon Feb 26 16:45:04 2018
  root_break = 1519659904

  # Sadly datetime doesn't know about leap seconds, so good 'ol POSIX time it is
  tnow = time.mktime(t.utctimetuple())

  result = ''
  if tnow >= root_break:
    # 'New' style root_id [0-9A-Z]*6
    dt = int(tnow - root_break)
    table = [str(i) for i in range(10)] + [chr(ord('A') + i) for i in range(26)]
    for i in range(6):
      result = table[dt%36] + result
      dt /= 36
  else:
    # 'Old' style root_id with epoch in 1979, [a-z]*6
    t79 = time.mktime(datetime(1979,1,1).utctimetuple()) # epoch

    # Time since epoch is in units of four seconds
    dt = int((tnow - t79) / 4)
    for i in range(6):
      result = chr('a' + dt%26) + result
      dt /= 26
  return result

def scale_string_float(value, scalefactor):
  # Scale floating point number stored as string without changing the precision
  n = len(value)
  v = value.lower().partition('.')
  expidx1 = v[0].find('e')
  expidx2 = v[2].find('e')
  predot = len(v[0][:expidx1])
  postdot = len(v[2][:expidx2])
  fmt = "{:%dg}"%(predot + postdot)
  return fmt.format(float(value) * scalefactor)

def fix_vex_for_hopps(vex):
  # Fix the vex file to comply with HOPPS
  
  # Add mk4_site_ID to $SITE
  for key in vex['STATION']:
    site = vex['STATION'][key]['SITE']
    st = vex['SITE'][site]['site_ID']
    vex['SITE'][site]['mk4_site_ID'] = STATIONMAP[st]

  # fix axis_offset in the $ANTENNA block (differs from the definition in vex 1.5)
  for key in vex['ANTENNA']:
    try:
      axis_offset = vex['ANTENNA'][key]['axis_offset']

      if type(axis_offset) != list:
        vex['ANTENNA'][key].pop('axis_offset')
        # difx2mark4 alway puts 'el' in the first field
        vex['ANTENNA'][key]['axis_offset'] = ['el', axis_offset]
    except KeyError:
      # Should we add axis_offset if it is missing?
      print 'No axis_offset for station', key
      pass
  
  # Fix the units for the clock offsets
  for key in vex['CLOCK']:
    for clock_early in vex['CLOCK'][key].getall('clock_early'):
      delay = clock_early[1].split()
      if (len(delay) == 2) and delay[1] == 'sec':
        delayusec = scale_string_float(delay[0], 1e6)
        clock_early[1] = "{} usec".format(delayusec)
      rate = clock_early[3].split()
      if len(rate) == 1:
        # SFXC defaults to usec / sec but HOPPS uses sec / sec
        scale = 1e-6
      else:
        # convert units to sec / sec
        unit = [x.strip() for x in ("".join(rate[1:])).split('/')]
        if len(unit) != 2:
          print 'Error, invalid unit for clock rate: "{}"'.format("".join(rate[1:]))
          exit(1)
        scale = 1e-6 if unit[0] == 'usec' else 1.
        if unit[1] == 'usec':
          scale *= 1e6
      ratesec = scale_string_float(rate[0], scale) 
      clock_early[3] = "{}".format(ratesec)

def hopsdate(t):
  yday = (t - datetime(t.year, 1, 1)).days + 1
  sec = t.second + t.microsecond / 1000000.
  return (t.year, yday, t.hour, t.minute, sec)

def create_t000(f):
  # type 000 header format:
  # <field>  <nbytes>  <description>
  # Type     3         Equal to 000
  # Version  2         Version number
  # Unused   3         Spaces
  # Date     16        " yyyyddd-hhmmss ", original file creation date.
  # Name     32        Filename relative to data directory, NULL-padded. 
  # Unused   8         NULLs.
  yday = (CREATIONDATE - datetime(CREATIONDATE.year, 1, 1)).days + 1
  t = "{:4}{:03}-{:02}{:02}{:02}".format(CREATIONDATE.year, yday, CREATIONDATE.hour, \
                                   CREATIONDATE.minute, CREATIONDATE.second)
  f.write(struct.pack("3s2s3s16s40s", "000", "01", "000", t, f.name[:31]))

def create_t100(f, bl, rootname, data):
  # type 100 header format:
  # <field>     <type> <nbytes> <description>
  # Type         ascii    3      100
  # Version      ascii    2      0-99
  # Unused       ascii    3      Spaces
  # procdate     date    12      Time of correlation
  # Baseline     ascii    2      Standard baseline id
  # rootname     ascii   34      Name of root, from exp. dir
  # Qcode        ascii    2      Correlation "quality", criteria TBD
  # Unused2      ascii    6      spaces for padding
  # Percent done r*4      4      0-100% of scheduled data processed
  #
  # start        date    12      Time of 1st AP
  # stop         date    12      Time of last AP
  # ndrec        i*4      4      Total no of data records in file
  # nindex       i*4      4      Number if index numbers present
  # nlags        i*2      2      Number of lags/type-120 record
  # nblocks      i*2      2      Number of blocks per index number

  # type, verion, unused
  f.write(struct.pack("3s2s3s", "100", "00", "000"))
  # procdate (creation date)
  ctime = datetime.utcfromtimestamp(os.stat(data.inputfile.name).st_ctime)
  yday = (ctime - datetime(ctime.year, 1, 1)).days + 1
  f.write(struct.pack('!4hf', ctime.year, yday, ctime.hour, ctime.minute, ctime.second))
  # NB setting Percent done to 0.0
  code = STATIONMAP[bl[0]] + STATIONMAP[bl[1]]
  f.write(struct.pack('!2s34s2s6sf', code, rootname, '\0\0', 6*' ', 0.0))
  # start date
  t = hopsdate(data.current_time())
  f.write(struct.pack("!4hf", *t))
  # stop date and ndrec will be filled in after the last AP was written
  f.write(struct.pack("!4hfi", 0, 0, 0, 0, 0.0, 0))
  # nindex, nlags, and nblocks
  nindex = len(data.vis[bl])
  nlags = data.global_header.nchan 
  f.write(struct.pack('!i2h', nindex, nlags, 1))

def create_t101(f, crosschan): 
  # type 101 header format:
  # <field>     <type> <nbytes> <description>
  # Type         ascii    3       101
  # Version      ascii    2       0-99 
  # Status       ascii    1       Currently unused, set to null
  # nblocks      i*2      2       Number of block table entries
  # Index        i*2      2       Index number
  # Primary      i*2      2       Index number of primary 101
  # Ref_chan_id  ascii    8       from vex, e.g. X1R
  # Rem_chan_id  ascii    8       from vex, e.g. X1L
  # Corr. board  i*2      2       Correlator board serial #
  # Corr. slot   i*2      2       Correlator board slot
  # Ref channel  i*2      2       SU output channel numbers
  # Rem channel  i*2      2       
  # Post mortem  i*4      4       Up to 32 1-bit flags
  # Block table  i*4      4*nblocks   One entry per block in snake

  # Type, version, status, nblocks
  f.write(struct.pack("!3s2ssh", "101", "00", "\0", 1))
  # Index, Primary Ref_chan_id, Rem_chan_id
  f.write(struct.pack("!2h8s8s", crosschan['index'], 0, crosschan['ref_chan_id'], crosschan['rem_chan_id']))
  # Remaining t101 fields are obsolete, set to zero
  f.write(struct.pack("!4h2i", *(6*[0])))

def write_t120s(data, scan, outfiles, t101_map, scalefactors, ap):
  # type 120 header format:
  # <field>     <type> <nbytes> <description>
  # Type        ascii     3       120
  # Version     ascii     2       0-99 
  # corrtype    ascii     1       Binary 1-4, lagdata struct type
  # nlags       i*2       2       Number of lags present
  # Baseline    ascii     2       Standard baseline id
  # rootcode    ascii     6       Standard root suffix
  # Index       i*4       4       Index into type 101 record
  # AP          i*4       4       Accumulation period number
  # flag_wgt    i*4/f*4   4       Lag data: Up to 32 correlation flags, 
  #                               Spectral: weight 0.0 - 1.0
  # Status      i*4       4       Up to 32 status bits
  # fr_delay    i*4       4       Mid-AP fractional delay (bits * 2^32)
  # delay_rate  i*4       4       Mid-AP delay rate (bits/sysclk * 2^32)
  # lagdata     array   variable  Correlation counts

  # Type, vesion, corrtype (equal to 5 for SFXC / DiFX), nlags
  integr_time = data.integration_time.total_seconds()
  hbase = struct.pack("!3s2sbh", "120", "00", 5, data.nchan)
  for bl in data.vis:
    if (ap * integr_time) >= min(scan['data_stop'][bl[0]], scan['data_stop'][bl[1]]):
      continue
    # Skip initial integrations which have no data 
    #TODO: should we do a weight cut-off for all visibilities?
    if (outfiles[bl]['nrec'] == 0) and (data.vis[bl].items()[0][1].weight == 0):
        continue

    if bl[0] == bl[1]:
      scalefac = 10000.
    else:
      scalefac = scalefactors[bl[0]] * scalefactors[bl[1]]
    blname = "{}{}".format(STATIONMAP[bl[0]], STATIONMAP[bl[1]])
    f = outfiles[bl]['file']
    outfiles[bl]['nrec'] += len(data.vis[bl])
    for ch in data.vis[bl]:
      f.write(hbase)
      # Baseline, rootcode, Index, ap
      t101 = t101_map[ch]
      f.write(struct.pack("!2s6s2i", blname, ROOTID, t101['index'], ap))
      # Normalize weight
      bw = scan['freq'][(ch.freqnr, ch.sideband, ch.pol1)]['bw']
      v = data.vis[bl][ch]
      weight = v.weight * 1e-6 / (2 * bw * integr_time)
      # Weight, Status,  fr_delay, delay_rate 
      f.write(struct.pack('!f3i', weight, 0, 0, 0))
      # The visibility data itself
      (v.vis[:-1]*scalefac).byteswap().tofile(f)

def create_t1map(scan, data):
  channels = {j:i for i,j in enumerate(sorted([x[:2] for x in scan["freq"]]))}
  crosschans = set()
  for bl in data.vis:
    for ch in data.vis[bl]:
      crosschans.add(ch)
  t1map = {}
  for ch in crosschans:
    key = (ch.freqnr, ch.sideband)
    # produce the same mapping as DiFX (LL = 1, RR = 2, LR = 3, RL = 4)
    index = channels[key]*10 + 2 + ch.pol1 + 2*ch.pol2 - 4*ch.pol1*ch.pol2
    t1map[ch] = {'index': index} 
    t1map[ch]['ref_chan_id'] = scan["freq"][(ch.freqnr, ch.sideband, ch.pol1)]['freq_id']
    t1map[ch]['rem_chan_id'] = scan["freq"][(ch.freqnr, ch.sideband, ch.pol2)]['freq_id']
  return t1map

def create_scaling_factors(vex, scan, data):
  # sqrt (10000 * Van Vleck correction)
  factors = [1.25331e2, 1.06448e2]
  nbits = get_nbits(vex, scan['mode'], data.stations)
  scalefactors = {}
  for n in nbits:
    stations = nbits[n]
    for station in stations:
      scalefactors[station] = factors[n-1]
  return scalefactors 

def initialise_next_scan(vex, exper, ctrl, data):
  scan = exper.get_scan(vex, ctrl, data)
  dirname = BASENAME + '/' + scan['name'].encode('ascii')
  if not os.path.exists(dirname):
    os.makedirs(dirname)
  rootname = create_rootfile(vex, ctrl, scan['name'], data.stations, dirname, ROOTID)
  scale = create_scaling_factors(vex, scan, data)
  t1map = create_t1map(scan, data)

  # create type 1 files
  outfiles = {}
  for bl in data.vis:
    code = STATIONMAP[bl[0]] + STATIONMAP[bl[1]]
    fname = dirname + '/' + code + '..' + ROOTID
    f = open(fname, 'w')
    outfiles[bl] = {'file': f, 'nrec': 0}
    create_t000(f)
    create_t100(f, bl, rootname, data)
    channels = sorted(data.vis[bl].keys())
    for ch in channels:
      create_t101(f, t1map[ch])
  return scan, outfiles, t1map, scale

def finalize_type1(outfiles, data, scan, nap):
 for bl in outfiles:
    f = outfiles[bl]['file']
    nrec = outfiles[bl]['nrec']
    stop = hopsdate(scan['start'] + nap * data.integration_time)
    # Write the remaining t100 fields: stop and ndrec
    f.seek(144)
    f.write(struct.pack('!4hfi', stop[0], stop[1], stop[2], stop[3], stop[4], nrec))
      
def parse_args():
  global CREATIONDATE
  CREATIONDATE = datetime.utcnow()
  parser = argparse.ArgumentParser(description='Convert SFXC output to mark4 format')
  parser.add_argument("vexfile", help='vex file')
  parser.add_argument("ctrlfile", help='SFXC control file used in the correlation')
  parser.add_argument('-c', "--create-ovex", help='Create toplevel ovex', action='store_true')
  parser.add_argument('-r', "--rootid", help='Manually specify rootid', default=create_root_id(CREATIONDATE))
  args = parser.parse_args()
  vex = Vex(args.vexfile)
  ctrl = json.load(open(args.ctrlfile, 'r'))
  return vex, ctrl, args.rootid, args.create_ovex

def process_job(vex, ctrl, rootid, create_ovex, basename="1234"):
  global ROOTID, BASENAME, STATIONMAP
  exper = experiment(vex)
  ROOTID = rootid
  BASENAME = basename
  stations = []
  for station_def in vex['STATION']:
    stations.append(station_def)
    continue
  stations.sort()
  sources = []
  try:
    scan = ctrl['scans'][0]
  except:
    start = ctrl['start']
    for scan in vex['SCHED']:
      if start == vex['SCHED'][scan]:
        break
  source = vex['SCHED'][scan]['source']
  sources.append(vex['SOURCE'][source]['source_name'])
  out_tuple = urlparse(ctrl['output_file'])
  output_file = out_tuple.netloc if out_tuple.path == '' else out_tuple.path
  data = SFXCData(output_file, stations, sources)
  STATIONMAP = create_one_letter_mapping(vex)
  fix_vex_for_hopps(vex)
  if create_ovex:
    try:
      setup_station = ctrl['setup_station']
    except:
      setup_station = ctrl['stations'][0]
    create_global_ovex(vex, setup_station)

  scan = {'stop': datetime(1,1,1)}
  outfiles = {}
  ap = 0
  while True:
    # Check if we need to move to next scan
    if data.current_time() >= scan["stop"]:
      # write end time and number of records for previous scan
      finalize_type1(outfiles, data, scan, ap)
      scan, outfiles, t101_map, scale = initialise_next_scan(vex, exper, ctrl, data)
      ap = 0

    write_t120s(data, scan, outfiles, t101_map, scale, ap)
    ap += 1

    if not data.next_integration():
      break
  # write end time and number of records for final scan
  finalize_type1(outfiles, data, scan, ap)
#########
########################## MAIN #################################3
########
if __name__ == "__main__":
  vex, ctrl, rootid, create_ovex = parse_args()
  process_job(vex, ctrl, rootid, create_ovex)
