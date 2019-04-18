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
from vex import Vex 

def ensure_list(x):
  if type(x) == list:
    return x
  return [x]

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

def get_nbits(vex, mode, stations):
  nbits = {}
  for station in stations:
    for das in vex['STATION'][station].getall('DAS'):
      try:
        rack_type = vex['DAS'][das]['electronics_rack_type']
      except KeyError:
        pass
      try:
        transport_type = vex['DAS'][das]['record_transport_type']
      except KeyError:
        pass
    datatype = 'VDIF'
    if (transport_type == 'Mark5A'):
      datatype = 'Mark5A'
    if (transport_type == 'Mark5B') and (rack_type not in ("DVP", "RDBE2", "WIDAR")):
      datatype = 'Mark5B'

    n = 0
    if datatype == 'Mark5B': 
      # Mark5B data
      for x in vex['MODE'][mode].getall('BITSTREAMS'):
        if station in x[1:]:
          n = 1
          for line in vex['BITSTREAMS'][x[0]].getall('stream_def'):
            if line[1] == 'mag':
              n = 2
              break
      if n == 0:
        print 'Warning: No $BITSTREAMS section for station: "{}", trying $TRACKS section'.format(station)
    elif datatype == 'VDIF':
      # VDIF data
      for x in vex['MODE'][mode].getall('THREADS'):
        if station in x[1:]:
          n = int(vex['THREADS'][x[0]]['thread'][5])
          break
      if n == 0:
        print 'Error: No THREADS section for station: "{}", trying $TRACKS section (is it really VDIF?)'.format(station)
    
    if n == 0:
      # Mark5A data or $BITSTREAMS/$THREADS not found
      for x in vex['MODE'][mode].getall('TRACKS'):
        if station in x[1:]:
          n = 1
          for line in vex['TRACKS'][x[0]].getall('fanout_def'):
            if line[2] == 'mag':
              n = 2
              break
    if n == 0:
      print 'Error: No TRACKS section for station: "{}"'.format(station)
      exit(1)
    # store result for station
    try:
      nbits[n].append(station)
    except KeyError:
      nbits[n] = [station]
  return nbits

def create_root(vex, ctrl, scan, data, dirname):
  rootname = dirname + '/' + data.source + '.' + ROOTID
  f = open(rootname, 'w')
  # $Head and $GLOBAL
  f.write('$OVEX_REV;\nrev = 1.5;\n$GLOBAL;\n')
  for i in ('EOP', 'EXPER'):
    f.write('    ref ${} = {};\n'.format(i, vex['GLOBAL'][i]))
  
  # $EXPER
  name = vex['GLOBAL']['EXPER']
  f.write('*\n$EXPER;\n  def {};\n'.format(name))
  for i in ('exper_name', 'exper_description', 'PI_name'):
    f.write('    {} = {};\n'.format(i, vex['EXPER'][name][i]))
  # Target correlator has to difx, otherwise HOPS assumes it is lagdata
  f.write('    target_correlator = difx;\n  enddef;\n')

  # $MODE
  try:
    setup_station = ctrl['setup_station']
  except:
    setup_station = ctrl['stations'][0]
    pass
  f.write('*\n$MODE;\n  def {};\n'.format(scan['mode']))
  for x in vex['MODE'][scan['mode']].getall('FREQ'):
    stations = set(x) - (set(x) - set(data.stations))
    if len(stations) > 0:
      f.write('    ref $FREQ = {}:{};\n'.format(x[0], ':'.join(stations)))
    if setup_station in x[1:]:
      freq = x[0]
  for x in vex['MODE'][scan['mode']].getall('BBC'):
    if setup_station in x[1:]:
      bbc = x[0]
  for x in vex['MODE'][scan['mode']].getall('IF'):
    if setup_station in x[1:]:
      if_ = x[0]
  stationstring = ':'.join(data.stations)
  f.write('    ref $BBC = {}:{};\n'.format(bbc, stationstring))
  f.write('    ref $IF = {}:{};\n'.format(if_, stationstring))
  nbits = get_nbits(vex, scan['mode'], data.stations)
  for n in nbits:
    f.write('    ref $TRACKS = tracks_{}bit:{};\n'.format(n, ":".join(nbits[n])))
  f.write('  enddef;\n')

  # STATION, SITE, ANTENNA, 'CLOCK', $TAPELOG_OBS, SOURCE
  blocks = [('STATION', 'all', ['$SITE', '$ANTENNA', '$CLOCK'])]
  blocks.append(('ANTENNA', 'all', ['antenna_diam', 'axis_type', 'axis_offset', 'pointing_sector']))
  blocks.append(('SITE', 'all', ['site_type', 'site_name', 'site_ID', 'site_position', \
                 'horizon_map_az', 'horizon_map_el', 'occupation_code', 'mk4_site_ID'])) 
  blocks.append(('TAPELOG_OBS', 'all', ['VSN']))
  blocks.append(('SOURCE', scan['source'], ['source_name', 'ra', 'dec', 'ref_coord_frame']))
  blocks.append(('CLOCK', 'all', ['clock_early']))
  blocks.append(('EOP', 'all', [])) # difx2mark4 leaves this empty
  blocks.append(('BBC', bbc, ['BBC_assign']))
  blocks.append(('IF', if_, ['if_def']))
  fmt = {'BBC_assign': '    BBC_assign = &{} : {} : &{};\n'}
  fmt['if_def'] = '    if_def = &{} : {} : {} : {} : {} : {} : {};\n'
  fmt['pointing_sector'] = '    pointing_sector = &{} : {} :  {} : {} : {} : {} : {};\n'
  for block, keys, items in blocks:
    keys = [k for k in vex[block].iterkeys()] if keys == 'all' else ensure_list(keys)
    f.write('*\n${};\n'.format(block))
    for k in keys:
      f.write('  def {};\n'.format(k))
      for item in items:
        pre = '    ref ' if item[0] == '$' else '    '
        for m in vex[block][k].getall(item.lstrip('$')):
          if item in fmt:
            if item == 'if_def' and len(m) < 6:
              m.append('0 MHz')
              pass
            if item == 'if_def' and len(m) < 7:
              m.append('0 MHz')
              pass
            f.write(fmt[item].format(*ensure_list(m)))
          else:
            f.write(pre + '{} = {};\n'.format(item, ' : '.join(ensure_list(m))))
      f.write('  enddef;\n')

  # TRACKS
  f.write('*\n$TRACKS;\n')
  for n in nbits:
    f.write('  def tracks_{}bit;\n    bits/sample = {};\n  enddef;\n'.format(n, n))
  
  # FREQ
  f.write('*\n$FREQ;\n')
  mode = scan['mode']
  for x in vex['MODE'][mode].getall('FREQ'):
    if setup_station in x[1:]:
      setupfreq = x[0]
      break
  for x in vex['MODE'][mode].getall('FREQ'):
    stations = list(set(x[1:]) - (set(x[1:]) - set(data.stations)))
    if len(stations) == 0:
      continue
    f.write('  def {};\n'.format(x[0]))
    for y in vex['FREQ'][setupfreq].getall('chan_def'):
      setupch = scan['allfreq'][setup_station][y[4]]
      sb = 0 if (setupch['sb'] == 'L') else 1
      smin = setupch['freq'] + (sb - 1) * setupch['bw']
      smax = setupch['freq'] + sb * setupch['bw']
      spol = setupch['pol']
      sid = setupch['freq_id']
      for chname in scan['allfreq'][stations[0]]:
        ch = scan['allfreq'][stations[0]][chname]
        sb = 0 if (ch['sb'] == 'L') else 1
        chmin = ch['freq'] + (sb - 1) * ch['bw']
        chmax = ch['freq'] + sb * ch['bw']
        chpol = ch['pol']
        # 1. Match case where band from setup station is contained in this subband.
        # 2. Match when this subband is contained in its entire in the setup station subband
        if (spol == chpol) and (((smin >= chmin) and (smax <= chmax)) or \
                                 (smin <= chmin) and (smax >= chmax)):
          f.write('    chan_def = {} : : {} : {} : {} : &{} : &{};\n'.format(sid, *y[1:]))
          break
    f.write('    sample_rate = {};\n  enddef;\n'.format(vex['FREQ'][setupfreq]['sample_rate']))

  # SCHED
  name = scan['name']
  sched = vex['SCHED'][name]
  f.write('*\n$SCHED;\n scan {};\n'.format(name))
  f.write('    start = {};\n'.format(sched['start']))
  f.write('    mode = {};\n'.format(sched['mode']))
  f.write('    source = {};\n'.format(sched['source']))
  for x in sched.getall('station'):
    if x[0] in data.stations:
      f.write('    station = {} : {} : {} : {} : {} : &{} : {};\n'.format(*x))
  f.write('  endscan;\n')
  
  # EVEX
  evex = "*\n$EVEX_REV;\n  rev = 1.0;\n$EVEX;\n  def 1234_std;\n" +\
         "    corr_exp#   = {};\n    ovex_file   = dummy;\n".format(BASENAME) +\
         "    cvex_file   = dummy;\n    svex_file   = dummy;\n" +\
         "    AP_length   = {:f} sec;\n    speedup_factor = 1.0;\n".format(ctrl["integr_time"]) +\
         "    ref $CORR_CONFIG = CDUM;\n    ref $SU_CONFIG  = SDUM;\n" +\
         "   enddef;\n"
  f.write(evex)

  # IVEX and LVEX
  ivex = "$IVEX_REV;\n  rev = 1.0;\n$CORR_INIT;\n  def INIT_DUMMY;\n" +\
         "    system_tempo = 1.00;\n    bocf_period =  8000000;\n" +\
         "    ref $PBS_INIT = PBS_DUMMY;\n  enddef;\n" +\
         "$PBS_INIT;\n  def PBS_DUMMY;\n  enddef;\n$LVEX_REV;\n" +\
         "  rev = 1.0;\n$LOG;\n  def log_dummy;\n  enddef;\n"
  f.write(ivex)
  return rootname

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

def write_t120s(data, scan, outfiles, t101_map, scalefactors, n):
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
      f.write(struct.pack("!2s6s2i", blname, ROOTID, t101['index'], n))
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
  rootname = create_root(vex, ctrl, scan, data, dirname)
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
  parser.add_argument('-r', "--rootid", help='Manually specify rootid', default=create_root_id(CREATIONDATE))
  args = parser.parse_args()
  vex = Vex(args.vexfile)
  ctrl = json.load(open(args.ctrlfile, 'r'))
  return vex, ctrl, args.rootid

def process_job(vex, ctrl, rootid, basename="1234"):
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
    scan = ctrl["scans"][0]
    source = vex['SCHED'][scan]['source']
    sources.append(vex['SOURCE'][source]['source_name'])
  except:
    pass
  out_tuple = urlparse(ctrl['output_file'])
  output_file = out_tuple.netloc if out_tuple.path == '' else out_tuple.path
  data = SFXCData(output_file, stations, sources)
  STATIONMAP = create_one_letter_mapping(vex)
  fix_vex_for_hopps(vex)

  scan = {'stop': datetime(1,1,1)}
  outfiles = {}
  n = 0
  while True:
    # Check if we need to move to next scan
    if data.current_time() >= scan["stop"]:
      # write end time and number of records for previous scan
      finalize_type1(outfiles, data, scan, n)
      scan, outfiles, t101_map, scale = initialise_next_scan(vex, exper, ctrl, data)
      n = 0

    write_t120s(data, scan, outfiles, t101_map, scale, n)
    n += 1

    if not data.next_integration():
      break
  # write end time and number of records for final scan
  finalize_type1(outfiles, data, scan, n)
#########
########################## MAIN #################################3
########
if __name__ == "__main__":
  vex, ctrl, rootid = parse_args()
  process_job(vex, ctrl, rootid)
