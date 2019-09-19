#!/usr/bin/env python
import argparse
from vex import Vex 
from experiment import experiment
import stationmap
DEFAULT_EXPER_NUM = 16383

def ensure_list(x):
  if hasattr(x, '__iter__'):
    return list(x)
  return [x]

def create_global_ovex(vex, setup_station):
  exper = experiment(vex)
  exp = vex['GLOBAL']['EXPER']
  try:
    exper_num = vex['EXPER'][exp]['exper_num']
  except KeyError:
    exper_num = DEFAULT_EXPER_NUM

  outfile = open(str(exper_num)+'.ovex', 'w')
  scans = [scan for scan in vex['SCHED']]
  stations = [station for station in vex['STATION']]
  write_ovex(vex, setup_station, scans, stations, exper, outfile)

def create_rootfile(vex, ctrl, scan, stations, dirname, rootid):
  exper = experiment(vex)
  source = vex['SCHED'][scan]["source"]
  rootname = dirname + '/' + source + '.' + rootid
  outfile = open(rootname, 'w')
  try:
    setup_station = ctrl['setup_station']
  except:
    setup_station = ctrl['stations'][0]

  scans = ensure_list(scan) 
  write_ovex(vex, setup_station, scans, stations, exper, outfile)
  write_evex(outfile, ctrl["integr_time"], dirname)
  write_ivex(outfile)
  return rootname

def write_ovex(vex, setup_station, scans, stations, exper, outfile):
  # $Head and $GLOBAL
  outfile.write('$OVEX_REV;\nrev = 1.5;\n$GLOBAL;\n')
  for i in ('EOP', 'EXPER'):
    outfile.write('    ref ${} = {};\n'.format(i, vex['GLOBAL'][i]))
  
  # $EXPER
  name = vex['GLOBAL']['EXPER']
  outfile.write('*\n$EXPER;\n  def {};\n'.format(name))
  for i in ('exper_num', 'exper_name', 'exper_description', 'PI_name'):
    outfile.write('    {} = {};\n'.format(i, vex['EXPER'][name][i]))
  # Target correlator has to difx, otherwise HOPS assumes it is lagdata
  outfile.write('    target_correlator = difx;\n  enddef;\n')

  # $MODE
  modes = set(vex['SCHED'][scan]['mode'] for scan in scans)
  bbcs = set()
  ifs = set()
  for mode in modes:
    outfile.write('*\n$MODE;\n  def {};\n'.format(mode))
    mode_stations = []
    for x in vex['MODE'][mode].getall('FREQ'):
      st = set(x) - (set(x) - set(stations))
      mode_stations += st
      if len(mode_stations) > 0:
        outfile.write('    ref $FREQ = {}:{};\n'.format(x[0], ':'.join(st)))
    for x in vex['MODE'][mode].getall('BBC'):
      if setup_station in x[1:]:
        bbc = x[0]
    for x in vex['MODE'][mode].getall('IF'):
      if setup_station in x[1:]:
        if_ = x[0]
    stationstring = ':'.join(mode_stations)
    outfile.write('    ref $BBC = {}:{};\n'.format(bbc, stationstring))
    bbcs.add(bbc)
    outfile.write('    ref $IF = {}:{};\n'.format(if_, stationstring))
    ifs.add(if_)
    nbits = get_nbits(vex, mode, mode_stations)
    for n in nbits:
      outfile.write('    ref $TRACKS = tracks_{}bit:{};\n'.format(n, ":".join(nbits[n])))
    outfile.write('  enddef;\n')

  # STATION, SITE, ANTENNA, 'CLOCK', $TAPELOG_OBS, SOURCE
  blocks = [('STATION', 'all', ['$SITE', '$ANTENNA', '$CLOCK'])]
  blocks.append(('ANTENNA', 'all', ['antenna_diam', 'axis_type', 'axis_offset', 'pointing_sector']))
  blocks.append(('SITE', 'all', ['site_type', 'site_name', 'site_ID', 'site_position', \
                 'horizon_map_az', 'horizon_map_el', 'occupation_code', 'mk4_site_ID'])) 
  blocks.append(('TAPELOG_OBS', 'all', ['VSN']))
  sources = set(vex['SCHED'][scan]['source'] for scan in scans)
  blocks.append(('SOURCE', sources, ['source_name', 'ra', 'dec', 'ref_coord_frame']))
  blocks.append(('CLOCK', 'all', ['clock_early']))
  blocks.append(('EOP', 'all', [])) # difx2mark4 leaves this empty
  blocks.append(('BBC', bbcs, ['BBC_assign']))
  blocks.append(('IF', ifs, ['if_def']))
  fmt = {'BBC_assign': '    BBC_assign = &{} : {} : &{};\n'}
  fmt['if_def'] = '    if_def = &{} : {} : {} : {} : {} : {} : {};\n'
  fmt['pointing_sector'] = '    pointing_sector = &{} : {} :  {} : {} : {} : {} : {};\n'
  for block, keys, items in blocks:
    keys = [k for k in vex[block].iterkeys()] if keys == 'all' else ensure_list(keys)
    outfile.write('*\n${};\n'.format(block))
    for k in keys:
      outfile.write('  def {};\n'.format(k))
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
            outfile.write(fmt[item].format(*ensure_list(m)))
          else:
            outfile.write(pre + '{} = {};\n'.format(item, ' : '.join(ensure_list(m))))
      outfile.write('  enddef;\n')

  # TRACKS
  outfile.write('*\n$TRACKS;\n')
  for n in nbits:
    outfile.write('  def tracks_{}bit;\n    bits/sample = {};\n  enddef;\n'.format(n, n))
  
  # FREQ
  outfile.write('*\n$FREQ;\n')
  for mode in modes:
    for x in vex['MODE'][mode].getall('FREQ'):
      if setup_station in x[1:]:
        setupfreq = x[0]
        break
    for x in vex['MODE'][mode].getall('FREQ'):
      scan_stations = list(set(x[1:]) & set(stations))
      if len(scan_stations) == 0:
        continue
      outfile.write('  def {};\n'.format(x[0]))
      for y in vex['FREQ'][setupfreq].getall('chan_def'):
        setupch = exper.get_freq(mode, setup_station)[y[4]]
        sb = 0 if (setupch['sb'] == 'L') else 1
        smin = setupch['freq'] + (sb - 1) * setupch['bw']
        smax = setupch['freq'] + sb * setupch['bw']
        spol = setupch['pol']
        sid = setupch['freq_id']
        chmap = exper.get_freq(mode, setup_station)
        for chname in chmap:
          ch = chmap[chname]
          sb = 0 if (ch['sb'] == 'L') else 1
          chmin = ch['freq'] + (sb - 1) * ch['bw']
          chmax = ch['freq'] + sb * ch['bw']
          chpol = ch['pol']
          # 1. Match case where band from setup station is contained in this subband.
          # 2. Match when this subband is contained in its entire in the setup station subband
          if (spol == chpol) and (((smin >= chmin) and (smax <= chmax)) or \
                                   (smin <= chmin) and (smax >= chmax)):
            outfile.write('    chan_def = {} : : {} : {} : {} : &{} : &{};\n'.format(sid, *y[1:]))
            break
      outfile.write('    sample_rate = {};\n  enddef;\n'.format(vex['FREQ'][setupfreq]['sample_rate']))

  # SCHED
  outfile.write('*\n$SCHED;\n')
  for scan in scans:
    sched = vex['SCHED'][scan]
    scan_stations = set(stations) & set(x[0] for x in sched.getall('station'))
    if len(scan_stations) < 2:
      continue
    outfile.write(' scan {};\n'.format(scan))
    outfile.write('    start = {};\n'.format(sched['start']))
    outfile.write('    mode = {};\n'.format(sched['mode']))
    outfile.write('    source = {};\n'.format(sched['source']))
    for x in sched.getall('station'):
      if x[0] in stations:
        outfile.write('    station = {} : {} : {} : {} : {} : &{} : {};\n'.format(*x))
    outfile.write('  endscan;\n')
  
def write_evex(outfile, integr_time, dirname):
  basename = dirname.split('/')[0] # defaults to 1234
  # EVEX
  evex = "*\n$EVEX_REV;\n  rev = 1.0;\n$EVEX;\n  def 1234_std;\n" +\
         "    corr_exp#   = {};\n    ovex_file   = dummy;\n".format(basename) +\
         "    cvex_file   = dummy;\n    svex_file   = dummy;\n" +\
         "    AP_length   = {:f} sec;\n    speedup_factor = 1.0;\n".format(integr_time) +\
         "    ref $CORR_CONFIG = CDUM;\n    ref $SU_CONFIG  = SDUM;\n" +\
         "   enddef;\n"
  outfile.write(evex)

def write_ivex(outfile):
  # IVEX and LVEX
  ivex = "$IVEX_REV;\n  rev = 1.0;\n$CORR_INIT;\n  def INIT_DUMMY;\n" +\
         "    system_tempo = 1.00;\n    bocf_period =  8000000;\n" +\
         "    ref $PBS_INIT = PBS_DUMMY;\n  enddef;\n" +\
         "$PBS_INIT;\n  def PBS_DUMMY;\n  enddef;\n$LVEX_REV;\n" +\
         "  rev = 1.0;\n$LOG;\n  def log_dummy;\n  enddef;\n"
  outfile.write(ivex)

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
  one_letter_codes = stationmap.one_letter_codes

  # Add mk4_site_ID to $SITE
  for key in vex['STATION']:
    site = vex['STATION'][key]['SITE']
    st = vex['SITE'][site]['site_ID']
    vex['SITE'][site]['mk4_site_ID'] = one_letter_codes[st]

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

def parse_args():
  parser = argparse.ArgumentParser(description='Converts a vex file to an ovex file.')
  parser.add_argument("vexfile", help='vex file')
  parser.add_argument("setup_station", help='Setup station for the experiment')
  parser.add_argument('-c', "--code-file", help='JSON file containing one letter station codes')
  args = parser.parse_args()
  vex = Vex(args.vexfile)
  setup_station = args.setup_station
  if setup_station not in vex['STATION']:
    parser.error('Setup station is not in vex file')
  return vex, setup_station, args.code_file

#########
########################## MAIN #################################3
########
if __name__ == "__main__":
  vex, setup_station, codefile = parse_args()
  stationmap.create_one_letter_mapping(vex, codefile)
  fix_vex_for_hopps(vex)
  create_global_ovex(vex, setup_station)
