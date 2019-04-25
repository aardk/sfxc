#!/usr/bin/env python
import re
import struct
import numpy as np
from datetime import datetime, timedelta
from collections import namedtuple

class experiment:
  bands = [['I',    100.0,    150.0], ['G',    150.0,    225.0], \
          ['P',    225.0,    390.0], ['L',    390.0,   1750.0], \
          ['S',   1750.0,   3900.0], ['C',   3900.0,   6200.0], \
          ['X',   6200.0,  10900.0], ['K',  10900.0,  36000.0], \
          ['Q',  36000.0,  46000.0], ['V',  46000.0,  56000.0], \
          ['W',  56000.0, 100000.0]]
  
  def __init__(self, vex):
    self.exper_name = vex['EXPER'][vex['GLOBAL']['EXPER']]['exper_name']
    modes = {}
    for mode in vex['MODE']:
      # Create per station frequency mapping for mode
      IFs = self._get_IFs(vex, mode)
      BBCs = self._get_BBCs(vex, mode)
      FREQs = self._get_FREQs(vex, mode, BBCs, IFs)
      modes[mode] = {'byname': {}, 'bykey': {}}
      for station in FREQs:
        modes[mode]['bykey'][station] = {} 
        modes[mode]['byname'][station] = {} 
        for key in FREQs[station]:
          f = FREQs[station][key]
          modes[mode]['byname'][station][f['name']] = f
          modes[mode]['bykey'][station][f['key']] = f
          f0, bw, sb, pol, freq_nr = (f['freq'], f['bw'], f['sb'], f['pol'], f['freq_nr'])
      self.modes = modes
        
  def _get_IFs(self, vex, mode):
    IFs = {}
    # Get all IFs
    for IFmode in vex['MODE'][mode].getall('IF'):
      IF = {}
      IFname = IFmode[0]
      for line in vex['IF'][IFname].getall('if_def'):
        IF[line[0]] = line[2]
      for station in IFmode[1:]:
        IFs[station] = IF
    return IFs

  def _get_BBCs(self, vex, mode):
    # Get all BBCs
    BBCs = {}
    for BBCmode in vex['MODE'][mode].getall('BBC'):
      BBC = {}
      BBCname = BBCmode[0]
      for line in vex['BBC'][BBCname].getall('BBC_assign'):
        BBC[line[0]] = line[2]
      for station in BBCmode[1:]:
        BBCs[station] = BBC
    return BBCs

  def _band(self, freq):
    for band in self.bands:
      if (freq >= band[1]) and (freq < band[2]):
        return band[0]
    # If band is not found return default band
    return 'B'

  def _get_FREQs(self, vex, mode, BBCs, IFs):
    # Now get all frequencies
    FREQs = {}
    for FREQmode in vex['MODE'][mode].getall('FREQ'):
      FREQ = {}
      FREQname = FREQmode[0]
      if len(FREQmode) == 1:
        continue
        
      freqs = []
      for line in vex['FREQ'][FREQname].getall('chan_def'):
        f0 = float(line[1].split()[0])
        bw = float(line[3].split()[0])
        ch = line[4]
        freqs.append((f0, ch))
        FREQ[ch] = {}
        FREQ[ch]['freq'] = f0
        FREQ[ch]['sb'] = line[2]
        FREQ[ch]['bw'] = bw
        FREQ[ch]['bbc'] = line[5]
        FREQ[ch]['name'] = ch
        S = FREQmode[1]
        IF = BBCs[S][line[5]]
        FREQ[ch]['pol'] = IFs[S][IF]
      freqs = sorted(freqs)
      fprev = freqs[0][0]
      i = 0
      for f in freqs:
        if f[0] != fprev:
          fprev = f[0]
          i += 1
        ch = f[1]
        sb = FREQ[ch]['sb']
        pol = FREQ[ch]['pol']
        FREQ[ch]['freq_nr'] = i
        freq_id = "{}{:02}{}{}".format(self._band(f[0]), i, sb, pol)
        FREQ[ch]['freq_id'] = freq_id
        # in the SFXC data channels are stored as freq_nr, sb, pol
        FREQ[ch]['key'] = (i, 0 if sb == 'L' else 1, 0 if pol == 'R' else 1)

      for station in FREQmode[1:]:
        FREQs[station] = FREQ
    return FREQs

  def _vex_time(self, timestring):
    year, yday, hour, minute, sec = [int(i) for i in re.split('[y|d|h|m|s]', timestring)[:-1]]
    secofday = 3600 * hour + 60 * minute + sec
    return datetime(year, 1, 1) + timedelta(yday-1, secofday) 

  def get_scan(self, vex, ctrl, data):
    if 'scans' in ctrl:
      scans = ctrl["scans"]
    else:
      scans = [scan for scan in vex['SCHED']]

    datatime = None
    for scan in scans:
      sourcefound = False
      for src in vex['SCHED'][scan].getall('source'):
        if vex['SOURCE'][src]['source_name'] == data.source:
          sourcefound = True
          break

      if sourcefound:
        tstart = self._vex_time(vex['SCHED'][scan]['start'])
        maxduration = 0
        data_stop = {}
        for line in vex['SCHED'][scan].getall('station'):
          if line[0] in data.stations:
            duration = int(line[2].split()[0])
            data_stop[line[0]] = duration
            maxduration = max(maxduration, duration)
        datatime = data.current_time()
        tstop = tstart + timedelta(0, maxduration)
        if (datatime >= tstart) and (datatime < tstop):
          try:
            setup_station = ctrl["setup_station"]
          except KeyError:
            setup_station = ctrl["stations"][0]
          print "scan = '{}', mode = '{}', start = '{}', stop = '{}', setup = '{}'".format(scan, vex['SCHED'][scan]['mode'], tstart, tstop, setup_station)
          result = {}
          mode = vex['SCHED'][scan]['mode']
          result['mode'] = mode
          result['start'] = tstart
          result['stop'] = tstop
          result['data_stop'] = data_stop
          result['name'] = scan
          result['source'] = data.source
          result['freq'] = self.modes[mode]['bykey'][setup_station]
          result['allfreq'] = self.modes[mode]['byname']
          return result

    print "Error, coulnd't find find scan for t =", datatime, ', source = ', data.source
    exit(1)
