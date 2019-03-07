#!/usr/bin/env python
import struct
import numpy as np
from datetime import datetime, timedelta
from collections import namedtuple

class SFXCData:	    
  timeslice_header_size = 16
  uvw_header_size = 32
  stat_header_size = 24
  baseline_header_size = 8

  Channel = namedtuple('Channel', 'freqnr, sideband, pol')
  CrossChannel = namedtuple('CrossChannel', 'freqnr, sideband, pol1, pol2')
  Visibility = namedtuple('Visibility', 'vis, weight')
  
  def __init__(self, corfilename, stations=[], sources=[]):
    try:
      self.inputfile = open(corfilename, 'rb')
    except:
      print >> sys.stderr, "Error : Could not open " + corfilename
      sys.exit(1)
    
    self._integration_byte_pos = []
    self.vis = {}
    self.uvw = {}
    self.stats = {}
    self.stations = []
    self.channels = []
    self.current_int = -1
    self.current_slice = -1
    self.exp_stations = stations
    self.exp_sources = sources

    # Parse global header
    self._parse_global_header()

    # Read in the first integration
    self.next_integration()

  def current_time(self):
    """ Returns the datetime corresponding to the current integration """
    return self.start_time + self.current_slice * self.integration_time

  def prev_integration(self):
    """ Read in the the previous integration, returns True on success """
    if self.current_int <= 0:
      return False
    
    self.current_int -= 1
    inputfile.seek(self._integration_byte_pos[self.current_int])
    return self._read_integration()

  def next_integration(self):
    """ Read in the next integration, returns True on success """
    oldpos = self.inputfile.tell()
    if self._read_integration():
       self._integration_byte_pos.append(oldpos)
       self.current_int += 1
       return True
    return False

  def _read_integration(self):
    """ Read the next integration from disk, return True on success """
    timeslice_header_size = self.timeslice_header_size
    uvw_header_size = self.uvw_header_size
    stat_header_size = self.stat_header_size
    baseline_header_size = self.baseline_header_size

    inputfile = self.inputfile
    tsheader_buf = inputfile.read(timeslice_header_size)
    if len(tsheader_buf) != timeslice_header_size:
      return False
    
    timeslice_header = struct.unpack('4i', tsheader_buf)
    current_slice = timeslice_header[0]
    nchan = self.nchan
    exp_stations = self.exp_stations
    stations = set()
    channels = set()
    uvw = {}
    stats = {}
    vis = {}
    while timeslice_header[0] == current_slice:
      nbaseline = timeslice_header[1]
      nuvw = timeslice_header[2]
      nstatistics = timeslice_header[3]
      for i in range(nuvw):
        uvw_buf = inputfile.read(uvw_header_size)
        uvw_header = struct.unpack('i2h3d', uvw_buf)
        station = exp_stations[uvw_header[0]]
        stations.add(station)
        source = self.exp_sources[uvw_header[1]]
        uvw[station] = uvw_header[-3:]

      for i in range(nstatistics):
        stat_buf = inputfile.read(stat_header_size)
        stat_header = struct.unpack('4B5i', stat_buf)
        station = exp_stations[stat_header[0]]
        channel = self.Channel(*stat_header[1:4])
        channels.add(channel)
        try:
          stats[station][channel] = stat_header[-4:]
        except KeyError:
          stats[station] = {channel: stat_header[-4:]}

      for i in range(nbaseline):
        bl_buf = inputfile.read(baseline_header_size)
        weight, s1, s2, c, dummy = struct.unpack('i4B', bl_buf)
        bl = (exp_stations[s1], exp_stations[s2])
        pol1 = c&1
        pol2 = (c&2) >> 1
        sb = (c&4)>>2
        freq = c >> 3
        ch = self.CrossChannel(freq, sb, pol1, pol2)
        data = np.fromfile(inputfile, np.complex64, nchan + 1)
        try:
          vis[bl][ch] = self.Visibility(data, weight) 
        except KeyError:
          vis[bl] = {ch: self.Visibility(data, weight)}

      tsheader_buf = inputfile.read(timeslice_header_size)
      if len(tsheader_buf) < timeslice_header_size:
        break
      timeslice_header = struct.unpack('4i', tsheader_buf)

    if len(tsheader_buf) == timeslice_header_size:
      inputfile.seek(-timeslice_header_size, 1)

    self.stations = list(stations)
    self.channels = list(channels)
    self.current_slice = current_slice
    self.uvw = uvw
    self.stats = stats
    self.vis = vis
    self.source = source

    return True
  
  def _parse_global_header(self):
    """ Read global header from correlator file. """
    inputfile = self.inputfile
    gheader_size_buf = inputfile.read(4)
    self.global_header_size = struct.unpack('i', gheader_size_buf)[0]
    inputfile.seek(0)
    buf = inputfile.read(self.global_header_size)
    if self.global_header_size >= 92:
      htype = namedtuple('Header', 'size, exper, start_year, start_day, ' + \
                         'start_time, nchan, integr_time, ' + \
                         'output_format_version,sfxc_version, pol_type, ' + \
                         'sfxc_branch, jobnr, subjobnr, n_stations, ' \
                         'stations_offset, n_sources, sources_offset')
      h = htype._make(struct.unpack('i32s2h5ib15s2i4h', buf[:92]))
    else:
      htype = namedtuple('Header', 'size, exper, start_year, start_day, ' + \
                         'start_time, nchan, integr_time, ' + \
                         'output_format_version,sfxc_version, pol_type, ' + \
                         'sfxc_branch, jobnr, subjobnr')
      h = htype._make(struct.unpack('i32s2h5ib15s2i', buf[:84]))
    self.nchan = h.nchan

    self.integration_time = timedelta(0, h.integr_time / 1000000, h.integr_time % 1000000)
    self.start_time = datetime(h.start_year, 1, 1) + timedelta(h.start_day - 1, h.start_time, 0)
    # get stations and sources in a way that won't break when fields get to 
    # the header format
    if self.global_header_size >= 92:
      splitted = buf[92:].split('\0')
      self.exp_stations = splitted[:h.n_stations]
      self.exp_sources = splitted[h.n_stations:(h.n_stations + h.n_sources)]
    self.global_header = h
