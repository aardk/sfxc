#!/usr/bin/env python
from pylab import *
import sys
import struct
from collections import namedtuple

class hopsdate:
  def __init__(self, year, day, hour, minute, second):
    self.year = year
    self.day = day
    self.hour = hour
    self.minute = minute
    self.second = second
  
  def __str__(self):
    return '{0:4d}y{1:03d}d{2:02d}h{3:02d}m{4:09.6f}s'.format(self.year, self.day, self.hour, self.minute, self.second)

def read_t000(f):
  recid, version, unused, date, name = struct.unpack("3s2s3s16s40s", f.read(64))
  print 'date = "'+date+'", name = "' + name + '"'
  print 'Found t000, recid = {}, version = {}, unused = {}, date = {}, name = {}'.format(recid, version, unused, date, name)

def read_t300(f):
  # Type            ascii           3       300
  # Version         ascii           2       0-99 
  # Unused          ascii           2       Spaces
  # SU_number       i*1             1       Station unit, filled by suman
  # Id              ascii           1       1-char vex letter code
  # Intl_id         ascii           2       2-char international station code
  # Name            ascii           32      Full station name
  # Unused          ascii           1       Padding for alignment
  # Model_date      date            12      Start time for 1st spline
  # Model interval  r*4             4       Spline interval in seconds
  # Nsplines        i*2             2       Number of splines in scan
  recid, version, unused, su = struct.unpack("3s2s2sB", f.read(8))
  ID, Intl_id, name, unused2 = struct.unpack("s2s32ss", f.read(36))
  model_date = hopsdate(*struct.unpack("!hhhhf", f.read(12)))
  model_interval, nsplines, unused3 = struct.unpack("!fhh", f.read(8))
  print 'Found t300, Type = {}, version = {}, unused = {}, su_nr = {}, Id = {}, Intl_id = {}'.format(recid, version, unused, su, ID, Intl_id)
  print '  name = {}, unused = {}, model_date = {}, model_interval = {}, Nsplines = {}'.format(name, unused2, model_date, model_interval, nsplines)
  return nsplines

def read_t301(f, recid):
  # Type            ascii           3       301
  # Version         ascii           2       0-99 
  # Unused          ascii           3       Spaces
  # Interval        i*2             2       Sequential model interval number
  # Chan_id         ascii           32      Frequency channel identifier
  # Unused          ascii           6       Padding for alignment
  # Delay_spline    r*8 x 6         48      Delay spline coefficients
  version, unused, interval, chan_id, unused2 = struct.unpack("!2s3sh32s6s", f.read(45))
  print 'Found t301, recid = {}, version = {}, unused = {}, interval = {}, chan_id = {}, unused = {}'.format(recid, version, unused, interval, chan_id, unused2)
  spline = struct.unpack("!6d", f.read(6 * 8))
  print 'delay spline: {}'.format(spline)

def read_t302(f, recid):
  # Type            ascii           3       302
  # Version         ascii           2       0-99 
  # Unused          ascii           3       Spaces
  # Interval        i*2             2       Sequential model interval number
  # Chan_id         ascii           32      Frequency channel identifier
  # Unused          ascii           6       Padding for alignment
  # Phase_spline    r*8 x 6         48      Phase spline coefficients
  version, unused, interval, chan_id, unused2 = struct.unpack("!2s3sh32s6s", f.read(45))
  print 'Found t302, recid = {}, version = {}, unused = {}, interval = {}, chan_id = {}, unused = {}'.format(recid, version, unused, interval, chan_id, unused2)
  spline = struct.unpack("!6d", f.read(6 * 8))
  print 'phase spline: {}'.format(spline)

def read_t303(f, recid):
  # char        record_id[3];           /* Standard 3-digit id */
  # char        version_no[2];          /* Standard 2-digit version # */
  # char        unused1[3];             /* Reserved space */
  # short   interval;                   /* Sequential model interval number */
  # char    chan_id[32];                /* Frequency channel identifier */
  # char    unused2[6];                 /* Padding */
  # double  azimuth[6];                 // Azimuth (deg) coefficients
  # double  elevation[6];               // Elevation (deg) coefficients
  # double parallactic_angle[6];        // Par. angle (deg CCW el line from RA line)
  # double u[6];                        // Baseline projections toward source (m)
  # double v[6];
  # double w[6];
  version, unused, interval, chan_id, unused2 = struct.unpack("!2s3sh32s6s", f.read(45))
  print 'Found t303, recid = {}, version = {}, unused = {}, interval = {}, chan_id = {}, unused = {}'.format(recid, version, unused, interval, chan_id, unused2)
  rec = struct.unpack("!36d", f.read(36 * 8))
  az, el, par, u, v, w = [rec[6*i:6*(i+1)] for i in range(6)]
  print 'Azimuth spline: {}'.format(az)
  print 'Elevation spline: {}'.format(el)
  print 'Parallactic angle spline: {}'.format(par)
  print 'u spline: {}'.format(u)
  print 'v spline: {}'.format(v)
  print 'w spline: {}'.format(w)


def read_t309(f, recid):
  # char        record_id[3];           // Standard 3-digit id
  # char        version_no[2];          // Standard 2-digit version #
  # char        unused1[3];
  # int         su;                     // SU
  # int         ntones;                 // number of tones [0..64]
  # double      rot;                    // ROT at start of AP
  # double      acc_period;             // in secs
  # struct ch1_tag
  #  {
  #  char      chan_name[8];
  #  double    freq;                   // tone frequency in Hz
  #  U32       acc[64][2];             // accumulators for 64 freqs x 2 quads (C..S)
  #  } chan[64];
  version, unused, su, ntones, rot, acc = struct.unpack("!2s3s2i2d", f.read(29))
  print 'Found t309, recid = {}, version = {}, unused = {}, su = {}, ntones = {}, rot = {}, acc_period = {}'.format(recid, version, unused, su, ntones, rot, acc)
  for i in range(64):
    name, freq = struct.unpack('!8sd', f.read(16))
    print ' ch_name = {}, freq = {}, accumulators:'.format(name, freq)
    for j in range(4):
      acc = struct.unpack('!32i', f.read(32*4))
      for k in range(16):
        print '({}, {}) '.format(acc[2*k], acc[2*k+1]),
      print

if len(sys.argv) != 2:
  print 'Usage: {} <type 3 file>'.format(sys.argv[0])
  exit(1)

f = open(sys.argv[1], 'r')
read_t000(f)
nsplines = read_t300(f)
recid = f.read(3)
while len(recid) == 3:
  if recid == '301': 
    read_t301(f, recid)
  elif recid == '302': 
    read_t302(f, recid)
  elif recid == '303': 
    read_t303(f, recid)
  elif recid == '309': 
    read_t309(f, recid)
  else:
    print 'Unknown record id:"{}"'.format(recid)
    break
  recid = f.read(3)
