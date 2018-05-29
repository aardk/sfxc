#!/usr/bin/env python
from math import sqrt
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

def read_t100(f):
  recid, version, unused = struct.unpack("3s2s3s", f.read(8))
  procdate = hopsdate(*struct.unpack("!hhhhf", f.read(12)))
  baseline, rootname, qcode, padding, percdone = struct.unpack('!2s34s2s6sf', f.read(48))
  print 'Found t100, recid = {}, version = {}, unused = {}, procdate = {}'.format(recid, version, unused, procdate)
  print '  baseline = {}, rootname = {}, qcode = {}, paddinf = {}, percentage done = {}'.format(baseline, rootname, qcode, padding, percdone)
  start = hopsdate(*struct.unpack("!hhhhf", f.read(12)))
  stop = hopsdate(*struct.unpack("!hhhhf", f.read(12)))
  ndrec, nindex, nlags, nblocks = struct.unpack('!2i2h', f.read(12))
  print '  start = {}, stop = {}, ndrec = {}, nindex = {}, nlags = {}, nblocks = {}'.format(start, stop, ndrec, nindex, nlags, nblocks)
  return ndrec, nindex, nlags

def read_t101(f):
  recid, version, status = struct.unpack("3s2ss", f.read(6))
  print 'Found t101, recid = {}, version = {}, status = {}'.format(recid, version, status)
  nblocks, index, primary = struct.unpack("!3h", f.read(6))  
  print '  nblocks = {}, index = {}, primary = {}'.format(nblocks, index, primary)
  ref_chan_id, rem_chan_id, board, slot = struct.unpack("!8s8s2h", f.read(20))
  print '  ref_chan_id = {}, rem_chan_id = {}, board = {}, slot = {}'.format(ref_chan_id, rem_chan_id, board, slot)
  ref_chan, rem_chan, post, blocks = struct.unpack("!2h2i", f.read(12))
  print '  ref_chan = {}, rem_chan = {}, post mortem= {}, blocks = {}'.format(ref_chan, rem_chan, post, blocks)

def read_t120(f):
  recid, version, corrtype, nlags, baseline, rootcode = struct.unpack("!3s2sbh2s6s", f.read(16))
  print 'Found t120, recid = {}, version = {}, corrtype = {}, nlags = {}, baseline = {}, rootcode = {}'.format(recid, version, corrtype, nlags, baseline, rootcode)
  index, ap = struct.unpack('!2i', f.read(8))
  print '  index = {}, ap = {},'.format(index, ap),
  if corrtype == 5:
    flag_wgt = struct.unpack('!f', f.read(4))
    print 'weight = {:.6},'.format(flag_wgt[0]),
  else:
    flag_wgt = struct.unpack('!i', f.read(4))
    print 'flag = {}'.format(flag_wgt[0]),
  status, fr_delay, delay_rate = struct.unpack('!3i', f.read(12))
  data = struct.unpack('!{}f'.format(2*nlags), f.read(2 * 4 * nlags))
  amp = sqrt(sum([x**2 for x in data]) / nlags)
  print 'status = {}, fr_delay = {}, delay_rate = {}, amp = {} '.format(status, fr_delay, delay_rate, amp)
#Type            ascii           3       120
#Version         ascii           2       0-99 
#corrtype        ascii           1       Binary 1-4, lagdata struct type
#
#nlags           i*2             2       Number of lags present
#Baseline        ascii           2       Standard baseline id
#rootcode        ascii           6       Standard root suffix
#
#Index           i*4             4       Index into type 101 record
#AP              i*4             4       Accumulation period number
#Flag            i*4             4       Up to 32 correlation flags
#Status          i*4             4       Up to 32 status bits
#fr_delay        i*4             4       Mid-AP fractional delay (bits * 2^32)
#delay_rate      i*4             4       Mid-AP delay rate (bits/sysclk * 2^32)
#lagdata         array        variable   Correlation counts

if len(sys.argv) != 2:
  print 'Usage: {} <type 1 file>'.format(sys.argv[0])
  exit(1)

f = open(sys.argv[1], 'r')
read_t000(f)
ndrec, nindex, nlags = read_t100(f)
for i in range(nindex):
  read_t101(f)

for i in range(ndrec):
  read_t120(f)
