#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys, struct
from numpy import *
from pylab import *
from time import sleep
from optparse import OptionParser
import pdb

timeslice_header_size = 16
uvw_header_size = 32
stat_header_size = 24
baseline_header_size = 8

plots_per_row = 2   # The number of pulse profiles per row in the plot window

def read_time_slice(inputfiles, visibilities, station_idx, ref_station, nchan):
  nbins = len(inputfiles)
  nif = visibilities.shape[2]
  numsb = visibilities.shape[3]
  npol = visibilities.shape[4]
  for bin in range(nbins):
    inputfile = inputfiles[bin]
    try:
      #get timeslice header
      tsheader_buf = inputfile.read(timeslice_header_size)
      timeslice_header = struct.unpack('4i', tsheader_buf)
    
      nuvw = timeslice_header[2]
      inputfile.seek(uvw_header_size * nuvw, 1)
    
      # Read the bit statistics
      nstatistics = timeslice_header[3]
      stat_buffer = inputfile.seek(stat_header_size * nstatistics, 1)
    
      nbaseline = timeslice_header[1]
      baseline_data_size = (nchan + 1) * 8 # data is complex floats
      baseline_buffer = inputfile.read(nbaseline * (baseline_header_size + baseline_data_size)) 
      fmt = str(2*(nchan+1)) + 'f'

      index = 0
      for b in range(nbaseline):
        bheader = struct.unpack('i4c', baseline_buffer[index:index + baseline_header_size])
        index += baseline_header_size
        station1 = struct.unpack('i', bheader[1]+'\x00\x00\x00')[0]
        station2 = struct.unpack('i', bheader[2]+'\x00\x00\x00')[0]
        byte = struct.unpack('i', bheader[3]+'\x00\x00\x00')[0]
        pol1 = byte&1
        pol2 = (byte>>1)&1
        sideband = (byte>>2)&1
        freq_nr = byte>>3
        #print 'pol1='+str(pol1) +  ', pol2='+str(pol2)+',sideband='+str(sideband)+',freq_nr='+str(freq_nr)+' ; s1='+str(station1)+', s2='+str(station2)
        J = complex(0,1)
        if (station1 == ref_station or station2 == ref_station) and (station1 != station2) and (pol1 == pol2):
          station = station_idx[station1] if (station2 == ref_station) else station_idx[station2]
          pol_idx = 0 if (npol == 1) else pol1
          sb_idx = 0 if (numsb == 1) else sideband

          buf = struct.unpack(fmt,baseline_buffer[index:index + baseline_data_size])
          #Skip over the first/last channel
          vreal = sum(buf[2:2*nchan:2])
          vim = sum(buf[3:2*nchan:2])
          if isnan(vreal)==False and isnan(vim)==False:
	    visibilities[station, bin, freq_nr, sb_idx, pol_idx] += vreal + J * vim
        index += baseline_data_size
    except struct.error:
      # triggered by EOF
      return False
      
  return True

def update_plots(profile, nint, prevprofile, prevnint, fig, plots):
  nstation = profile.shape[0]
  weight = nint * 1. / (nint + prevnint)
  newprofile = profile * weight + prevprofile * (1. - weight)
  interactive(True)
  #pdb.set_trace()
  for i in range(nstation):
    divider = max(newprofile[i,:].max(), 1)
    data = newprofile[i,:] / divider
    data = data - data.min()
    plots[i].set_ydata(data)
  fig.canvas.draw()
  fig.canvas.flush_events()

def get_station_list(vexfile):
  stations = []
  try:
    vex = open(vexfile, 'r')
  except:
    print 'Error : could not open vexfile : ' + vexfile
    sys.exit(1)
  line = vex.readline()
  while (line != '') and (line.lstrip().startswith('$STATION') == False):
    line = vex.readline()
  if(line == ''):
    print 'Couldn\'t find station block in vexfile'
    sys.exit(1)

  line = vex.readline()
  while (line != '') and (line.lstrip().startswith('$') == False):
    line = line.lstrip()
    if(line.startswith('def')):
      end = line.find(';')
      if(end < 0):
	print 'Error parsing vex file, station definition doesn\'t end with ;'
        print '  : ' + line
        exit(1)
      station = line[4:end].lstrip()
      stations.append(station)
      print 'found station #%d : %s'%(len(stations), station)
    line = vex.readline()
  vex.close()
  stations.sort()
  return stations

def initialize(base_file_name, nbins, station_list):
  inputfiles = []
  # open all input files and read global header
  for bin in range(1,nbins+1):
    filename = base_file_name + '.bin'+str(bin)
    try:
      inputfile = open(filename, 'rb')
    except:
      print "Error : Could not open " + filename
      sys.exit(1)
    gheader_size_buf = inputfile.read(4)
    global_header_size = struct.unpack('i', gheader_size_buf)[0]
    inputfile.seek(0)
    gheader_buf = inputfile.read(global_header_size)
    global_header = struct.unpack('i32s2h5i4c',gheader_buf[:64])
    nchan = global_header[5]
    inputfiles.append(inputfile)
  # determine parameters from first bin
  inputfile = inputfiles[0]
  #get timeslice header
  tsheader_buf = inputfile.read(timeslice_header_size)
  timeslice_header = struct.unpack('4i', tsheader_buf)
  integration_slice = timeslice_header[0]
  stations_found = zeros(len(station_list))
  nsubint = 0
  pol=0
  sb=0
  nif=0
  while (integration_slice == 0) and (len(tsheader_buf) == timeslice_header_size):
    # get the uvw buffer
    nuvw = timeslice_header[2]
    inputfile.seek(uvw_header_size * nuvw, 1)
    
    # Read the bit statistics
    nstatistics = timeslice_header[3]
    stat_buffer = inputfile.seek(stat_header_size * nstatistics, 1)
    
    nbaseline = timeslice_header[1]
    baseline_data_size = (nchan + 1) * 8 # data is complex floats
    baseline_buffer = inputfile.read(nbaseline * (baseline_header_size + baseline_data_size)) 

    index = 0
    for b in range(nbaseline):
      bheader = struct.unpack('i4c', baseline_buffer[index:index + baseline_header_size])
      index += baseline_header_size
      station1 = struct.unpack('i', bheader[1]+'\x00\x00\x00')[0]
      station2 = struct.unpack('i', bheader[2]+'\x00\x00\x00')[0]
      stations_found[station1] = 1
      stations_found[station2] = 1
      byte = struct.unpack('i', bheader[3]+'\x00\x00\x00')[0]
      pol |=  byte&3
      sb |= ((byte>>2)&1) + 1 
      nif = max((byte>>3) + 1, nif)
      #print 's1=%d, s2=%d, pol=%d, sb=%d, nif=%d, if_found=%d, sb_found=%d'%(station1, station2, pol, sb, nif, byte>>3, (byte>>2)&1)
      index += baseline_data_size
    tsheader_buf = inputfile.read(timeslice_header_size)
    if len(tsheader_buf) == timeslice_header_size:
        timeslice_header = struct.unpack('4i', tsheader_buf)
        integration_slice = timeslice_header[0]
    nsubint += 1
  
  inputfile.seek(global_header_size)
  stations_in_job = []
  for i in range(stations_found.size):
    if stations_found[i] == 1:
      stations_in_job.append(i)
  numsb = 2 if (sb == 3) else 1
  npol = 2 if (pol == 3) else 1
  return (inputfiles, nchan, nif, numsb, npol, stations_in_job, nsubint)

def create_plot_window(station_list, stations_in_job, ref_station, nbins):
  plots = []
  interactive(True)
  fig = figure()
  suptitle('Pulse profiles for station %s'%(station_list[ref_station]) , fontsize=12)
  nstation = len(stations_in_job) - 1
  nrows = int(ceil(nstation*1./ plots_per_row))
  f = 1
  empty_data = zeros([nbins])

  for i in range(len(stations_in_job)):
    station = stations_in_job[i]
    if  station != ref_station:
      plt = nrows * 100 + plots_per_row * 10 + f
      subplot(plt)
      p, = plot(empty_data)
      title('%s - %s'%(station_list[station], station_list[ref_station]))
      axis([0, nbins, 0, 1])
      draw()
      plots.append(p)
      f += 1
  fig.canvas.draw()
  fig.canvas.flush_events()
  return fig, plots

def get_options():
  parser = OptionParser('%prog <vex file> <base_file_name 1> ... <base_file_name N> <nbins> <ref_station>')
  (options, args) = parser.parse_args()

  N = len(args)
  if N < 4:
    parser.error('Invalid number of arguments')
  base_file_names = args[1:-2]
  return(args[0], base_file_names, int(args[-2]), args[-1])

######### MAIN CODE
(vexfile, base_file_names, nbins, ref_station_name) = get_options()
station_list = get_station_list(vexfile)
try:
  ref_station = station_list.index(ref_station_name)
except:
  print 'Ref station %s is not in vex file'%(ref_station_name)
  sys.exit(1)

scans = []
stations_in_job = set()
for base_file_name in base_file_names:
  scan = {}
  (inputfiles, nchan, nif, numsb, npol, stations_in_scan, nsubint) = initialize(base_file_name, nbins, station_list)
  scan["files"] = inputfiles
  scan["nchan"] = nchan
  scans.append(scan)
  stations_in_job = stations_in_job.union(stations_in_scan)
stations_in_job = list(stations_in_job)

station_idx = zeros(max(stations_in_job)+1,dtype=int)
#pdb.set_trace()
index = 0
for i in range(len(stations_in_job)):
  if stations_in_job[i] != ref_station:
    station_idx[stations_in_job[i]] = index
    index += 1
nstation = len(stations_in_job)

prev_profile = zeros([nstation - 1, nbins])
fig, plots = create_plot_window(station_list, stations_in_job, ref_station, nbins)
prev_nint = 0
update_freq = 1 # show the first integration immediately
for scan in scans:
  time_slice = 0
  visibilities = zeros([nstation - 1, nbins, nif, numsb, npol], dtype=complex128)
  while read_time_slice(scan["files"], visibilities, station_idx, ref_station, scan["nchan"]):
    time_slice += 1
    if time_slice % (nsubint * update_freq) == 0:
      profile = abs(visibilities).sum(2).sum(2).sum(2)
      nint = time_slice / nsubint
      update_plots(profile, nint, prev_profile, prev_nint, fig, plots)
      update_freq = 10 # Update plot once every 10 integrations
  # Scalar average profiles from different scans
  nint = time_slice / nsubint
  weight = nint * 1. / (nint + prev_nint)
  profile = abs(visibilities).sum(2).sum(2).sum(2)
  prev_profile = prev_profile * (1. - weight) + profile * weight
  prev_nint += nint

# Block until window is closed
print 'Done reading data.'
update_plots(zeros([nstation - 1, nbins]), 0, prev_profile, prev_nint, fig, plots)
show(block=True)
