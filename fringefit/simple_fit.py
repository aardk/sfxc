#!/usr/bin/env python
from pylab import *
from scipy.optimize import minimize
import numpy.polynomial as polynomial
import sys, struct, pdb
from datetime import datetime, timedelta
import parameters, vex_time
import vex as Vex
import json
import signal
from sfxcdata_utils import *
from optparse import OptionParser
import time

def p_good(V, n, max_int, dx):
  if(V > 13): # for V>13 the probility of error at the machine precision
    return 1;
  Z = arange(0,max_int,dx)
  return trapz(Z*exp(-(Z**2+V**2)/2)*i0(Z*V)*(1-exp(-(Z**2)/2))**(n-1),dx=dx)

def goto_integration(corfile, start_pos, timeout):
  corfile.seek(0, 2)
  pos = corfile.tell()
  if (timeout > 0):
    oldpos = pos
    oldtime = datetime.utcnow()
    while pos < start_pos:
      time.sleep(1)
      corfile.seek(0, 2)
      pos = corfile.tell()
      newtime = datetime.utcnow()
      if pos != oldpos:
        oldpos = pos
        oldtime = newtime
      elif (newtime-oldtime).seconds > timeout:
        break

  if pos < start_pos:
    return False

  corfile.seek(start_pos)
  return True

def read_integrations(inputfile, data, int_read, param, scan_param, n_integrations, ref_station, timeout):
  global_header_size = param.global_header_size
  timeslice_header_size = parameters.timeslice_header_size
  uvw_header_size = parameters.uvw_header_size
  stat_header_size = parameters.stat_header_size
  baseline_header_size = parameters.baseline_header_size

  stations_in_job = scan_param['stations']
  n_stations = len(stations_in_job)
  channels = scan_param['channels']
  n_channels = len(channels)
  nsubint = scan_param['nsubint']
  nchan = param.nchan

  n_baseline = n_stations*(n_stations-1)/2
  #print stations_in_job
  for i in range(n_integrations):
    #print `i`+'/'+`n_integrations`
    for j in range(nsubint):
      #print `j`+'/'+`nsubint`
      tsheader_buf = read_data(inputfile, timeslice_header_size, timeout)
      timeslice_header = struct.unpack('4i', tsheader_buf)
      # skip over headers
      nbaseline = timeslice_header[1]
      nuvw = timeslice_header[2]
      nstatistics = timeslice_header[3]
      inputfile.seek(uvw_header_size * nuvw + stat_header_size * nstatistics, 1)

      baseline_data_size = (nchan + 1) * 8 # data is complex floats
      baseline_buffer = read_data(inputfile, nbaseline * (baseline_header_size + baseline_data_size), timeout)
      index = 0
      baselines = []
      for b in range(nbaseline):
        #print `b`+'/'+`nbaseline`
        bheader = struct.unpack('i4b', baseline_buffer[index:index + baseline_header_size])
        index += baseline_header_size
        station1 = stations_in_job.index(bheader[1])
        station2 = stations_in_job.index(bheader[2])
        #print 's1='+`bheader[1]`+', s2='+`bheader[2]`, 'station1',station1,', station2',station2, ', ref_station = ', ref_station
        pol1 = bheader[3]&1
        pol2 = (bheader[3]&2) >> 1
        #print 'pol1 = ' + `pol1` + ', pol2= ' + `pol2`
        if ((station1 == ref_station) ^ (station2 == ref_station)) and (pol1 == pol2):
          pol = bheader[3]&1
          sb = (bheader[3]&4)>>2
          freq = bheader[3] >> 3
          #print 'pol = ' + `pol` + ', sb = ' + `sb` + ', freq = ' + `freq`
          chan = channels.index([freq, sb, pol])
          size_str = `2*nchan+2`+'f'
          #pdb.set_trace()
          values = array(struct.unpack(size_str, baseline_buffer[index:index + baseline_data_size]), dtype='f8')
          station = station1 if station2 == ref_station else station2
          # Flip phase if needed
          if station1 == ref_station:
            data[chan, station, i, :] =  values[0:2*nchan+2:2] + complex64(1j) * values[1:2*nchan+2:2]
          else:
            data[chan, station, i, :] =  values[0:2*nchan+2:2] - complex64(1j) * values[1:2*nchan+2:2]
        index += baseline_data_size
    int_read[0] += 1

def lag_offsets(data, n_station, delays, rates):
  #pdb.set_trace()
  for chan in range(data.shape[0]):
    for station in range(n_station):
      # First determine peak in lag space
      lag = irfft2(data[chan, station, :, :])
      rl = abs(lag)
      max_idx = rl.argmax()
      idx = max_idx % rl.shape[1]
      idy = max_idx / rl.shape[1]
      x = idx if idx < rl.shape[1]/2 else idx - rl.shape[1]
      y = idy if idy < rl.shape[0]/2 else idy - rl.shape[0]
      delays[chan,station] = x
      rates[chan,station] = y

def apply_model(data, station1, station2, delays, rates, param, scan_param):
  channels = scan_param['channels']
  vex_freqs = param.freqs
  n_station = len(scan_param['stations'])
  n_baseline = n_station * (n_station-1) / 2

  N = data.shape[2] - 1
  T = data.shape[1]
  f = arange(0, N+1) /(2.*N)
  t = arange(0, T) / float(T)
  bldata = zeros(data.shape, dtype='c8')
  for ch in range(len(channels)):
    rate1 = rates[ch, station1]
    rate2 = rates[ch, station2]
    delay1 = delays[ch, station1]
    delay2 = delays[ch, station2]
    delay = (delay1 - delay2)
    rate = (rate1 - rate2)
    if (ch == 0) and (station1==-1):
      pdb.set_trace()
    bldata[ch] = data[ch] * exp(-2j*pi*f*delay)
    bldata[ch] = (bldata[ch].T * exp(-2j*pi*t*rate)).T
  return bldata

# FIXME: Create single cost function for delay and rate
def cost_function_delay(delay, vis):
    N = (vis.size - 1) * 2
    v = exp(-2j * pi * delay * arange(N/2 + 1) / N) * vis
    vsum = v[N/10:9*N/10].sum()
    return 1. / (abs(vsum) + 1e-8)

def cost_function_rate(rate, vis):
    N = vis.size
    v = exp(-2j * pi * rate * arange(N) / N) * vis
    vsum = v.sum()
    return 1. / (abs(vsum) + 1e-8)

def refine_offsets(data, station, delays, rates):
  nchan = data.shape[0]
  nfreq = data.shape[2]^1
  # Only use the inner 80% of the band
  fmin = nfreq / 10
  fmax = nfreq+1 - nfreq/10
  f = arange(fmin, fmax) * pi / nfreq
  #pdb.set_trace()
  for chan in range(nchan):
    # Find the peak lag, for coarse delay compensation
    fringe = irfft2(data[chan, :, :])
    max_idx = abs(fringe).argmax()
    idx = max_idx % fringe.shape[1]
    idy = max_idx / fringe.shape[1]
    vis = rfft(fringe[idy,:])
    result_delay = minimize(cost_function_delay, 0., args=(vis,), bounds=[(-2, 2)], method='L-BFGS-B')
    delay = -result_delay.x[0]
    if station == -4:
        pdb.set_trace()
    vis_rate = sum(data[chan, :, fmin:fmax]*exp(-1j * f * delay), axis=1)
    result_rate = minimize(cost_function_rate, 0., args=(vis_rate,), bounds=[(-2, 2)], method='L-BFGS-B')
    rate = -result_rate.x[0]
    # Store results
    delays[chan,station] = delay
    rates[chan,station] = rate

def compute_snr(data, snr, station, sample_rate, integration_time):
  nchan = data.shape[0]
  nint = data.shape[1]
  nlag = (data.shape[-1] - 1) * 2
  # Noise estimate based on perfect sampler stats and 2bit sampling
  # NB: This requires normalized visibilities
  noise = sqrt(nint) / (0.881 * sqrt(sample_rate * 1e6 * integration_time))
  for chan in range(nchan):
    # Post fringe fitting
    fringe = abs(irfft(data[chan, :, :]).sum(axis=0))
    snr[chan, station] = fringe.max() / noise

def get_options():
  parser = OptionParser('%prog [options] <vex file> <reference station>')
  parser.add_option('-c', '--ctrl', dest='controlfile',
                    help="Use json control file")
  parser.add_option('-f', '--file',
                    dest='corfile',
                    help='Correlator output file')
  parser.add_option('-b', '--begin-time', dest='begin_time',
                    help='Start time of clock search [Default : first integration]')
  parser.add_option('-e', '--end-time', dest='end_time',
                    help='End time of clock search [Default : last integration]')
  parser.add_option('-t', '--timeout', dest='timeout', type='int', default = 0,
                    help='Determines after how many seconds of inactivity we assume that the correlator job ended. [Default : 0 seconds]')
  parser.add_option('-s', '--setup-station', help='Setup station (if not specified in control file)')
  (options, args) = parser.parse_args()
  if len(args) != 2:
    parser.error('invalid number of arguments')
  if not (options.controlfile != None) ^ (options.corfile !=None):
    parser.error('either a control file or a correlator output file needs to be specified.')
  vex_name = args[0]
  ref_station = args[1]
  ctrl = None
  setup_station = options.setup_station
  if options.controlfile != None:
    try:
      ctrl = json.load(open(options.controlfile, "r"))
      if setup_station == None:
          try:
            setup_station = ctrl["setup_station"]
          except KeyError:
            setup_station = ctrl["stations"][0]
      # TODO : Add file extensions for pulsar binning / multiple phase center
      if ctrl["output_file"][:7] != 'file://':
        raise StandardError('Correlator output_file should start with file://')
      corfile = ctrl["output_file"][7:]
      try:
        if ctrl["pulsar_binning"]:
          corfile += ".bin1"
      except KeyError:
        # No pulsar_binning keyword in control file
        pass
    except StandardError, err:
      print >> sys.stderr, "Error loading control file : " + str(err)
      sys.exit(1)
  else:
    corfile = options.corfile
    if setup_station == None:
        setup_station = ref_station

  try:
    vex = Vex.parse(open(vex_name, "r").read())
  except StandardError, err:
    print >> sys.stderr, "Error loading vex file : " + str(err)
    sys.exit(1)
  try:
    inputfile = open(corfile, 'rb')
  except:
    print >> sys.stderr, "Error : Could not open " + corfile
    sys.exit(1)
  if options.timeout > 0:
    # We are running concurrently with a correlator job
    if (options.end_time == None) and (options.controlfile == None):
        parser.error("When running concurrently with a correlator job either the end_time option should\
         be set or the control file should be used to drive this program.")
  param = parameters.parameters(vex, corfile, setup_station, options.timeout)
  return (vex, ctrl, inputfile, ref_station, options, param)

def write_clocks(vex, param, scan_param, delays, rates, snr, ref_station, begin_time, n_integrations):
  vex_stations = [s for s in vex['STATION']]
  vex_stations.sort()
  # First compute channel weights
  W = zeros(snr.shape)
  for c in range(n_channels):
    for b in range(n_stations):
      P = p_good(snr[c,b], param.nchan, 20, 0.01)
      W[c, b] = P**4/(P**4 + 16.*(1-P)**4)
  N = snr.shape[1]
  tot_weights = sum(W,axis=0)
  tot_weights[ref_station] = 1
  weights = W / (tot_weights + 1e-6)
  weights[:,ref_station] = 1
  stations = [vex_stations[i] for i in scan_param['stations']]
  channels = [param.channel_names[param.vex_channels.index(c)] for c in scan_param['channels']]
  print '{'
  print '\"stations\" : [',
  for s in range(n_stations-1):
    print '\"'+stations[s]+'\", ',
  print '\"'+stations[n_stations-1]+'\"],'
  print '\"channels\" : [',
  for c in range(n_channels-1):
    print '\"'+channels[c]+'\", ',
  print '\"'+channels[n_channels-1]+'\"],'
  tmid = begin_time + timedelta(seconds = param.integration_time*n_integrations/2)
  tmid_tuple = tmid.utctimetuple()
  tmid_string = '%dy%dd%dh%dm%ds' %(tmid_tuple[0], tmid_tuple[7], tmid_tuple[3], tmid_tuple[4], tmid_tuple[5])
  print '\"epoch\" : \"'+tmid_string+'\",'
  print '\"clocks\" : {'
  for s in range(n_stations):
    print '    \"' + stations[s] +'\" : {'
    for c in range(n_channels):
      base_freq = param.freqs[scan_param['channels'][c][0]]
      sb =  -1 if scan_param['channels'][c][1] == 0 else 1
      delay = -delays[c,s] / (param.sample_rate)
      rate = -sb * rates[c,s]/(param.integration_time*(base_freq + sb * param.sample_rate/4) * n_integrations)
      snr_tot = snr[c,s]
      weight = weights[c,s]
      if (c < n_channels - 1):
        print '        \"' + channels[c] + '\" : [' + `delay` + ', ' + `rate` + ', ' + `snr_tot` + ', ' + `weight` + '],'
      else:
        print '        \"' + channels[c] + '\" : [' + `delay` + ', ' + `rate` + ', ' + `snr_tot` + ', ' + `weight` + ']'
    if(s < n_stations-1):
      print '    },'
    else:
      print '    }'
      print '}'
  print '}'

def clock_search(data, delays, rates, snr, param, options, scan_param, refstation):
  n_stations = len(scan_param['stations'])

  #Get initial estimate
  lag_offsets(data, n_stations, delays, rates)

  # Refine clocksearch
  ddelays = zeros(delays.shape)
  drates = zeros(rates.shape)
  for it in range(2):
    for station in range(n_stations):
      bldata = apply_model(data[:,station,:], refstation, station, delays, rates, param, scan_param)
      refine_offsets(bldata, station, ddelays, drates)
    delays += ddelays
    rates += drates

  # Compute SNR
  for station in range(n_stations):
    bldata = apply_model(data[:,station,:], refstation, station, delays, rates, param, scan_param)
    compute_snr(bldata, snr, station, param.sample_rate, param.integration_time)

def goto_scan(scan, inputfile, param, timeout):
  # Wait the correlation to reach the requested scan
  scan_param = param.get_scan_param(scan)
  if (timeout > 0):
    inputfile.seek(0, 2)
    oldpos = inputfile.tell()
    oldtime = datetime.utcnow()
    while scan_param == None:
      time.sleep(1)
      inputfile.seek(0, 2)
      newpos = inputfile.tell()
      newtime = datetime.utcnow()
      if newpos != oldpos:
        oldpos = newpos
        oldtime = newtime
      elif (newtime-oldtime).seconds > timeout:
        break
      scan_param = param.get_scan_param(scan)
  return scan_param

def get_parameters(vex, ctrl, corfile, scan_param, begin_time, end_time, timeout):
  scan = scan_param["name"]
  start_of_scan = vex_time.get_time(vex['SCHED'][scan]["start"])
  scan_duration = int(vex["SCHED"][scan]["station"][2].partition('sec')[0])
  if scan == param.get_scan_name(param.starttime):
    diff_time = begin_time - param.starttime
    scan_duration -= (param.starttime - start_of_scan).seconds
  else:
    diff_time = begin_time - start_of_scan
  scan_duration -= diff_time.seconds
  if end_time != None:
    scan_duration = min(scan_duration, (end_time - begin_time).seconds)
  slice_in_scan = int(diff_time.seconds / param.integration_time)
  start_slice = slice_in_scan + scan_param["start_slice"]
  start_fpos = scan_param["fpos"] + slice_in_scan * scan_param["slice_size"]
  if timeout > 0:
    # We are running concurrently with a correlator job
    if (end_time == None):
      end_time = vex_time.get_time(ctrl["stop"])
    diff_time = end_time - begin_time
    nseconds = min((diff_time.days*86400+diff_time.seconds), scan_duration)
    nslice = int(nseconds / param.integration_time)
  else:
    corfile.seek(0,2)
    pos = corfile.tell()
    nslice = min((pos - start_fpos + 1)/scan_param['slice_size'],
                 int(scan_duration / param.integration_time))
  return start_fpos, start_slice, nslice

def sighandler(signum, frame):
  # Data reading was interupted, compute delays and rates.
  # All variables are from global scope
  data_read = data[:,:,0:int_read[0],:]
  if int_read[0] > 0:
    clock_search(data_read, delays, rates, snr, param, options, scan_param, ref_index)
    write_clocks(vex, param, scan_param, delays, rates, snr, ref_index, begin_time, n_integrations)
  sys.exit(0)

####
######################## MAIN ##############################
####
(vex, ctrl, corfile, ref_station, options, param) = get_options()
timeout = options.timeout

# Determine which scan to clock search
if options.begin_time != None:
  begin_time = max(param.starttime, vex_time.get_time(options.begin_time))
  scan = param.get_scan_name(begin_time)
else:
  begin_time = param.starttime
  scan = param.get_scan_name(param.starttime)

# Go to current scan
scan_param = goto_scan(scan, corfile, param, timeout)
if scan_param == None:
  print >> sys.stderr, "Error : Premature end of correlation file"
  exit(1)

end_time = vex_time.get_time(options.end_time) if options.end_time != None else None
start_fpos, start_slice, n_integrations = get_parameters(vex, ctrl, corfile, scan_param, begin_time, end_time, timeout)

# initialize variables
channels = scan_param["channels"]
stations = scan_param["stations"]

vex_stations = [s for s in vex['STATION']]
vex_stations.sort()
try:
  ref_station_nr = vex_stations.index(ref_station)
except:
  print >> sys.stderr, 'Error : reference station ' + ref_station + ' is not found in vex file'
  sys.exit(1)
try:
  ref_index = stations.index(ref_station_nr)
except:
  print >> sys.stderr, 'Error : reference station ' + ref_station + ' is not found in the correlator output file'
  sys.exit(1)

n_channels = len(channels)
n_stations = len(stations)
delays = zeros([n_channels, n_stations])
rates = zeros([n_channels, n_stations])
snr = zeros([n_channels, n_stations])
#print [vex_stations[i] for i in param.stations]

# Skip to first integration
if not goto_integration(corfile, start_fpos, timeout):
  print >> sys.stderr, 'Error : Data ended prematurely'
  sys.exit(1)

# Register signal handler, in case the run is aborted
data = zeros([n_channels, n_stations, n_integrations, param.nchan + 1], dtype='c8')
int_read = [0]
signal.signal(signal.SIGINT, sighandler)
# Read the baseline data
try:
  read_integrations(corfile, data, int_read, param, scan_param, n_integrations, ref_index, timeout)
except EndOfData:
  # data file ended prematurely
  if  int_read[0] < 2:
    print >> sys.stderr, 'Error : Data ended prematurely'
    sys.exit(1)
  else:
    print >> sys.stderr, "Warning : data ended before the requested time. Proceeding with ", int_read[0], " datapoints."
  n_integrations = int_read[0]

# compute delays and rates
data = data[:,:,:n_integrations,:]
clock_search(data, delays, rates, snr, param, options, scan_param, ref_index)

############ Print output in JSON format ################
write_clocks(vex, param, scan_param, delays, rates, snr, ref_index, begin_time, n_integrations)
