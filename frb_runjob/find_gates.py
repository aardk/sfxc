#!/usr/bin/env python3
import filterbank as fb
import json
import os
import numpy as np
import re
import subprocess

from argparse import ArgumentParser
from datetime import datetime, timedelta
from urllib.parse import urlsplit

from cor2filterbank import create_filterbank
from matplotlib import pyplot as pl
from numpy.polynomial.chebyshev import chebfit, chebval
from run_frb_jobs import run_job
from vex import Vex

def get_delay(f, year, day, hour, minute, second):
  line = f.readline()
  search = "{}y{:03d}d{:02d}h{:02d}m{:02d}".format(year, day, hour, minute, int(second))
  while line != "":
    if line.startswith(search):
      start = line.find('(') + 1
      end = start + line[start:].find(',')
      return float(line[start:end])
    line = f.readline()

  return 0

def bandpass(data, nif, thresh_power=2, thresh_std=1, debug = False):
    data /= np.median(data) # Normalize data first
    freq = data.mean(axis=0)
    stds = data.std(axis=0)
    nchan = len(freq) // nif
    # Mask band edges
    guard = max(1, nchan // 20)
    mask = np.array(
            [False if ((i%nchan) < guard) or ((i%nchan) >= (nchan - guard)) 
                   else True 
                   for i in range(len(freq))])
    ch = np.arange(nchan)
    for ifreq in range(nif):
        bchan = ifreq * nchan
        echan = (ifreq+1) * nchan
        band = freq[bchan:echan]
        med = np.median(band[guard:(nchan-guard)])
        mask2 = np.logical_and(mask[bchan:echan], np.less(band, med*1.5))
        cheb = chebfit(ch[mask2], band[mask2], 6)
        bp = chebval(ch, cheb)
        if (debug):
            pl.plot(band, label='data')
            pl.plot(bp, label='fit')
        band /= bp
        stds[bchan:echan] /= bp
        if (debug):
            pl.plot(band, label='bandpassed')
            pl.legend()
            pl.show()
    bandpass_med = np.median(freq[mask])
    # Filter out outliers first before computing std
    mask_outliers = np.logical_and(mask, np.less(freq, bandpass_med * 1.5))
    # Now filter out according to threshold
    bandpass_std = freq[mask_outliers].std()

    stds_med = np.median(stds[mask])
    stds_std = stds[mask].std()

    mask = np.logical_and(mask, np.less(freq, bandpass_med + bandpass_std * thresh_power))
    mask = np.logical_and(mask, np.greater_equal(freq, bandpass_med - bandpass_std * thresh_power))
    mask = np.logical_and(mask, np.less(stds, stds_med + stds_std * thresh_std))
    mask = np.logical_and(mask, np.greater_equal(stds, stds_med - stds_std * thresh_std))

    data[:, np.logical_not(mask)] = 1. # The median is equal to 1. due to how we normalize
    data /= data.mean(axis=0)

def mjd2date(mjdtime:float):
    t0 = datetime(1858, 11, 17)
    days = np.floor(mjdtime)
    fsec = (mjdtime - days) * 86400
    sec = int(np.floor(fsec))
    us = int(1e6 * (fsec - sec))
    return t0 + timedelta(days, sec, us)

def parse_polyco(fname):
    # TODO do more complete parse of polyco (with error handling)
    polyco = {}
    with open(fname, 'r') as f:
        line = f.readline().split()
        polyco['dm'] = float(line[4])
        line = f.readline().split()
        polyco['fmid'] = float(line[5])
    return polyco
  
def get_url_path(p):
    """ Because we abuse the URL format (by making file:// optional) sometimes a path will be
    in urlsplit.path and sometimes in urlsplit.netloc
    """
    x = urlsplit(p)
    return x.path if x.path != '' else x.netloc

def ensure_delay_file(vex_file_name, ctrl, station):
    """ Ensure that delay file exists, we need it to convert topocentric arrival times to geocentric time
    """
    exper = ctrl['exper_name']
    delay_dir = get_url_path(ctrl['delay_directory'])
    delay_file = os.path.join(delay_dir, f"{exper}_{station}.del")
    if (not os.path.isfile(delay_file)) or (os.stat(delay_file).st_size < 16):
        args = ['generate_delay_model', vex_file_name, station, delay_file]
        subprocess.run(args, check=True)
    return delay_file

def ensure_filterbank(vex, vex_file_name, ctrl, ctrl_file_name, station, machines, disable_create=False):
    outname = get_url_path(ctrl['output_file'])
    corfilename = outname + f"_{station}"
    filterbankname= os.path.splitext(outname)[0] + '.fb'

    exists = os.path.isfile(filterbankname)
    if exists:
        stat = os.stat(filterbankname)
    non_zero_exists = (exists) and (stat.st_size > 16)
    if (disable_create) and (not non_zero_exists):
        print(f'Error: Filterbank file creation is disabled but filterbank file "{filterbankname}" doesn\'t exist or is too small')
        exit(1)
    # if control file is newer than filterbank file we also recreate the filterbank file
    if (not non_zero_exists) or (stat.st_mtime < os.stat(ctrl_file_name).st_mtime):
        # Run correlator
        run_job(vex_file_name, ctrl_file_name, False, machines)
        setup_station = ctrl['setup_station']
        with open(corfilename, 'rb') as f:
            fbout = open(filterbankname, 'wb')
            create_filterbank(vex, [f], fbout, 0, -1, 3, setup_station)
            fbout.close()
    
    return filterbankname

def vextime(tstring):
    year, doy, hour, minute, second = (int(x) for x in re.split('y|d|h|m|s', tstring)[:-1])
    nsec = 60 * (60 * hour + minute) + second
    t0 = datetime(year, 1, 1)
    return t0 + timedelta(days=doy-1, seconds=nsec)

def get_scan(ctrl, vex):
    integr_start = vextime(ctrl['start'])
    for scanname, scan in vex['SCHED'].items():
        tstart = vextime(scan['start'])
        maxsec = 0
        for station in scan.getall('station'):
            sec = int(station[2].split()[0])
            maxsec = max(maxsec, sec)
        tend = tstart + timedelta(seconds=sec)
        if (integr_start >= tstart) and (integr_start < tend):
            return scanname
    print(f"Error: couldn't find scan belonging to time {integr_start}")
    exit(1)

def get_if(ctrl, vex, scan):
    integr_start = vextime(ctrl['start'])
    setup_station = ctrl['setup_station']
    mode = vex['SCHED'][scan]['mode']
    try:
        freq_id, = (freq[0] for freq in vex['MODE'][mode].getall('FREQ') if setup_station in freq[1:])
    except ValueError:
        print(f"Error: could'n find $FREQ definition of station {setup_station} in mode {mode}")
        exit(1)

    # create pairs (Freq_edge, sideband)
    ifs = set((line[1].split()[0], line[2]) for line in vex['FREQ'][freq_id].getall('chan_def'))
    print(f"ifs = {ifs}")
    return len(ifs)

def get_pulse(data, noise_power, noise_std, istart):
    N = len(data)
    W = np.arange(1, N // 4)
    pulse = np.zeros(len(W)) 
    snr = np.zeros(len(W))
    for i in range(len(W)):
        w = W[i]
        #  Get boxcar average data, we should realy use convolve for this
        d = np.array([data[k:k+w].sum() for k in range(0, N-w)])
        j = d.argmax()
        pulse[i] = j
        snr[i] = (d[j] - noise_power * w) / (np.sqrt(w) * noise_std)
    best = snr.argmax()
    #pl.plot(snr);
    #pl.plot(pulse);
    #pl.plot(W); pl.show()
    pulse_start = int(round(istart + pulse[best]))
    pulse_end = int(round(pulse_start + W[best]))
    pulse_snr = snr[best]
    return pulse_start, pulse_end, pulse_snr

def get_search_window(vex, ctrl, header, source, delay_file_name):
    bursttime = mjd2date(ctrl['burst_mjd'])
    polyco = parse_polyco(get_url_path(ctrl['pulsars'][source]['polyco_file']))

    DM = polyco['dm']
    fmid = polyco['fmid']
    # FIXME get top frequency from vexfile
    disp_delay = 4.15e3 * DM * (1 / fmid**2 - 1/1510.**2)

    # FIXME make this conditional
    delay_plot_name = os.path.splitext(delay_file_name)[0] + '.txt'
    args = ['plot_delay_table', delay_file_name, delay_plot_name, '-n', '0']
    subprocess.run(args, check=True)

    delays = open(delay_plot_name, 'r')
    year = bursttime.year
    day = (bursttime - datetime(year, 1, 1)).days + 1
    WIDTH = 0.1

    # Determine a priory pulse location from arrival times
    hour, minute, second = bursttime.hour, bursttime.minute, bursttime.second + bursttime.microsecond / 1000000.
    delay = get_delay(delays, year, day, hour, minute, second)
    new_second = second - delay + disp_delay

    # Define rectangular window around pulse in which we conduct our search
    # Start of data in seconds, lets not worry about scans around midnight here
    tstart = round((header['tstart'] - np.floor(header['tstart'])) * 86400)
    dt = header['tsamp']
    N = int(WIDTH / dt)
    print(f"WIDTH = {WIDTH}, dt = {dt}, N={N}")
    # The pulse location within the data set
    pulse_location = hour * 3600 + minute *60 + new_second
    pulse_offset = pulse_location - tstart # [seconds]
    W0 = max(int(pulse_offset / dt) - N // 2, 0)

    return W0, N

def gatesearch(vex_file, ctrl_file, interactive, machines, bif, eif, down_freq, down_time, disable_create, write_ctrl, gatesearch_filename=''):
    vex = Vex(vex_file)
    with open(ctrl_file, 'r') as f:
        ctrl = json.load(f)

    search_ctrl_file = get_url_path(ctrl['search_ctrl_file'])
    with open(search_ctrl_file, 'r') as f:
        search_ctrl = json.load(f)

    station = search_ctrl['stations'][0]
    if len(search_ctrl['stations']) > 1:
        print(f"Warning: Only doing pulse search for first station {station}")
    
    scan = get_scan(ctrl, vex)
    
    # Ensure delay model exists
    delay_file_name = ensure_delay_file(vex_file, search_ctrl, station)

    # Ensure filterbank files exist
    fb_file_name = ensure_filterbank(vex, vex_file, search_ctrl, search_ctrl_file, station, machines, disable_create)
    header, bheader, data = fb.loaddata(fb_file_name)

    # NB: We assume all IF from vex file are in the filterbank file (NB: frequency axis is reversed)
    nif = get_if(ctrl, vex, scan)
    nchan = data.shape[-1] // nif
    bchan = 0 if eif == -1 else (nif - eif - 1) * nchan
    echan = (nif - bif) * nchan

    # Bandpass and blank RFI channels
    bandpass(data, nif, debug=False)
    
    if gatesearch_filename == '':
        gatesearch_filename = os.path.splitext(fb_file_name)[0] + '.pdf'

    # Find time range around pulse (taking into acount geometric delay and dispersion) to search
    source = vex['SCHED'][scan]['source']
    win_start, Nwin = get_search_window(vex, search_ctrl, header, source, delay_file_name)
    win_end = win_start + Nwin
    # Do the actual search
    d2 = data[:, bchan:echan].mean(axis=1)
    pulse_begin, pulse_end, pulse_snr = get_pulse(d2[win_start:win_end], d2.mean(), d2.std(), win_start)

    integr_start = vextime(ctrl['start'])
    tsamp = header['tsamp']
    # Fake pulsar has period 2 seconds
    if integr_start.second % 2 == 0:
        gate_begin = pulse_begin * tsamp / 2
        gate_end = pulse_end * tsamp / 2 
    else:
        gate_begin = pulse_begin * tsamp / 2 + 0.5
        gate_end = pulse_end * tsamp / 2 + 0.5
        if gate_begin >= 1.0:
            gate_begin = (gate_begin % 1.)
            gate_end = (gate_end % 1.)

    print(f'Pulse_start = {pulse_begin,} , pulse_end = {pulse_end,}, pulse_snr = {pulse_snr}')
    print(f'            "interval": [')
    print(f'                           {gate_begin},')
    print(f'                           {gate_end}')
    print(f'                        ],')
    ctrl['pulsars'][source]['interval'] = [gate_begin, gate_end]
    if (write_ctrl):
        with open(ctrl_file, 'w') as f:
            json.dump(ctrl, f, indent=4)
     
    # plot results
    create_plots(data, tsamp, pulse_begin, pulse_end, pulse_snr, bchan, echan, down_freq, down_time, interactive, gatesearch_filename)

def create_plots(data, tsamp, pulse_begin, pulse_end, pulse_snr, bchan, echan, down_freq, down_time, interactive, gatesearch_filename, pixels_pulse = 4, Ntime = 128):
    # Down sample pulse profile, this makes weak detections easier to spot. Default to at most 4 pixels over burst
    # NB: down_freq == 0, down_time == 0 signals auto-scale
    if down_time == 0:
        Ndown = max(1, (pulse_end - pulse_begin) // pixels_pulse)
    else:
        Ndown = down_time
    
    if down_freq == 0:
        Mdown = max(1, (echan - bchan) // 128)
    else:
        Mdown = down_freq
    M = (echan - bchan) // Mdown
    
    if (pulse_begin <  Ndown * Ntime // 2):
        W0 = 0
        W1 = min(len(data), Ndown * Ntime)
    elif (len(data) - pulse_begin < Ndown * Ntime // 2):
        W1 = len(data)
        W0 = max(0, W1 - Ndown * Ntime)
    else:
        W0 = pulse_begin - Ndown * Ntime // 2
        W1 = W0 +  Ndown * Ntime
    
    N = (W1 - W0) // Ndown
    t = np.arange(W0, W0 + N*Ndown, Ndown) * tsamp + Ndown/2 * tsamp
    
    data2 = data[W0:(W0+N*Ndown), bchan:(bchan+M*Mdown)].reshape([N, Ndown, M, Mdown]).mean(axis=3).mean(axis=1)

    fig, (ax1, ax2) = pl.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [4,1]}, figsize=[5.8, 8.3])
    # Waterfall plot
    ax1.imshow(data2.T, extent=[t[0], t[-1], bchan, bchan+M*Mdown])
    ax1.set_aspect('auto')
    ax1.set_title(f'SNR = {pulse_snr:.2f}, downsample = {Ndown}, {Mdown}')
    # Pulse profile integrated over frequency
    profile = data2.mean(axis = 1)
    ax2.plot(t, profile)
    ax2.set_title(f'Pulse profile')
    # Add Guiding lines around burst in all three plots
    ax1.axvline(x = pulse_begin * tsamp, color='r', linestyle=':')
    ax1.axvline(x = pulse_end * tsamp, color='r', linestyle=':')
    ax2.axvline(x = pulse_begin * tsamp, color='r')
    ax2.axvline(x = pulse_end * tsamp, color='r')
    ax2.axhline(y = profile.mean() + 3*profile.std(), color='g')
    # We always save a copy of the plot
    pl.savefig(gatesearch_filename)
    if (interactive):
        pl.show()
    pl.close()

if __name__ == "__main__":
    p = ArgumentParser(description="Perform gate seacrh for FRB, assuming ctrl file was produced by create_frb_jobs.py")
    p.add_argument("-i", "--interactive",
                    default = False,
                    action = "store_true",
                    help="Perform interactive gate search (work in progress, currently only shows plots)")
    p.add_argument("-d", "--disable-create",
                    default = False,
                    action = "store_true",
                    help="Disable creation of new filterbank files when control file is updated")
    p.add_argument("-f", "--down-freq",
                    type=int,
                    default = 0,
                    help="Downsample frequency for the pulse plots, default=0 (auto)")
    p.add_argument("-t", "--down-time",
                    type=int,
                    default = 0,
                    help="Downsample frequency for the pulse plots, default=0 (auto)")
    p.add_argument("-b", "--bif",
                    type=int,
                    default = 0,
                    help="First IF (as is defined in the vex file) to be used in pulse search")
    p.add_argument("-e", "--eif",
                    type=int,
                    default = -1,
                    help="Last IF (as is defined in the vex file) to be used in pulse search")
    p.add_argument("-m", "--machines",
                   default="k", type=str,
                   help="Machines to run correlator nodes on, default='k', allowed values are k,l,m,n,o ; it is also possible to list individual machines e.g. k1,k2")
    p.add_argument("-s", "--simulate",
                    default = False,
                    action = "store_true",
                    help="Don't write results to control file")
    p.add_argument("vex", type=str, help="VEX file of experiment")
    p.add_argument("ctrlfiles", type=str, nargs='+', help="SFXC control files generated by create_frb_jobs.py")
    args = p.parse_args()
    machines = args.machines.split(',')
    for ctrlname in args.ctrlfiles:
        gatesearch(args.vex, ctrlname, args.interactive, machines, args.bif, args.eif, args.down_freq, args.down_time, args.disable_create, not args.simulate)
