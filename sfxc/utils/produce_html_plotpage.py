#!/usr/bin/env python
import sys
import struct
import argparse
import json
import time
import os
import numpy as np
# NB: ppgplot also works with Giza which is preferred over classic PGPLOT, 
# Giza uses the CAIRO backend produces which creates much nicer looking plots.
import ppgplot as pp
from multiprocessing import Process, Pool
from datetime import datetime, timedelta
from urlparse import urlparse
from shutil import copyfile
from vex import Vex 
# These two modules are part of sfxc2mark4
from sfxcdata import SFXCData
from experiment import experiment

NUMBER_THREADS = 8

def write_html(vexfile, ctrl, scan, exper, global_header, integr, stats, tstart, ap):
  stations = sorted(set([x[0] for x in integr.keys()]))
  baselines = sorted(filter(lambda x:x[0] != x[1], integr.keys()))
  setup_station = ctrl["setup_station"] if "setup_station" in ctrl else ctrl["stations"][0]
  exper_name = global_header.exper.rstrip('\0')

  channels, bbcs = create_channel_mapping(exper, scan["mode"], setup_station)
  try:
      dir_tuple = urlparse(ctrl['html_output'])
      dirname = dir_tuple.netloc if dir_tuple.path == '' else dir_tuple.path
  except KeyError:
      dirname = '.'
  imgdir = dirname + '/' + scan["name"]
  if not os.path.exists(imgdir):
    os.makedirs(imgdir)
  htmlname = imgdir + "/{}_{}.html".format(exper_name, scan["name"])
  html = open(htmlname, 'w')
  newvex = imgdir + "/" + os.path.basename(vexfile)
  copyfile(vexfile, newvex) 

  integr_time =  timedelta(0, 
                           global_header.integr_time // 1000000, 
                           global_header.integr_time % 1000000)
  plots = generate_plots(integr, channels, imgdir, ap, integr_time.total_seconds())
  
  # HTML HEAD
  html.write("<html><head>\n" \
    + "  <title>SFXC output - {} </title>\n".format(exper_name) \
    + "  <style> BODY,TH,TD{font-size: 10pt } \n" \
    + "#fimage {\n" \
    + "position: fixed;\n" \
    + "top: 50%;\n" \
    + "border: 1px solid black;\n" \
    + "margin-left: 7px;\n" \
    + "}\n" \
    + "</style>\n" \
    + "</head> <body>\n")
  # Jay's code to make images float
  html.write("<script language=\"JavaScript\"><!--\n" \
    + "function show(imageSrc) {\n" \
    + "  if (document.images) document.images['plot_image'].src = imageSrc;\n" \
    + "}\n" \
    + "//--></script>\n\n")

  # Print preamble
  #html.write("<a href='{}'>Vex file</a> -- \n".format(imgdir + '/' + newvex) \
  html.write("<a href='{}'>Vex file</a> -- \n".format(os.path.basename(vexfile)) \
    + "Scan name = {}, total averaging time = {} sec<br>".format(scan["name"], (ap * integr_time).total_seconds()) \
    + "Timerange: {}-{}\n".format(tstart, (tstart + ap * integr_time).time()))
  
  # Print table
  html.write("<table border=1 bgcolor='#dddddd' cellspacing=0>\n")
  # First row
  html.write("<tr>\n" \
    + "  <th rowspan=2> {} </th>\n".format(exper_name) \
    + "  <th colspan={}>Auto correlations (BBC number)</th>\n".format(len(stations)) \
    + "  <th colspan={}>Cross correlations (SNR, lag offset)</th>\n".format(len(baselines)) \
    + "</tr>\n")
  #second row
  html.write("<tr>\n")
  for station in stations:
    html.write("<th>{}</th>".format(station))
  for bl in baselines:
    html.write("<th>{}-{}</th>".format(*bl))
  # floating image
  bl = (stations[0], stations[0])
  filename = plots[bl].values()[0]["basename"] + ".png"
  html.write("<td rowspan=99><div id='fimage'>"
    + '<img src="{}" name="plot_image">'.format(filename)
    + '</div></td>\n'
    + '</tr>\n')

  # Create sorted list of channels, sorting order polarizations: RCP-RCP, RCP-LCP, LCP-LCP, LCP-RCP
  def srt_chan(a,b):
    c = cmp(a[:-1], b[:-1])
    if c==0:
      if a[-2] == 0:
        return cmp(a[-1], b[-1])
      else:
        return cmp(b[-1], a[-1])
    return c

  channels_in_data = set()
  for bl in baselines:
    for ch in integr[bl].keys():
      channels_in_data.add(ch)
  print channels_in_data
  channels_in_data = sorted(channels_in_data, cmp=srt_chan)

  # print all (auto-) correlations
  for ch in channels_in_data:
    channel1 = (ch.freqnr, ch.sideband, ch.pol1)
    html.write("<tr>\n")
    html.write("<th>")
    freq = channels[channel1]["freq"]
    sb = 'LSB' if ch.sideband == 0 else 'USB'
    pol1 = 'RCP' if ch.pol1 == 0 else 'LCP'
    pol2 = 'RCP' if ch.pol2 == 0 else 'LCP'
    html.write("{:.2f}MHz, {}, {}-{}".format(freq, sb, pol1, pol2))
    html.write("</th>\n")
    # Auto correlations
    if pol1 == pol2:
      for station in stations:
        bl = (station, station)
        try:
          basename = plots[bl][ch]["basename"]
          html.write("  <td>" 
            + "<a href = '{name}_large.png' OnMouseOver=\"show('{name}.png');\">".format(name=basename)
            + "{}</a></td>\n".format(bbcs[station][channel1]))
        except KeyError:
          # Station didn't observe channel
          html.write("  <td> </td>\n")
    else:
      html.write("  <td colspan='{}'>Cross hands</td>\n".format(len(stations)))

    # Cross correlations
    for bl in baselines:
      try:
        basename = plots[bl][ch]["basename"]
        snr = plots[bl][ch]["snr"]
        offset = plots[bl][ch]["lag"]
        if pol1 == pol2:
          # For parallel-hands the background colour indicates SNR
          if snr >= 6:
            # SNR 6: dark green; SNR >=9 light green
            green = min(int(100 + 155 * (snr-6) / 3.), 255)
            colour = "#00{:02x}00".format(green)
          else:
            # SNR 5: dark red; SNR <=3 light red
            red = min(int(100 + 155 * (6-snr) / 3.), 255)
            colour = "#{:02x}0000".format(red)
          html.write("  <td bgcolor='{}'>".format(colour))
        else:
          html.write("  <td>" )
        html.write("<a href = '{name}_lag_large.png' OnMouseOver=\"show('{name}_lag.png');\">".format(name=basename)
            + "{:.1f}</a>".format(snr)
          + "<a href = '{name}_ampl_large.png' OnMouseOver=\"show('{name}_ampl.png');\">".format(name=basename)
          + "A</a>"
          + "<a href = '{name}_ph_large.png' OnMouseOver=\"show('{name}_ph.png');\">".format(name=basename)
          + "P</a><br><font size=-2>offset: {}</font>".format(offset)
          + "</td>\n")
      except KeyError:
        # Channel wasn't observed for this baseline
        html.write("  <td> <br> </td>\n")
    html.write("</tr>\n")
  html.write("</table>\n")

  # Sampler statistics
  html.write("<h1> Sampler statistics </h1><br>\n")
  for station in stations:
    html.write("<table border=1 bgcolor='#dddddd' cellspacing=0>\n"
      + "<tr>\n"
      + "  <th>{}</th>\n".format(station)
      + "  <th> - - </th>\n"
      + "  <th> - + </th>\n"
      + "  <th> + - </th>\n"
      + "  <th> + + </th>\n"
      + "  <th> invalid </th>"
      + "</tr>\n")
    for ch in sorted(stats[station].keys()):
      freq = channels[ch]["freq"]
      sb = 'LSB' if ch.sideband == 0 else 'USB'
      pol = 'RCP' if ch.pol == 0 else 'LCP'
      html.write("<tr>\n"
        + "  <th>{:.2f}MHz, {}, {}</th>".format(freq, sb, pol))
      for level in stats[station][ch]:
        html.write(" <td>{:.2f}%</td> ".format(100 * level / ap))
      html.write("\n</tr>\n")
    html.write("</table><br>\n")
                                             
def create_channel_mapping(experiment, mode, setup_station):
  """ Create dict which gives the frequencies, and bandwidths for all
  channels in the experiment. Als a dict is generated which contains which 
  BBC corresponds to each channel of the reference station."""
  channels = {}
  bbcs = {}
   
  for station in experiment.modes[mode]['bykey']:
    bbcs[station] = {}
    for setup_ch_key in experiment.modes[mode]['bykey'][setup_station]:
      setup_ch = experiment.modes[mode]['bykey'][setup_station][setup_ch_key]
      channels[setup_ch_key] = {"freq": setup_ch['freq'], "bw": setup_ch['bw']}
      sb = 1 if setup_ch['sb'] == 'U' else 0
      setup_fmin = setup_ch['freq'] + (sb-1) * setup_ch['bw']
      setup_fmax = setup_ch['freq'] + sb * setup_ch['bw']
      for chkey in experiment.modes[mode]['bykey'][station]:
        ch = experiment.modes[mode]['bykey'][station][chkey]
        if setup_ch['pol'] != ch['pol']:
          continue
        sb = 1 if ch['sb'] == 'U' else 0
        fmin = ch['freq'] + (sb-1) * ch['bw']
        fmax = ch['freq'] + sb * ch['bw']
        # Two possible matches: 
        # 1. The channel is entirely contained within one of the 
        #     setup station's channels
        # 2. The setup channel is entirely contained within the
        #    current channel 
        if ((fmin >= setup_fmin) and (fmax <= setup_fmax)) or \
            ((setup_fmin>= fmin) and (setup_fmax <= fmax)):
          bbcs[station][setup_ch_key] = ch["bbc_nr"]
          break
  return channels, bbcs

def generate_plots(integr, channels, outdir, ap, integr_time):
  # Create sorted list of stations and baselines
  stations = sorted(set([x[0] for x in integr.keys()]))
  baselines = sorted(filter(lambda x:x[0] != x[1], integr.keys()))
  pool = Pool(processes = NUMBER_THREADS)
 
  jobs = []
  for bl in integr:
    for ch in integr[bl]:
      key = ch[:3]
      channel = channels[key]
      job = {'bl': bl, 'freqnr': ch.freqnr, "sb": ch.sideband}
      job["vis"] = integr[bl][ch]
      job["freq"] = channel["freq"]
      job["bw"] = channel["bw"]
      job["ap"] = ap
      job["integr_time"] = integr_time
      job["outdir"] = outdir
      if len(ch) == 3:
        job["pol1"] = ch.pol
      else:
        job["pol1"] = ch.pol1
        job["pol2"] = ch.pol2
      jobs.append(job)
  results = pool.map(plot_thread, jobs)
  plots = {}
  i = 0
  for bl in integr:
    plots[bl] = {}
    for ch in integr[bl]:
      plots[bl][ch] = results[i]
      i += 1
  return plots

#def plot_thread(bl, ch, freq, bw, vis, outdir, ap, integr_time):
def plot_thread(job):
  bl = job["bl"]
  # auto-correlation
  if bl[0] == bl[1]:
    st = bl[0]
    sb = 'lsb' if job["sb"] == 0 else 'usb'
    pol = 'rcp' if job["pol1"] == 0 else 'lcp'
    data = abs(job["vis"]) / max(1, job["ap"])
    n = data.size
    # plot spectrum (amplitude)
    f0 = job["freq"]
    f = np.linspace(f0, f0 + job["bw"], n)
    basename = "{st}_{pol}_freq{freq}_{sb}".format(st=st, pol=pol, freq=job["freqnr"], sb=sb)
    name = basename + '.png'
    name_large = basename + '_large.png'
    plot = {'basename': basename}
    label = "{st}-{st}-f{freq}-{sb}".format(st=st, freq=job["freqnr"], sb=sb)
    pp.pgopen(job["outdir"] + '/' + name_large + '/png')
    pp.pgpap(11, 0.75)
    pp.pgenv(f[0], f[-1], min(data), max(data), 0, 1)
    pp.pglab('Freq [MHz]', 'Amplitude', '{}'.format(bl[0]))
    pp.pgline(f, data)
    pp.pgclos()
    pp.pgopen(job["outdir"] + '/' + name + '/png')
    pp.pgpap(5, 0.75)
    pp.pgenv(0, n, min(data), max(data), 0, 1)
    pp.pglab('Channel', 'Amplitude', '{}'.format(bl[0]))
    pp.pgline(np.arange(n), data)
    pp.pgclos()
  else:
    # cross-correlations
    channel1 = (ch.freqnr, ch.sideband, ch.pol1)
    sb = 'lsb' if job["sb"] == 0 else 'usb'
    pol1 = 'rcp' if job["pol1"] == 0 else 'lcp'
    pol2 = 'rcp' if job["pol2"] == 0 else 'lcp'
    data = abs(job["vis"]) / max(1, job["ap"])
    n = data.size
    basename = "{st1}_{pol1}_{st2}_{pol2}_freq{freq}_{sb}".format(st1=bl[0], st2=bl[1], pol1=pol1, pol2=pol2, freq=job["freqnr"], sb=sb)
    plot = {'basename': basename}
    # plot spectrum (amplitude)
    f0 = job["freq"]
    f = np.linspace(f0, f0 + job["bw"], n)
    name = basename + '_ampl.png'
    name_large = basename + '_ampl_large.png'
    bl_str = "{}-{}".format(*bl)
    label = "{st1}-{st2}-f{freq}-{sb}".format(st1=bl[0], st2=bl[1], freq=job["freqnr"], sb=sb)
    pp.pgopen(job["outdir"] + '/' + name_large + '/png')
    pp.pgpap(11, 0.75)
    pp.pgenv(f[0], f[-1], min(data), max(data), 0, 1)
    pp.pglab('Freq [MHz]', 'Amplitude', bl_str)
    pp.pgline(f, data)
    pp.pgclos()
    pp.pgopen(job["outdir"] + '/' + name + '/png')
    pp.pgpap(5, 0.75)
    pp.pgenv(0, n, min(data), max(data), 0, 1)
    pp.pglab('Channel', 'Amplitude', bl_str)
    pp.pgline(np.arange(n), data)
    pp.pgclos()

    # plot spectrum (phase)
    f0 = job["freq"]
    f = np.linspace(f0, f0 + job["bw"], n)
    name = basename + '_ph.png'
    name_large = basename + '_ph_large.png'
    data = job["vis"]
    phase = np.arctan2(np.imag(data), np.real(data)) * 180 / np.pi 
    label = "{st1}-{st2}-f{freq}-{sb}".format(st1=bl[0], st2=bl[1], freq=job["freqnr"], sb=sb)
    pp.pgopen(job["outdir"] + '/' + name_large + '/png')
    pp.pgpap(11, 0.75)
    pp.pgenv(f[0], f[-1], min(phase), max(phase), 0, 1)
    pp.pglab('Freq [MHz]', 'Phase', bl_str)
    pp.pgline(f, phase)
    pp.pgclos()
    pp.pgopen(job["outdir"] + '/' + name + '/png')
    pp.pgpap(5, 0.75)
    pp.pgenv(0, n, min(phase), max(phase), 0, 1)
    pp.pglab('Channel', 'Phase', bl_str)
    pp.pgline(np.arange(n), phase)
    pp.pgclos()

    # plot fringe (lag spectrum)
    t = np.arange(-(n-1), n-1)
    name = basename + '_lag.png'
    name_large = basename + '_lag_large.png'
    lags = abs(np.fft.fftshift(np.fft.irfft(data)))
    lag = lags.argmax() - n + 1
    plot["lag"] = lag
    # To estimate SNR compute theoreticalnoise in case of perfect sampler
    # stats. Also in the correction for non-centerred fringes we assume a 
    # triangular lag window
    # FIXME: Handle 1bit, and 1bit against 2bit case
    sample_rate = 2 * job["bw"] * 1e6
    # correct for non-centered fringe
    lag_offset = 1. - abs(lag)/float(n-1) if (lag < n-1) else 1.
    # factor 2 in denomenator accounts for window function
    noise = np.sqrt(ap) / (2 * 0.881 * np.sqrt(sample_rate * job["integr_time"] * lag_offset))
    plot["snr"] = lags.max() / noise
    label = "{st1}-{st2}-f{freq}-{sb}".format(st1=bl[0], st2=bl[1], freq=job["freqnr"], sb=sb)
    pp.pgopen(job["outdir"] + '/' + name_large + '/png')
    pp.pgpap(11, 0.75)
    pp.pgenv(t[0], t[-1], min(lags), max(lags), 0, 0)
    pp.pglab('Lag', 'Ampl', bl_str)
    pp.pgline(t, lags)
    pp.pgclos()
    pp.pgopen(job["outdir"] + '/' + name + '/png')
    pp.pgpap(5, 0.75)
    pp.pgenv(t[0], t[-1], min(lags), max(lags), 0, 0)
    pp.pglab('Lag', 'Ampl', bl_str)
    pp.pgline(t, lags)
    pp.pgclos()
  return plot
  
def parse_args():
  parser = argparse.ArgumentParser(description='Convert SFXC output to mark4 format')
  parser.add_argument("vexfile", help='vex file')
  parser.add_argument("ctrlfile", help='SFXC control file used in the correlation')
  args = parser.parse_args()
  vex = Vex(args.vexfile)
  ctrl = json.load(open(args.ctrlfile, 'r'))
  return vex, args.vexfile, ctrl

#########
########################## MAIN #################################3
########
if __name__ == "__main__":
  vex, vexfilename, ctrl = parse_args()
  exper = experiment(vex)
  out_tuple = urlparse(ctrl['output_file'])
  output_file = out_tuple.netloc if out_tuple.path == '' else out_tuple.path
  data = SFXCData(output_file)
  nchan = data.nchan

  # Dummy scan for first iteration
  scan = {"stop": datetime(1,1,1)}
  ap = 0
  while True:
    # Check if we need to move to next scan
    if data.current_time() >= scan["stop"]:
      if ap > 0:
        write_html(vexfilename, ctrl, scan, exper, data.global_header, integr, stats, tstart, ap)
      # initialise scan
      scan = exper.get_scan(vex, ctrl, data)
      tstart = data.current_time()
      integr = {}
      for bl in data.vis:
        integr[bl] = {ch:np.zeros(nchan+1, dtype='c16') for ch in data.vis[bl]}
      stats = {}
      for station in data.stats:
        stats[station] = {ch:np.zeros(5, dtype=float) for ch in data.stats[station]}
      ap = 0

    for station in data.stats:
      for ch in data.stats[station]:
          st = data.stats[station][ch]
          n = float(sum(st))
          stats[station][ch] += np.array(st) / n
    for bl in data.vis:
      for ch in data.vis[bl]:
        integr[bl][ch] += data.vis[bl][ch].vis
    ap += 1

    if not data.next_integration():
      break
  # write final scan
  if ap > 0:
    write_html(vexfilename, ctrl, scan, exper, data.global_header, integr, stats, tstart, ap)

