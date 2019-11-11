#!/usr/bin/env python
import sys
import struct
import argparse
import json
import time
import os
import numpy as np
import gnuplotlib as pg
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
    dirname = dir_tuple.netloc + dir_tuple.path
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
  html.write("<html><head>\n"
    + "  <title>SFXC output - {} </title>\n".format(exper_name)
    + "  <style> BODY,TH,TD{font-size: 10pt } \n"
    + "    .popup_img a { position:relative; }\n"
    + "    .popup_img a img { position:absolute; display:none; top:20; height:200; z-index:99;}\n"
    + "    .popup_img a:hover img { display:block; }\n"
    + "</style>\n"
    + "</head> <body>\n")
  # Print preamble
  html.write("<a href='{}'>Vex file</a> -- \n".format(os.path.basename(vexfile)) \
    + "Scan name = {}, total averaging time = {} sec<br>".format(scan["name"], (ap * integr_time).total_seconds()) \
    + "Timerange: {}-{}\n".format(tstart, (tstart + ap * integr_time).time()))
   
  html.write("<div class='popup_img'>\n")
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
            + "<a href='{name}_large.png'>".format(name=basename)
            + "<img src='{name}.png' />".format(name=basename)
            + "{bbc_nr}</a></td>\n".format(bbc_nr=bbcs[station][channel1]))
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
        html.write("<a href='{name}_lag_large.png'>".format(name=basename)
          + "<img src='{name}_lag.png' />".format(name=basename)
          + "{:.1f}</a>".format(snr)
          + "<a href='{name}_ampl_large.png' >".format(name=basename)
          + "<img src='{name}_ampl.png' />".format(name=basename)
          + "A</a>"
          + "<a href='{name}_ph_large.png' >".format(name=basename)
          + "<img src='{name}_ph.png' />".format(name=basename)
          + "P</a><br><font size=-2>offset: {}</font>".format(offset)
          + "</td>\n")
      except KeyError:
        # Channel wasn't observed for this baseline
        html.write("  <td> <br> </td>\n")
    html.write("</tr>\n")
  html.write("</table>\n")
  html.write("</div>\n")

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
  channels in the experiment. Also a dict is generated which contains which 
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


def plot_thread(job):
  opts_small = {'terminal': 'png font ",9" size 300,200', '_with':'lines', 'unset': 'grid', '_set': 'format y "%1.0e"'}
  opts_small_phase = opts_small.copy()
  opts_small_phase['_set'] = 'format y "%3.0f"'
  opts_large = {'terminal': 'png font ",12" size 1024,768', '_with':'lines'}

  bl = job["bl"]
  if job["sb"] == 0: 
    vis = job["vis"][-1::-1] # Flip band for lower sideband
  else:
    vis = job["vis"]
  # frequencies
  f0 = job["freq"] + (job['sb'] - 1) * job['bw']
  f1 = f0 + job['bw']
  n = vis.size
  f = np.linspace(f0, f1, n)

  # auto-correlation
  if bl[0] == bl[1]:
    st = bl[0]
    sb = 'lsb' if job["sb"] == 0 else 'usb'
    pol = 'rcp' if job["pol1"] == 0 else 'lcp'
    ampl = abs(vis) / max(1, job["ap"])
    # plot spectrum (amplitude)
    basename = "{st}_{pol}_freq{freq}_{sb}".format(st=st, pol=pol, freq=job["freqnr"], sb=sb)
    name = job["outdir"] + '/' + basename + '.png'
    name_large = job["outdir"] + '/' + basename + '_large.png'
    plot = {'basename': basename}
    p = pg.gnuplotlib(output=name_large, xlabel='Freq [MHz]', 
                      ylabel='Amplitude', title=bl[0], **opts_large)
    p.plot(f, ampl)
    p = pg.gnuplotlib(output=name, xlabel='Channel', ylabel='Amplitude', 
                      title=bl[0], **opts_small)
    p.plot(ampl)
  else:
    # cross-correlations
    channel1 = (ch.freqnr, ch.sideband, ch.pol1)
    sb = 'lsb' if job["sb"] == 0 else 'usb'
    pol1 = 'rcp' if job["pol1"] == 0 else 'lcp'
    pol2 = 'rcp' if job["pol2"] == 0 else 'lcp'
    ampl = abs(vis) / max(1, job["ap"])
    basename = "{st1}_{pol1}_{st2}_{pol2}_freq{freq}_{sb}".format(st1=bl[0], st2=bl[1], pol1=pol1, pol2=pol2, freq=job["freqnr"], sb=sb)
    plot = {'basename': basename}
    # plot spectrum (amplitude)
    name = job["outdir"] + '/' + basename + '_ampl.png'
    name_large = job["outdir"] + '/' + basename + '_ampl_large.png'
    bl_str = "{}-{}".format(*bl)
    p = pg.gnuplotlib(output=name_large, xlabel='Freq [MHz]', 
                      ylabel='Amplitude', title=bl_str, **opts_large)
    p.plot(f, ampl)
    p = pg.gnuplotlib(output=name, xlabel='Channel', ylabel='Amplitude', 
                      title=bl_str, **opts_small)
    p.plot(ampl)

    # plot spectrum (phase)
    f0 = job["freq"]
    f = np.linspace(f0, f0 + job["bw"], n)
    name = job["outdir"] + '/' + basename + '_ph.png'
    name_large = job["outdir"] + '/' + basename + '_ph_large.png'
    phase = np.arctan2(np.imag(vis), np.real(vis)) * 180 / np.pi 
    p = pg.gnuplotlib(output=name_large, xlabel='Freq [MHz]', 
                      ylabel='Phase', title=bl_str, **opts_large)
    p.plot(f, phase)
    p = pg.gnuplotlib(output=name, xlabel='Channel', 
                      ylabel='Phase', title=bl_str, **opts_small_phase)
    p.plot(phase)

    # plot fringe (lag spectrum)
    t = np.arange(-(n-1), n-1)
    name = job["outdir"] + '/' + basename + '_lag.png'
    name_large = job["outdir"] + '/' + basename + '_lag_large.png'
    lags = abs(np.fft.fftshift(np.fft.irfft(vis)))
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
    p = pg.gnuplotlib(output=name_large, xlabel='Lag', ylabel='Ampl',
                      xmin=-(n-1), xmax=(n-2),
                      title=bl_str, **opts_large)
    p.plot(t, lags)
    p = pg.gnuplotlib(output=name, xlabel='Lag', ylabel='Ampl',
                      xmin=-(n-1), xmax=(n-2),
                      title=bl_str, **opts_small)
    p.plot(t, lags)
  return plot

def parse_args():
  parser = argparse.ArgumentParser(description='Create diagnostic HTML plotpage from SFXC output')
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
  output_file = out_tuple.netloc + out_tuple.path
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

