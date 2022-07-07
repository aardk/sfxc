#!/usr/bin/python

# Copyright (c) 2009 Joint Institute for VLBI in Europe (Netherlands)
# All rights reserved.
#  
# Author(s): Aard Keimpema <keimpema@jive.nl>
# 

import sys, os
import subprocess
import multiprocessing
from optparse import OptionParser
from vex import Vex 

import json

def gen_delay_tables(vex_file, nprocs, delay_directory, stations, exper_name):
  nstations = len(stations)
  if nprocs != None:
    nprocs = min(nstations, nprocs)
  else:
    # If nprocs is not specified make it equal to the number of cpu cores
    try:
      ncpu = multiprocessing.cpu_count()
      nprocs = min(nstations, ncpu)
    except NotImplemented:
      nprocs = nstations
  procs = [(None, None) for i in range(nprocs)]
  for i in range(nprocs):
    station = stations[i]
    procs[i] = (station, gen_delay(delay_directory, exper_name, vex_file, station))
  for i in range(nprocs, nstations):
    station = stations[i]
    done = False
    # TODO this should be done using multiprocessing
    while not done:
      for j in range(nprocs):
        status = procs[j][1].poll()
        if status != None:
          if status < 0:
            print "Error failed generating delay file for station", procs[j][0]
            sys.exit(1)
          else:
            done = True
            procs[j] = (station, gen_delay(delay_directory, exper_name, vex_file, station))
            break
  # Wait for processes to finish
  for j in range(nprocs):
    status = procs[j][1].wait()
    if status < 0:
      print "Error failed generating delay file for station", procs[j][0]

def gen_delay(delay_directory, exper_name, vex_file, station):
  delay_file = delay_directory + "/" + exper_name + "_" + station + ".del"
  cmd = ["generate_delay_model", vex_file, station , delay_file]
  print cmd
  proc = subprocess.Popen(cmd)
  print "Generating delay model for station " + station
  return proc


######### MAIN CODE
if __name__ == "__main__":
  usage = "Usage: %prog [options] <vex file> [control file]"
  parser = OptionParser(usage=usage)
  parser.add_option("-p", "--number-process", dest="nprocs", 
                    help="Specify the number of generate jobs to be run simultaneously; Default: Number of cpu cores",
                    action="store", type="int")
  parser.add_option("-d", "--dir", dest="delay_directory", 
                    help="Output directory for delay files; Default: \"delay_directory\" from ctrl file or if no ctrl file is specified current directory",
                    action="store", type="str")
  parser.add_option("-s", "--stations", dest="stations", 
                    help="Comma separated list of station two letter codes; Default: \"stations\" from ctrl file or if no ctrl file is specified all stations in vex file",
                    action="store", type="str")
  parser.add_option("-e", "--exp", dest="exper_name", 
                    help="Experiment name; Default: \"exper_name\" from ctrl file or when absent it is taken from vex file",
                    action="store", type="str")
  opts, args = parser.parse_args()
  if len(args) == 1:
    ctrl = None
  elif len(args) == 2:
    with open(args[1], 'r') as f:
      ctrl = json.load(f)
  else:
    parser.error("Invalid number of arguments")

  vex = Vex(args[0])

  if opts.stations != None:
     stations = opts.stations.split(',') 
  elif ctrl != None:
    try:
      stations = ctrl["stations"]
    except KeyError:
      print "Error: Control file does not have a stations field"
      exit(1)
  else:
     stations = [x for x in vex['STATION']]
  
  if opts.delay_directory != None:
    delay_directory = opts.delay_directory
  elif ctrl != None:
    try:
      delay_directory = ctrl["delay_directory"]
      if delay_directory[:7] != "file://":
        print "WARNING : delay_directory should start with file://"
      elif delay_directory == "file://":
        delay_directory = "."
      else:
        delay_directory = delay_directory[7:]
    except KeyError:
      # /tmp is default dir where SFXC looks for delay files
      delay_directory = "/tmp"
  else:
      delay_directory = "."
  
  if opts.exper_name != None:
    exper_name = opts.exper_name
  elif ctrl != None:
    try:
      exper_name = ctrl["exper_name"]
    except KeyError:
      exper_name = vex['GLOBAL']['EXPER']
  else:
    exper_name = vex['GLOBAL']['EXPER']
   
  gen_delay_tables(args[0], opts.nprocs, delay_directory, stations, exper_name)
