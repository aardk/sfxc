#! /usr/bin/python

# Standard Python modules.
import os, re, stat, sys, time
import subprocess
import urlparse
import optparse

# The json module is new in Python 2.6; fall back on simplejson if it
# isn't available.
try:
    import json
except:
    import simplejson as json
    pass

# JIVE Python modules.
from vex import Vex

def vex2time(str):
    tupletime = time.strptime(str, "%Yy%jd%Hh%Mm%Ss");
    return time.mktime(tupletime)

def time2vex(secs):
    tupletime = time.gmtime(secs)
    return time.strftime("%Yy%jd%Hh%Mm%Ss", tupletime)

usage = "usage: %prog [options] vexfile controlfiles"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-n", "--nodes", dest="number_nodes",
                  default=256, type="int",
                  help="Number of correlator nodes",
                  metavar="N")
parser.add_option("-d", "--no-delays", dest="gen_delays",
                  default = True,
                  action = "store_false",
                  help="Disable generating new delays")
parser.add_option("-m", "--machines", dest="machines",
                  default="a,b,c,d,e,f,g,h,i,j", type="string",
                  help="Machines to run correlator nodes on",
                  metavar="LIST")
(options, args) = parser.parse_args()

if len(args) < 2:
    parser.error("incorrect number of arguments")
    pass

vex_file = args[0]
# Parse the VEX file.
vex = Vex(vex_file)
exper = vex['GLOBAL']['EXPER']

# Proper time.
os.environ['TZ'] = "UTC"

mk5s = ['10.88.1.' + str(200+x) for x in range(17)]
mk5s += ['10.88.1.' + str(200+x) for x in range(20,26)]
manager_node = "head.sfxc"
output_node = "head.sfxc"
log_node = "head.sfxc"

# Generate a list of all media used for this experiment.
media = {}
for station in vex['STATION']:
    tapelog_obs = vex['STATION'][station]['TAPELOG_OBS']
    media[station] = []
    for vsn in vex['TAPELOG_OBS'][tapelog_obs].getall('VSN'):
        media[station].append({'vsn': vsn[1], 'start': vex2time(vsn[2]),
                               'stop': vex2time(vsn[3])})
        continue
    continue

# Generate a list of machines to use.
machines = []
for machine in options.machines.split(','):
    if machine in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']:
        for unit in [0, 1, 2, 3]:
            machines.append("sfxc-" + machine + str(unit) + ".sfxc")
            continue
        pass
    else:
        machines.append(machine)
        pass
    continue

for ctrl_file in args[1:]:
    basename = os.path.splitext(os.path.basename(ctrl_file))[0]
    machine_file = basename + ".machines"
    rank_file = basename + ".ranks"
    log_file = basename +".log"

    fp = open(ctrl_file, 'r')
    json_input = json.load(fp)
    fp.close()
    if 'file_parameters' not in json_input:
      json_input['file_parameters'] = {}
    start = vex2time(json_input['start'])

    stations = json_input['stations']
    stations.sort()

    # Make sure the delay files are up to date, and generate new ones
    # if they're not.
    if options.gen_delays:
      procs = {}
      success = True
      for station in stations:
          path = urlparse.urlparse(json_input['delay_directory']).path
          delay_file = path + '/' +  exper + '_' + station + '.del'
          if not os.access(delay_file, os.R_OK) or \
             os.stat(delay_file).st_mtime < os.stat(vex_file).st_mtime:
              args = ['generate_delay_model', vex_file, station, delay_file]
              procs[station] = subprocess.Popen(args, stdout=subprocess.PIPE)
              pass
          continue
      for station in procs:
          output = procs[station].communicate()[0]
          procs[station].wait()
          if procs[station].returncode != 0:
              print "Delay model couldn't be generated for " + station + ":"
              print output
              path = urlparse.urlparse(json_input['delay_directory']).path
              delay_file = path + '/' +  exper + '_' + station + '.del'
              os.remove(delay_file)
              success = False
              pass
          continue
      if not success:
          sys.exit(1)
          pass

    # Figure out the directory where the input data files are located.
    # Should be the same for all stations otherwise things will become
    # way too complicated.
    data_dir = None
    for station in stations:
      if station not in json_input['file_parameters']:
        data_source = json_input['data_sources'][station][0]
        if urlparse.urlparse(data_source).scheme == 'file':
            path = urlparse.urlparse(data_source).path
            if not data_dir:
                data_dir = os.path.dirname(path)
            assert data_dir == os.path.dirname(path)

    # Check if the input data files are there.  Do this in a loop that
    # gets repeated until all files have been found.
    missing = True
    while missing:
        if data_dir:
            # For every Mark5, generate a list of files present.  This is
            # faster than checking each file individually.
            data_dir_list = {}
            for mk5 in mk5s:
                args = ['/usr/bin/ssh', mk5, '/bin/ls ' + data_dir]
                p = subprocess.Popen(args, stdout=subprocess.PIPE)
                output = p.communicate()[0]
                p.wait()
                if p.returncode == 0:
                    data_dir_list[mk5] = output.split()
                else:
                    data_dir_list[mk5] = []
                    pass
                continue

            # For each station, figure out on which Mark5 the input data
            # is located.
            input_nodes = {}
            for station in stations:
                data_source = json_input['data_sources'][station][0]
                path = urlparse.urlparse(data_source).path
                file = os.path.basename(path)
                for mk5 in mk5s:
                    if file in data_dir_list[mk5]:
                        input_nodes[station] = mk5
                        break
                    continue
                continue
            pass
        else:
            # For every Mark5, generate a list of VSNs present.  This is
            # faster than checking each VSN individually.
            vsn_list = {}
            for mk5 in mk5s:
                args = ['/usr/bin/ssh', mk5, 'bin/vsnread']
                p = subprocess.Popen(args, stdout=subprocess.PIPE)
                output = p.communicate()[0]
                p.wait()
                if p.returncode == 0:
                    vsn_list[mk5] = output.split()
                else:
                    vsn_list[mk5] = []
                    pass
                continue

            # For each station, figure out on which Mark5 the input data
            # is located.
            input_nodes = {}
            for station in stations:
                data_source = json_input['data_sources'][station][0]
                path = urlparse.urlparse(data_source).path
                path = path.lstrip('/')
                vsn = path.split(':')[0]
                for mk5 in mk5s:
                    if vsn in vsn_list[mk5]:
                        input_nodes[station] = mk5
                        break

        # Check if we found them all.  If not, give the operator a
        # chance to mount the missing media.
        missing = False
        for station in stations:
          if station not in json_input['file_parameters']:
            if not station in input_nodes:
                if not missing:
                    print "Please mount media with the following VSNs:"
                    missing = True
                    pass
                for medium in media[station]:
                    if start >= medium['start'] and start < medium['stop']:
                        print medium['vsn']
                        pass
                    pass
                pass
            continue

        if missing:
            raw_input("Ready? ")
            pass
        continue

    # Initialize ranks
    ranks = {}
    # Create a MPI machine file for the job.
    fp = open(machine_file, 'w')
    print >>fp, manager_node, "slots=4" #Assume manager, output, and log node on the same machine
    #print >>fp, output_node
    #print >>fp, log_node
    for station in stations:
      if station not in json_input['file_parameters']:
        #ifhn = "ifhn="+input_nodes[station]
	print >>fp, " #", station
        print >>fp, input_nodes[station], " slots=4"
        ranks[input_nodes[station]]  = 4
    for station in json_input['file_parameters']:
      source = json_input['file_parameters'][station]['sources'][0]
      machine = source.partition(':')[0]
      if machine.startswith('sfxc-'):
        if not machine.endswith('.sfxc'):
          machine += '.sfxc'
          if machine not in machines:
            print >>fp, machine, " slots=8"
            ranks[machine] = 8
      elif machine.startswith('aribox') or (machine == '10.88.1.230'):
        print >>fp, machine, " slots=16"
        ranks[machine] = 16

    #for i in range(8):
    for machine in machines:
        print >>fp, machine, " slots=8"
        ranks[machine] = 8
    fp.close()

    # Create a MPI rank file for the job.
    fp = open(rank_file, 'w')
    print >>fp, "rank 0=", manager_node, "slot=0"
    print >>fp, "rank 1=", output_node, "slot=1"
    print >>fp, "rank 2=", log_node, "slot=2,3"
    rank=2
    # Create ranks
    for station in stations:
      rank += 1
      if station in json_input['file_parameters']:
        source = json_input['file_parameters'][station]['sources'][0]
        node = source.partition(':')[0]
        if node.startswith('sfxc-') and not node.endswith('.sfxc'):
            node += '.sfxc'
        nthread = json_input['file_parameters'][station]['nthreads']
        slots = 'slot=' + ','.join([str(i) for i in range(ranks[node]-nthread, ranks[node])])
        ranks[node] -= nthread
        print >>fp, "rank", str(rank), "=", node, slots
      else:
        print >>fp, "rank", str(rank), "=", input_nodes[station], "slot=0,1,2"

    for i in range(8):
        for machine in machines:
          if ranks[machine] > 0:
            rank += 1
            print >>fp, "rank", str(rank), "=", machine, "slot=", str(ranks[machine]-1)
            ranks[machine] -= 1
    fp.close()

    # Start the job.
    number_nodes = 3 + len(stations) + options.number_nodes
    sfxc = "`which sfxc`"
    #cmd = "mpirun -machinefile " + machine_file + " " \
    #    + "-n " + str(number_nodes) + " " \
    #    + sfxc + " " + ctrl_file + " " + vex_file \
    #    + " 2>&1 | tee " + log_file
    cmd = "mpirun --mca btl_tcp_if_include bond0,eth2.4,eth0 " \
        + "--mca oob_tcp_if_exclude eth1,eth2,eth3 "\
        +"-machinefile " + machine_file + " " \
        "--rankfile " + rank_file + " " \
        + "-n " + str(number_nodes) + " " \
        + sfxc + " " + ctrl_file + " " + vex_file \
        + " 2>&1 | tee " + log_file
    print cmd
    os.system(cmd)

    continue
