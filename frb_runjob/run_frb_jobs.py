#! /usr/bin/env python3

# Standard Python modules.
from argparse import ArgumentParser
import json
import os
import re
import subprocess
from functools import partial
from multiprocessing import Pool
from shutil import which
from urllib.parse import urlsplit
from vex import Vex

# Number of cores in sfxc nodes
SFXC_CORES = {'out':20, 'k': 16, 'l': 16, 'm': 20, 'n': 20, 'o': 20}
FB_CORES = {'flexbuf0': 12, 'flexbuf1': 16, 'flexbuf2': 16, 'flexbuf3': 12,
            'flexbuf4': 16, 'flexbuf5': 16, 'flexbuf6': 16, 'flexbuf7': 20,
            'flexbuf8': 20, 'flexbuf9': 16, 'flexbuf10': 20, 'flexbuf11': 20,
            'flexbuf12': 16, 'flexbuf13': 20, 'flexbuf14': 20, 'flexbuf15': 20,
            'flexbuf16': 20, 'flexbuf17': 20, 'flexbuf18': 20, 'aribox': 8}

def generate_delay_file(vex_file, x):
    """ Create delay file, this is intended to be executed using Pool.map()
    """
    station, filename = x
    args = ['generate_delay_model', vex_file, station, filename]
    subprocess.run(args, check=True)

def run_job(vex_file, ctrl_file, generate_delays, machines):
    vex = Vex(vex_file)
    with open(ctrl_file, 'r') as f:
        ctrl = json.load(f)
    exper = ctrl['exper_name']

    basename = os.path.splitext(os.path.basename(ctrl_file))[0]
    machine_file = basename + ".machines"
    rank_file = basename + ".ranks"
    log_file = basename +".log"

    # Generate a list of machines to use.
    head_node = 'out'
    stations = ctrl['stations']
    inputmachines = [ctrl['flexbufs'][st] for st in stations]
    cormachines = []
    # FIXME: the machine name checking could be done with a more sophisticated regex
    machine_re = re.compile(r'([k-o])([0-9]{0,2})')
    for machine in machines:
        m = machine_re.match(machine)
        if m == None:
            print(f'Error: Unknown machine "{machine}"')
            exit(1)
        machine_frame, machine_nr = m.groups()
        first_nr, last_nr = (0,3) if machine_frame in ('k', 'l') else (1, 14)
        if len(machine_nr) == 0:
            for i in range(first_nr, last_nr+1):
                cormachines.append(f'sfxc-{machine}{i}')
        else:
            nr = int(machine_nr)
            if (nr < first_nr) or (nr > last_nr):
                print(f'Error: machine "{machine}" doesn\'t exist.')
                exit(1)
            cormachines.append(f'sfxc-{machine}')
    
    
    # Make sure the delay files exist or regenerate them if vexfile was updated
    x = urlsplit(ctrl['delay_directory'])
    delay_dir = x.path if x.path != '' else x.netloc
    delays_to_generate = []
    for station in stations:
        delay_file = os.path.join(delay_dir, f"{exper}_{station}.del")
        exists = os.path.isfile(delay_file)
        if exists:
            stat = os.stat(delay_file)
        non_zero_exists = (exists) and (stat.st_size > 16)
        if (not generate_delays) and (not non_zero_exists):
            print(f'Error: Delay file creation is disabled but "{delay_file}" doesn\'t exist or is too small')
            exit(1)
        # if vexfile is newer than delay file also generate delay
        if (not non_zero_exists) or (stat.st_mtime < os.stat(vex_file).st_mtime):
            delays_to_generate.append((station, delay_file))
    
    if (generate_delays) and (len(delays_to_generate) > 0):
        f = partial(generate_delay_file, vex_file)
        with Pool() as p:
            p.map(f, delays_to_generate)

    # Initialize ranks and create an openmpi machine file
    ranks = {}
    with open(machine_file, 'w') as fp:
        #Assume manager, output, and log node on the same machine
        ranks[head_node] = SFXC_CORES[head_node]
        fp.write(f'{head_node} slots={ranks[head_node]}\n')
        for machine in set(inputmachines):
            ranks[machine] = FB_CORES[machine]
            fp.write(f'{machine} slots={ranks[machine]}\n')
        ncorrelation = 0
        for machine in cormachines:
            machine_frame = machine.split('-')[1][0]
            ranks[machine] = SFXC_CORES[machine_frame]
            ncorrelation += ranks[machine]
            fp.write(f'{machine} slots={ranks[machine]}\n')
    
    # Create openmpi rank file
    with open(rank_file, 'w') as fp:
        fp.write(f'rank 0={head_node} slot=0\n') # Manager node
        fp.write(f'rank 1={head_node} slot=1\n') # Log node
        ranklist = [str(x) for x in range(2, SFXC_CORES[head_node])]
        fp.write(f'rank 2={head_node} slot={",".join(ranklist)}\n') # output node
        rank = 3
        for station, fb in zip(stations, inputmachines):
            fbcount = inputmachines.count(fb)
            ncores = min(4, FB_CORES[fb] // fbcount)
            ranks[fb] = max(0, ranks[fb] - ncores)
            ranklist = [str(x) for x in range(ranks[fb], ranks[fb]+ncores)]
            fp.write(f'rank {rank}={fb} slot={",".join(ranklist)} # {station}\n')
            rank += 1
        maxrank = rank + ncorrelation
        while rank < maxrank:
            for machine in cormachines:
                if ranks[machine] > 0:
                    ranks[machine] -= 1
                    fp.write(f'rank {rank}={machine} slot={ranks[machine]}\n')
                    rank += 1

    nprocess = rank
    sfxc = which('sfxc')
    args = ['mpirun', 
            '--mca', 'btl', '^openib', 
            '--mca', 'oob', '^ud', 
            '--mca', 'btl_tcp_if_include', '10.88.0.0/24', 
            '--mca', 'oob_tcp_if_include', '10.88.0.0/24', 
            '--mca', 'orte_keep_fqdn_hostnames', '1', 
            '--mca', 'orte_hetero_nodes', '1',
            '-machinefile', machine_file,
            '--rankfile', rank_file,
            '-n', str(nprocess),
            sfxc, ctrl_file, vex_file]
    print(args)

    with open(log_file, 'w') as fp:
        p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        for line in p.stdout:
            print(line, end='')
            fp.write(line)
        p.wait()

if __name__ == "__main__":
    p = ArgumentParser(description="Run SFXC job using ctrlfile made by create_frb_jobs.py")
    p.add_argument("-d", "--no-generate-delays",
                   default = False,
                   action = "store_true",
                   help="Don't create new delays even if vexfile was updated")
    p.add_argument("-m", "--machines",
                   default="k", type=str,
                   help="Machines to run correlator nodes on, default='k', allowed values are k,l,m,n,o ; it is also possible to list individual machines e.g. k1,k2")
    p.add_argument("vex", type=str, help="VEX file of experiment")
    p.add_argument("ctrlfiles", type=str, nargs='+', help="SFXC control files generated by create_frb_jobs.py")
    args = p.parse_args()
    machines = args.machines.split(',')
    for ctrlname in args.ctrlfiles:
        run_job(args.vex, ctrlname, not args.no_generate_delays, machines)