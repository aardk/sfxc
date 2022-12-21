#!/usr/bin/env python3
import json
import os
import re
import struct
from argparse import ArgumentParser
from collections import namedtuple
from datetime import datetime, timedelta
from dis import dis
from multiprocessing.pool import ThreadPool
from stat import S_ISDIR

import paramiko
import pylab as pl

from vex import Vex

FLEXBUFS = ['aribox'] + [f"flexbuf{i}" for i in range(0,19)]
FB_BASE_PATH = '/mnt/'
VBSFILE = namedtuple('VBSFile', ['filename', 'index'])
BURST = namedtuple('Burst', ['scan', 'mjd', 'bursttime', 'topostation', 'dm'])

def mjd2date(mjdtime:float):
    t0 = datetime(1858, 11, 17)
    days = pl.floor(mjdtime)
    fsec = (mjdtime - days) * 86400
    sec = int(pl.floor(fsec))
    us = int(1e6 * (fsec - sec))
    return t0 + timedelta(days, sec, us)

def vextime(tstring):
    year, doy, hour, minute, second = (int(x) for x in re.split('y|d|h|m|s', tstring)[:-1])
    nsec = 60 * (60 * hour + minute) + second
    t0 = datetime(year, 1, 1)
    return t0 + timedelta(days=doy-1, seconds=nsec)

def datetime2vex(dt, integer_sec = True):
    doy = (dt - datetime(dt.year, 1, 1)).days + 1
    if integer_sec:
        return f"{dt.year}y{doy:03}d{dt.hour:02}h{dt.minute:02}m{dt.second:02}s"
    else:
        return f"{dt.year}y{doy:03}d{dt.hour:02}h{dt.minute:02}m{dt.second:02}.{dt.microsecond:06}s"


def get_scan_name(vex, bursttime, mjdtime):
    for scanname, scan in vex['SCHED'].items():
        tstart = vextime(scan['start'])
        maxsec = 0
        for station in scan.getall('station'):
            sec = int(station[2].split()[0])
            maxsec = max(maxsec, sec)
        tend = tstart + timedelta(seconds=sec)
        if (bursttime >= tstart) and (bursttime < tend):
            return scanname
    # More graceful exit would be nice
    print(f"Error: burst time {bursttime} (mjd={mjdtime}) doesn't match any scans in the vex file.")
    exit(1)

def create_burst_list(exper, vex, ctrl):
    bursts = {}
    with open(ctrl['burstfile'], 'r') as f:
        lines = f.readlines()
        for linenr, line in enumerate(lines):
            items = line.partition('#')[0].split()
            if len(items) == 0:
                continue
            if len(items) == 3:
                mjdtime = float(items[0])
                topostation = items[1]
                dm = float(items[2])
                bursttime = mjd2date(mjdtime)
                scan = get_scan_name(vex, bursttime, mjdtime)
                burst = BURST(scan, mjdtime, bursttime, topostation, dm)
                try:
                    bursts[scan].append(burst)
                except KeyError:
                    bursts[scan] = [burst]
            else:
                print(f'Error reading at line {linenr} of burstfile {ctrl["burstfile"]}\nContents of bad line: "{line}"')
                exit(1)
    return bursts

def vdif_start_time(sftp, filename):
    with sftp.open(filename) as f:
        h = struct.unpack('4I', f.read(16))
        sec_from_epoch = h[0] & 0x3fffffff
        ref_epoch = (h[1] >> 24) & 0x3f
        ref_year = 2000 + ref_epoch // 2
        ref_month = 1 + 6 * (ref_epoch & 1)
        t = datetime(ref_year, ref_month, 1) + timedelta(seconds=sec_from_epoch)
        return t
    # FIXME fail more gracefully here
    print(f'Error reading {filename} over sftp')
    exit(1)

def get_file_list(exper_:str):
    exper = exper_.lower()
    pool = ThreadPool()
    f = lambda x: get_fb_list(exper, x)
    results = pool.map(f, FLEXBUFS)
    flexbufs = {fb: results[i] for i, fb in enumerate(FLEXBUFS) if len(results[i]) > 0}
    return flexbufs

def get_fb_list(exper:str, fb:str):
    number_suffix = re.compile(r'\w+_\w+_\w+\.([0-9]+)')
    vbs_dir = re.compile(r'([a-zA-Z0-9]+)_([a-zA-Z0-9][a-zA-Z0-9])_(\w+)')
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(fb)
    sftp = ssh.open_sftp()
    # Loop over disks
    vbsscans = {}
    disks = [x for x in sftp.listdir_attr(FB_BASE_PATH) if S_ISDIR(x.st_mode)]
    for disk in disks:
        diskname = os.path.join(FB_BASE_PATH, disk.filename)
        dirs = [x for x in sftp.listdir_attr(diskname) if S_ISDIR(x.st_mode)]
        for d in dirs:
            m = vbs_dir.match(d.filename)
            if m != None:
                exp, station, scan = (g.lower() for g in m.groups())

                if exp == exper:
                    expdir = os.path.join(diskname, d.filename)
                    files =  [VBSFILE(os.path.join(expdir, x.filename), 
                                        int(number_suffix.match(x.filename).groups()[0])) 
                                for x in sftp.listdir_attr(expdir)
                                if (x.st_size > 32) and (number_suffix.match(x.filename) != None)]
                    if len(files) > 0:
                        if station not in vbsscans:
                            vbsscans[station] = {}
                        if scan not in vbsscans[station]:
                            vbsscans[station][scan] = {'flexbuf': fb, 'files': files}
                        else:
                            vbsscans[station][scan]['files'] += files
    # create sorted list of files for each station / scan
    for station, scans in vbsscans.items():
        for scanname, item in scans.items():
            item['files']  = sorted(item['files'], key = lambda x: x.index)
            item['tfirst'] = vdif_start_time(sftp, item['files'][0].filename)
            item['tlast'] = vdif_start_time(sftp, item['files'][-1].filename)
    sftp.close()
    ssh.close()
    return vbsscans

def create_polyco(polycofile, source_name, fmid, burst):
     with open(polycofile, 'w') as f:
        date = burst.bursttime.strftime('%d-%b-%y')
        utc = '000000.00' # Start midnight because it can be exactly represented in floating point
        tmid = pl.floor(burst.mjd)
        f.write(f"{source_name[:9]:9} {date:9} {utc:11} {tmid:20} {burst.dm:20}  0.000  0.000\n")
        f.write(f"0.000000            0.500000000000    COE  1600 3    {fmid:.2f}\n")
        f.write(f"0.00000000000000000e+00  0.00000000000000000e+00  0.00000000000000000e+00\n")

def create_file_list(vbsscans, burst):
    tguard = timedelta(seconds=2) # how much time around burst to include
    for name, vbsscan in vbsscans.items():
        if (burst.bursttime >= vbsscan['tfirst']) and (burst.bursttime < vbsscan['tlast']):
            files = vbsscan['files']
            n_fragment = files[-1].index - files[0].index + 1
            if n_fragment > 1:
                dt = (vbsscan['tlast'] - vbsscan['tfirst']) / (n_fragment - 1)
            else:
                dt = timedelta(seconds=2) # if there is only one fragment the value doesn't matter
            first = max(int(pl.floor((burst.bursttime - vbsscan['tfirst'] - tguard) / dt)) + files[0].index, 
                    files[0].index)
            last = min(int(pl.ceil((burst.bursttime - vbsscan['tfirst'] + tguard) / dt)) + files[0].index, 
                    files[-1].index)

            return ["file://" + f.filename for f in files if (f.index >= first) and (f.index <= last)]
    return []

def create_ctrl_files(vex, ctrl, bursts, flexbufs):
    try:
        number_channels = ctrl["number_channels"]
    except KeyError:
        number_channels = 64
        print(f"Info: number_channels not defined, defaulting to {number_channels}")

    try:
        integr_time = ctrl["integr_time"]
    except KeyError:
        integr_time = 2.0
        print(f"Info: integr_time not defined, defaulting to {integr_time}")

    try:
        sub_integr_time = ctrl["sub_integr_time"]
    except KeyError:
        sub_integr_time = 64
        print(f"Info: sub_integr_time (for pulse search) not defined, defaulting to {sub_integr_time} [µs]")

    for scan, burstlist in bursts.items():
        try:
            stations = ctrl['stations']
        except KeyError:
            print("Info: No station list defined in ctrl file, using stations from VEX file.")
            stations = [station[0] for station in vex['SCHED'][scan].getall('station')]
    
        try:
            setup_station = ctrl['setup_station']
        except KeyError:
            setup_station = 'Ef' if 'Ef' in stations else stations[0]
            print(f"Info: No setup station defined in ctrl file defaulting to {setup_station}")

        mode = vex['SCHED'][scan]['mode']
        try:
            freq_id, = (freq[0] for freq in vex['MODE'][mode].getall('FREQ') if setup_station in freq[1:])
        except ValueError:
            print(f"Error: could'n find $FREQ definition of station {setup_station} in mode {mode}")
            exit(1)
        
        try:
            channels = ctrl['channels']
        except KeyError:
            channels = [line[4] for line in vex['FREQ'][freq_id].getall('chan_def')]
        
        fmid = max([float(line[1].split()[0]) + float(line[3].split()[0]) * (-1 if line[2] == 'L' else 1) / 2
                    for line in vex['FREQ'][freq_id].getall('chan_def') 
                    if line[4] in channels])
        scanstart = vextime(vex['SCHED'][scan]['start'])
        scanduration = max([int(x[2].split()[0]) for x in vex['SCHED'][scan].getall('station') if x[0] in stations])
        scanend = scanstart + timedelta(seconds=scanduration)

        for burstnr, burst in enumerate(burstlist):
            if burst.topostation in stations:
                search_station = burst.topostation
            else:    
                search_station = setup_station
                print(f"Warning: Burst was defined for station {burst.topostation} but it is not in the list of stations, using {search_station} for burst search instead")

            data_sources = {}
            available_stations = []
            fbmap = {}

            for station_ in stations:
                station  = station_.lower()
                # get list of all flexbufs that have data for this station
                fbs = [fb for fb, item in flexbufs.items() if station in item.keys()]
                # get list of files covering the burst for station
                for fb in fbs:
                    filelist = create_file_list(flexbufs[fb][station], burst)
                    if len(filelist) > 0:
                        fbmap[station_] = fb
                        break
                if (len(fbs) == 0) or (len(filelist) == 0):
                    print (f"Warning: no files found for station {station_}, dropping station from correlation!")
                else:
                    available_stations.append(station_)
                    data_sources[station_] = filelist

            outctrl = {}
            outctrl["channels"] = channels
            outctrl["window_function"] = "NONE"
            outctrl["cross_polarize"] = True
            outctrl["stations"] = available_stations
            outctrl["setup_station"] = setup_station
            outctrl["data_sources"] = data_sources
            outctrl["flexbufs"] = fbmap
            exper_block = vex["GLOBAL"]["EXPER"]
            outctrl["exper_name"] = vex["EXPER"][exper_block]["exper_name"]
            outctrl["integr_time"] = integr_time
            outctrl["message_level"] = 1
            outctrl["number_channels"] = number_channels
            outctrl["delay_directory"] = "file://delays"
            # We'll use 2.048s integrations for the pulse search
            dt = max(3, int(pl.ceil(integr_time)))
            nsec = timedelta(seconds=dt)
            if burst.bursttime + nsec < scanend:
                start_time = datetime2vex(burst.bursttime)
                stop_time = datetime2vex(burst.bursttime + nsec)
            else:
                tstart = burst.bursttime - timedelta(seconds=int(pl.floor(integr_time)))
                start_time = datetime2vex(tstart)
                stop_time = datetime2vex(tstart + nsec)
            outctrl["start"] = start_time
            outctrl["stop"] = stop_time
            basefilename = f"{outctrl['exper_name'].lower()}_{scan}_burst{burstnr}"
            outctrl['output_file'] = "file://" + basefilename + '.cor'
            polycofile = basefilename + '.polyco'
            ctrlfile = basefilename + '.ctrl'
            pulsars = {}
            sourcename = vex['SCHED'][scan]['source']
            create_polyco(polycofile, sourcename, fmid, burst)
            pulsars[sourcename] = {"nbins": 1, 
                                   "interval": [0,1], 
                                   "coherent_dedispersion": True,
                                   "polyco_file": "file://" + polycofile}
            outctrl['burst_mjd'] = burst.mjd
            outctrl['pulsars'] = pulsars
            outctrl['filterbank'] = False
            outctrl['pulsar_binning'] = True
            searchctrlfile = 'fb_' + basefilename + '.ctrl'
            outctrl['search_ctrl_file'] = "file://" + searchctrlfile
            with open(ctrlfile, 'w') as f:
                json.dump(outctrl, f, indent=4)
            print(outctrl)
            # Now create ctrl file used for pulse search
            outctrl['output_file'] = "file://fb_" + basefilename + '.cor'
            outctrl['integr_time'] = 2.048
            outctrl['sub_integr_time'] = sub_integr_time
            outctrl['stations'] = [search_station]
            outctrl['filterbank'] = True
            outctrl['pulsar_binning'] = False
            with open(searchctrlfile, 'w') as f:
                json.dump(outctrl, f, indent=4)

def main(vex, ctrl):
    exper = vex['GLOBAL']['EXPER']
    flexbufs = get_file_list(exper)
    bursts = create_burst_list(exper, vex, ctrl)
    create_ctrl_files(vex, ctrl, bursts, flexbufs)
    # Make sure delay directory exists
    os.makedirs("delays", exist_ok=True)

if __name__ == "__main__":
    p = ArgumentParser(description="Create jobs for FRB bursts")
    p.add_argument("vex", type=str, help="VEX file of experiment")
    p.add_argument("ctrl", type=str, help="JSON control file containing the bursts")
    args = p.parse_args()
    vex = Vex(args.vex)
    with open(args.ctrl, 'r') as f:
        ctrl = json.load(f)
    main(vex, ctrl)
