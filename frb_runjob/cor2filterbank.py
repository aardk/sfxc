#!/usr/bin/env python3
import sys
import struct
import pdb
import optparse
import datetime
import vex
import re
import signal
from numpy import *
from time import sleep
from threading import Thread, current_thread
from multiprocessing import Queue, Process
timeslice_header_size = 16
uvw_header_size = 32
stat_header_size = 24
baseline_header_size = 8
nskip = 0
NSLICE_PER_INT = 32  # Divide integrations in NSLICE_PER_INT chuncks


def get_configuration(vexfile, corfile, setup_station):
    # TODO: This should also determine which subbands have been correlated
    cfg = {}
    corfile.seek(0)
    global_header_size, = struct.unpack('i', corfile.read(4))
    cfg['global_header_size'] = global_header_size
    print(f'global = {global_header_size}')
    corfile.seek(0)
    gheader_buf = corfile.read(global_header_size)
    global_header = struct.unpack('i32s2h5ib15s', gheader_buf[:76])
    cfg["nchan"] = global_header[5]
    cfg["inttime"] = global_header[6] / 1000000.
    if global_header[9] < 2:
        cfg["npol"] = 1
    elif global_header[9] == 2:
        cfg["npol"] = 2
    else:
        cfg["npol"] = 4
    # FIXME This should include nskip
    cfg["start_time"] = get_time(global_header[2], global_header[3], global_header[4])

    tsheader_buf = corfile.read(timeslice_header_size)
    timeslice_header = struct.unpack('4i', tsheader_buf)
    nsubint = timeslice_header[1]
    if cfg["npol"] == 4:
        nsubint //= 4
    cfg["nsubint"] = nsubint
    print(global_header[6]//nsubint, nsubint)
    if nsubint * (global_header[6]//nsubint) != global_header[6]:
        print("Error: There should be an integer number of subintegrations per integration")
        sys.exit(1)

    nsubband, maxfreq, bw = get_freq(vexfile, cfg["start_time"], setup_station)
    cfg["nsubband"] = nsubband
    cfg["maxfreq"] = maxfreq
    cfg["bandwidth"] = bw
    src_name, src_raj, src_dec = get_source(vexfile, cfg["start_time"])
    cfg["src_name"] = src_name
    cfg["src_raj"] = src_raj
    cfg["src_dec"] = src_dec

    cfg["mjd"] = mjd(global_header[2], global_header[3], global_header[4])
    print("jday=", cfg["mjd"], ", year = ", global_header[2],
          ", day=", global_header[3], ", sec = ", global_header[4])
    return cfg


def write_header(cfg, outfile, bif, eif, npol_out, decimate):
    nchan = cfg["nchan"] // 2 if decimate else cfg["nchan"]
    bw = cfg["bandwidth"]
    nsubband = cfg["nsubband"]

    header = [["HEADER_START"]]
    hlen = struct.pack('i', len(sys.argv[1]))
    header.append(["rawdatafile", hlen, sys.argv[1]])

    h = cfg["src_name"]
    hlen = struct.pack('i', len(h))
    header.append(["source_name", hlen, h])
    h = struct.pack('i', 0)
    header.append(["machine_id", h])
    h = struct.pack('i', 0)
    header.append(["telescope_id", h])
    h = struct.pack('d', cfg["src_raj"])
    header.append(['src_raj', h])
    h = struct.pack('d', cfg["src_dec"])
    header.append(['src_dej', h])
    h = struct.pack('d', 0.0)
    header.append(['az_start', h])
    header.append(['za_start', h])
    data_type = 1  # 0: raw data, 1 : filterbank
    h = struct.pack('i', data_type)
    header.append(['data_type', h])
    foff = -bw / nchan
    h = struct.pack('d', foff)
    header.append(['foff', h])
    fch1 = cfg["maxfreq"] - (nsubband-eif-1) * bw + foff / 2.
    print(fch1, bif, eif, bw)
    h = struct.pack('d', fch1)
    header.append(['fch1', h])
    h = struct.pack('i', nchan * (eif-bif+1))
    header.append(['nchans', h])
    h = struct.pack('i', 1)
    header.append(['nbeams', h])
    header.append(['ibeam', h])
    h = struct.pack('i', 32)
    header.append(['nbits', h])
    h = struct.pack('d', cfg["mjd"])
    header.append(['tstart', h])
    h = struct.pack('d', cfg["inttime"] / cfg["nsubint"])
    header.append(['tsamp', h])
    h = struct.pack('i', npol_out)
    header.append(['nifs', h])
    header.append(["HEADER_END"])

    for h in header:
        hlen = struct.pack('i', len(h[0]))
        outfile.write(hlen)
        outfile.write(h[0].encode())
        print(h)
        for i in h[1:]:
            if type(i) == str:
                outfile.write(i.encode())
            else:
                outfile.write(i)


def get_source(vexfile, start_time):
    scan = get_scan(vexfile, start_time)
    source = vexfile['SCHED'][scan]['source']
    name = vexfile["SOURCE"][source]["source_name"]
    sra = vexfile["SOURCE"][source]["ra"]
    tra = re.split('h|m|s', sra)
    ra = (int(tra[0])*100 + int(tra[1]))*100. + float(tra[2])
    sdec = vexfile["SOURCE"][source]["dec"]
    tdec = re.split('d|\'|\"', sdec)
    dec = (int(tdec[0])*100 + int(tdec[1]))*100 + float(tdec[2])
    return name, ra, dec


def get_freq(vexfile, start_time, setup_station):
    scan = get_scan(vexfile, start_time)
    mode = vexfile['SCHED'][scan]['mode']
    if setup_station == '':
        freq = vexfile['MODE'][mode]['FREQ'][0]
    else:
        freq = ''
        for f in vexfile['MODE'][mode].getall('FREQ'):
            if setup_station in f[1:]:
                freq = f[0]
        if freq == '':
            raise KeyError(setup_station)
    channels = set()
    for channel in vexfile['FREQ'][freq].getall('chan_def'):
        f0 = float(channel[1].split()[0])
        sb = 0 if (channel[2].strip().upper() == 'L') else 1
        bw = float(channel[3].split()[0])
        print(f0 + sb*bw, bw)
        channels.add(f0 + sb*bw)
    nsubband = len(channels)
    maxfreq = max(channels)
    print('bw = ', bw)
    return nsubband, maxfreq, bw


def get_scan(vexfile, start_time):
    for scan in vexfile['SCHED'].items():
        t = scan[1]['start']
        t = [int(x) if x != '' else 0 for x in re.split('y|d|h|m|s', t)]
        scan_start = get_time(t[0], t[1], t[2]*3600+t[3]*60+t[4])
        scan_len = datetime.timedelta(seconds=int(
            scan[1]['station'][2].rstrip("sec")))
        scan_end = scan_start + scan_len
        if start_time < scan_end:
            return scan[0]
    print("Could not find scan for t = ", start_time)
    sys.exit(1)


def mjd(year, day_of_year, sec_of_day):
    y = int(year) + 4799
    jdn = 365*y + (y//4) - (y//100) + (y//400) - 31738 - 2400000.5

    jdn = int(jdn) + day_of_year - 1
    return jdn + sec_of_day / 86400.


def print_global_header(infile):
    infile.seek(0)
    global_header_size, = struct.unpack('i', infile.read(4))
    infile.seek(0)
    gheader_buf = infile.read(global_header_size)
    global_header = struct.unpack('i32s2h5i4b', gheader_buf[:64])
    hour = global_header[4] // (60*60)
    minute = (global_header[4] % (60*60))//60
    second = global_header[4] % 60
    print("Experiment %s, SFXC version = %s, date = %dy%dd%dh%dm%ds, nchan = %d, int_time = %d, pol = %s" % (
        global_header[1].decode(), global_header[8], global_header[2], global_header[3], hour, minute, second, global_header[5], global_header[6], int(global_header[9])))


def parse_args():
    usage = "Usage : %prog [OPTIONS] <vex file> <cor file 1> ... <cor file N> <output_file>"
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-i", "--ifs", dest='ifs', type='string',
                      help='Range of sub-bands to correlate, format first:last. ' +
                      'For example -i 0:3 will write the first 4 sub-bands.')
    parser.add_option("-s", "--setup-station", dest='setup_station', type='string',
                      help="Define setup station, the frequency setup of this station is used in the conversion. " +
                           "The default is to use the first station in the vexfile",
                      default="")
    parser.add_option("-b", "--bandpass", dest='bandpass', default=False, action="store_true",
                      help='Apply bandpass')
    parser.add_option("-z", "--zerodm", dest='zerodm', default=False, action="store_true",
                      help='Apply zerodm subtraction')
    parser.add_option("-d", "--decimate-freq", dest='decimate', default=False, action="store_true",
                      help='Compensate for zeropadding by removing the odd numbered frequency points, default=no')
    parser.add_option("-p", "--pol", dest='pol', type='string', default='I',
                      help='Which polarization to use: R, L, I (=R+L), or F(=R, Re(RL), Im(RL), L), default=I')
    (opts, args) = parser.parse_args()

    if opts.ifs == None:
        bif = 0
        eif = -1
    else:
        ifs = opts.ifs.partition(':')
        bif = int(ifs[0])
        if ifs[1] != '':
            eif = int(ifs[2])
        else:
            eif = bif

    if opts.pol.upper() == 'L':
        pol = 2
    elif opts.pol.upper() == 'R':
        pol = 1
    elif opts.pol.upper() == 'I':
        pol = 3
    else:
        pol = 4

    infiles = []
    nargs = len(args)
    if nargs < 3:
        parser.error('Invalid number of arguments')
    try:
        for i in range(1, nargs - 1):
            infiles.append(open(args[i], 'rb'))
    except:
        print("Could not open file : " + args[i])
        sys.exit()

    vexfile = vex.Vex(args[0])
    outfile = open(args[-1], 'wb')
    return vexfile, infiles, outfile, bif, eif, pol, opts.setup_station, opts.zerodm, opts.bandpass, opts.decimate


def get_time(year, day, seconds):
    t = datetime.datetime(year, 1, 1)
    t += datetime.timedelta(days=day-1)
    t += datetime.timedelta(seconds=seconds)
    return t


def read_integration(infile, cfg):
    nchan = cfg['nchan']
    nsubband = cfg['nsubband']
    nsubint = cfg['nsubint']
    npol = cfg['npol']
    results = []
    tsheader_buf = infile.read(timeslice_header_size)
    if len(tsheader_buf) != timeslice_header_size:
        return results
    # get timeslice header
    timeslice_header = struct.unpack('4i', tsheader_buf)
    current_slice = timeslice_header[0]
    channel = 0
    nread = 0
    while (len(tsheader_buf) == timeslice_header_size) and (channel < nsubband*npol):
        # Read UVW
        nuvw = timeslice_header[2]
        infile.read(uvw_header_size * nuvw)
        # Read the bit statistics
        nstatistics = timeslice_header[3]
        infile.read(stat_header_size * nstatistics)
        # Read the baseline data
        nbaseline = timeslice_header[1]
        slice_size = nbaseline * (baseline_header_size + (nchan + 1) * 4)
        data = infile.read(slice_size)
        if len(data) != slice_size:
            return []
        results.append((timeslice_header, data))
        # Get next time slice header
        if nbaseline > nsubint:
            channel += 4  # If data has cross-polls then all polarization in each timeslice
        else:
            channel += 1

        if channel < nsubband*npol:
            tsheader_buf = infile.read(timeslice_header_size)
            if len(tsheader_buf) == timeslice_header_size:
                timeslice_header = struct.unpack('4i', tsheader_buf)
    return results


def parse_integration(indata, cfg, polarization):
    nchan = cfg['nchan']
    nsubband = cfg['nsubband']
    nsubint = cfg['nsubint']
    npol_in = 1 if (cfg['npol'] < 4) else 4
    npol_out = 1 if (polarization < 4) else 4

    data = zeros([nsubint, npol_out, nsubband * nchan], dtype=float32)
    for time_slice_header, time_slice in indata:
        for ipol in range(npol_in):
            size_slice = nsubint * (baseline_header_size + (nchan + 1) * 4)
            if len(time_slice) != size_slice * npol_in:
                break

            index = ipol * size_slice
            for b in range(nsubint):
                bheader = struct.unpack('i4B', time_slice[index:index+8])
                index += baseline_header_size
                baseline = frombuffer(time_slice, count=(
                    nchan+1), offset=index, dtype='float32')
                index += (nchan + 1) * 4
                weight, station1, station2, byte = bheader[:4]
                pol1 = byte & 1
                pol2 = (byte >> 1) & 1
                sideband = (byte >> 2) & 1
                freq_nr = byte >> 3
                channel_nr = 2*freq_nr + sideband
                outpol = -1
                if (pol1 == pol2) and ((pol1 + 1) & polarization):
                    # --pol = L, R or I
                    outpol = 0
                elif polarization == 4:
                    # --pol = F
                    outpol = abs(3 * pol1 - 2 * pol2)

                # print integration, bheader
                if (outpol >= 0):
                    if sideband == 0:
                        # outpol=3 is Im(RL), for LSB we need the complex conjugate
                        if outpol == 3:
                            vreal = -baseline[1:(nchan+1)]
                        else:
                            vreal = baseline[1:(nchan+1)]
                    else:
                        vreal = baseline[nchan-1::-1]
                    if isnan(vreal).any() == False:
                        # We write data in order of decreasing frequency
                        inv_ch = nsubband - channel_nr - 1
                        data[b, outpol, (inv_ch*nchan)
                                         :((inv_ch+1)*nchan)] += vreal
                    else:
                        print("b=("+repr(station1)+", "+repr(station2)+"), freq_nr = " +
                              repr(freq_nr)+",sb="+repr(sideband)+",pol="+repr(pol))
                        print("invalid data (not a number)")
    return data


def pad_zeros(outfile, npad, nsubint, nchan, npol):
    pad = zeros([nsubint, npol, nchan*nsubband], dtype=float32)
    for i in range(npad):
        pad.tofile(outfile)


def get_bandpass(infile, cfg, polarization):
    # First determine the size of one subint
    # FIXME: This is broken for cross-polls
    infile.seek(cfg['global_header_size'])
    time_slice_buffer = read_integration(infile, cfg)
    data = parse_integration(time_slice_buffer, cfg, polarization)
    pos = infile.tell()
    size = pos - cfg['global_header_size']
    infile.seek(0, 2)
    endpos = infile.tell()

    # Use 10 seconds in the middle of the scan
    n = endpos // size
    inttime = cfg["inttime"]
    if (n <= 3) or (n*inttime <= 3):
        start = 0
        toread = n
    elif (n <= 5) or (n*inttime <= 5):
        p = int(ceil(1./inttime))
        start = p
        toread = n - p - start
    else:
        m = int(ceil(5./inttime))
        start = n//2 - m
        toread = n//2 + m - start

    # Make the bandpass
    infile.seek(cfg['global_header_size'] + start*size)
    time_slice_buffer = read_integration(infile, cfg)
    data = parse_integration(time_slice_buffer, cfg, polarization)
    bandpass = data.sum(axis=0)
    for i in range(1, toread):
        time_slice_buffer = read_integration(infile, cfg)
        data = parse_integration(time_slice_buffer, cfg, polarization)
        bandpass += data.sum(axis=0)

    # Now normalize the bandpass per IF using median and smoothing RFI lines
    nsubband = cfg["nsubband"]
    nchan = cfg["nchan"]
    for i in range(nsubband):
        guard = int(round(nchan/10.))
        inner = bandpass[(i*nchan + guard):((i+1)*nchan - guard)]
        p = poly1d(polyfit(arange(len(inner)), inner, 8, w=1./(inner + 1e-6)))
        inner[:] = p(arange(len(inner)))
    maxnorm = median(bandpass)
    bandpass /= max(1., maxnorm)

    return bandpass


def start_next_input(infile, cfg, npol_out, decimate, nwritten, total_written, old_start_time):
    decimate_frac = 2 if decimate else 1
    infile.seek(0)
    global_header_size = struct.unpack('i', infile.read(4))[0]
    infile.seek(0)
    gheader_buf = infile.read(global_header_size)
    global_header = struct.unpack('i32s2h5i4c', gheader_buf[:64])
    nchan = global_header[5]
    inttime = global_header[6] // 1000000
    start_time = get_time(global_header[2], global_header[3], global_header[4])
    tsheader_buf = infile.read(timeslice_header_size)
    timeslice_header = struct.unpack('4i', tsheader_buf)
    nsubint = timeslice_header[1]
    if total_written > 0:
        # Check input and pad gap beteen scans with zeros
        error = False
        if cfg["nchan"] != nchan:
            print('Error: number of channels not constant between files')
            error = True
        if cfg["nsubint"] != nsubint:
            print('Error: number of subints per integration differs between files')
            error = True
        if cfg["inttime"] != inttime:
            print('Error : interation time differs between files')
            error = True
        if start_time <= old_start_time:
            print('Error : Inout files not in ascending time order')
            error = True
        else:
            dt = start_time - old_start_time
            diff = dt.days*86400 + dt.seconds
            if (diff % inttime) != 0:
                print(
                    'Error: consequtive input files have to be an integer number of integration times appart')
                error = True
        if error:
            print('Current file :', infile.name)
            sys.exit(1)
        npad = diff / inttime - nwritten + nskip
        print('padding ', npad, 'integrations')
        print('diff, inttime,nwritten,nskip = ',
              diff, inttime, nwritten, nskip)
        pad_zeros(outfile, npad, nsubint, nchan // decimate_frac, npol_out)
    else:
        npad = 0
    return start_time, npad


def reader_thread(outqueue, infile, cfg):
    thread = current_thread()
    thread.is_running = True
    timeslice_buffer = read_integration(infile, cfg)
    i = 0
    n = len(outqueue)
    while (len(timeslice_buffer) > 0) and (thread.is_running):
        while not outqueue[i].empty():
            sleep(0.1)
        for buf in timeslice_buffer:
            outqueue[i].put(buf)
        timeslice_buffer = read_integration(infile, cfg)
        i = (i+1) % n
    for i in range(n):
        outqueue[i].put(tuple())


def worker_thread(inqueue, outqueue, cfg, polarization, zerodm=False, bandpass=[]):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    nchan = cfg["nchan"]
    nsubband = cfg['nsubband']
    nsubint = cfg['nsubint']
    npol = cfg['npol']

    while(True):
        indata = []
        ch = 0
        while ch < nsubband * npol:
            d = inqueue.get()
            if len(d) == 0:
                outqueue.put([])
                return
            indata.append(d)
            nbaseline = d[0][1]
            if nbaseline > nsubint:
                ch += 4  # If data has cross-polls then all polarization in each timeslice
            else:
                ch += 1

        data = parse_integration(indata, cfg, polarization)
        if len(bandpass) > 0:
            # Don't blow up band edges too much
            bp = [x if x > 1e-2 else 1. for x in bandpass]
            data /= bp
        if zerodm:
            for i in range(start//nchan, end//nchan):
                cdata = data[:, i*nchan:(i+1)*nchan]
                for j in range(data.shape[0]):
                    cdata[j] -= cdata[j].sum() // nchan
        m = nsubint // NSLICE_PER_INT
        for i in range(NSLICE_PER_INT - 1):
            outqueue.put(data[i*m:(i+1)*m])
        outqueue.put(data[(NSLICE_PER_INT-1) * m:])


def writer_thread(inqueue, outfile, cfg, bif, eif, decimate_frac):
    nchan = cfg["nchan"]
    # NB: we ordered subbands in reverse order
    start = (cfg["nsubband"] - eif - 1) * nchan
    end = (cfg["nsubband"] - bif) * nchan
    slicenr = 0
    i = 0
    n = len(inqueue)
    data = inqueue[slicenr].get()
    while len(data) > 0:
        if len(data.shape) == 3:
            data[:, :, start:end:decimate_frac].tofile(outfile)
        else:
            data[:, start:end:decimate_frac].tofile(outfile)
        i += 1
        if i == NSLICE_PER_INT:
            i = 0
            slicenr += 1
            print("Written time slice", slicenr)
        data = inqueue[slicenr % n].get()


def create_filterbank(vexfile, infiles, outfile, bif, eif, pol, setup_station, zerodm=False, dobp=False, decimate=False):
    nthreads = 8  # For now hardcode the number of converter threads

    # Read global header
    global_header_size = struct.unpack('i', infiles[0].read(4))[0]
    cfg = get_configuration(vexfile, infiles[0], setup_station)
    if eif == -1:
        eif = cfg["nsubband"] - 1
    start_time = cfg["start_time"]
    
    if (pol < 4):
        npol_out = 1
    else:
        npol_out = 4
        if cfg["npol"] != 4:
            print(
                "Error: Full polarization requested, but correlator output doesn't contain cross-polls")
            exit(1)
    # Write filterbank header
    write_header(cfg, outfile, bif, eif, npol_out, decimate)
    start_time = cfg["start_time"]
    nwritten = 0
    total_written = 0
    for infile in infiles:
        print_global_header(infile)
        # Check next input and pad zeros if necessary
        start_time, npad = start_next_input(
            infile, cfg, npol_out, decimate, nwritten, total_written, start_time)
        total_written += npad
        if zerodm or dobp:
            bandpass = get_bandpass(infile, cfg, pol)
        else:
            bandpass = []
        nwritten = 0
        nchan = cfg["nchan"]
        infile.seek(global_header_size)
        nsubslice = cfg["nsubband"] * cfg["npol"]
        # Start reader thread
        inqueue = [Queue(nsubslice) for i in range(nthreads)]
        reader = Thread(target=reader_thread, args=(inqueue, infile, cfg))
        reader.daemon = True
        reader.start()
        # Start output thread
        outqueue = [Queue() for i in range(nthreads)]
        decimate_frac = 2 if decimate else 1
        writer = Thread(target=writer_thread, args=(
            outqueue, outfile, cfg, bif, eif, decimate_frac))
        writer.start()
        # Worker processes (Multiprocessing module)
        workers = [Process(target=worker_thread, args=(
            inqueue[i], outqueue[i], cfg, pol, zerodm, bandpass)) for i in range(nthreads)]
        for worker in workers:
            worker.start()
        # Wait for threads to close using polling in order to keep ctrl-c working
        try:
            while writer.is_alive():
                writer.join(0.5)
        except KeyboardInterrupt:
            reader.is_running = False
        writer.join()
        for worker in workers:
            worker.join()
        reader.join()


############################### Main program ##########################
#########
if __name__ == "__main__":
    vexfile, infiles, outfile, bif, eif, pol, setup_station, zerodm, dobp, decimate = parse_args()
    create_filterbank(vexfile, infiles, outfile, bif, eif,
                      pol, setup_station, zerodm, dobp, decimate)
