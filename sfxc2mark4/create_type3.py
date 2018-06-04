#!/usr/bin/env python

import json
import math
import os
import os.path
import struct
import sys
import time
import urlparse

import numpy as np
import scipy as sp
from scipy import fftpack, signal

from vex import Vex

from stationmap import station_codes

os.environ['TZ'] = 'UTC'
time.tzset()

class ScanInfo:
    def __init__(self, vex, station, scan):
        self.frequencies = []
        self.channels = []
        self.max_bandwidth = 0
        self.upper = False
        self.lower = False

        self.station_id = str(station)
        self.station_name = str(station.upper())
        self.scan = scan
        
        if station in station_codes:
            self.station_code = station_codes[station]
        else:
            self.station_code = 'x'
            pass

        mode = vex['SCHED'][scan]['mode']
        freqs = vex['MODE'][mode].getall('FREQ')
        for freq in freqs:
            if station in freq[1:]:
                break
            continue
        bbcs = vex['MODE'][mode].getall('BBC')
        for bbc in bbcs:
            if station in bbc[1:]:
                break
            continue
        _ifs = vex['MODE'][mode].getall('IF')
        for _if in _ifs:
            if station in _if[1:]:
                break
            continue
        value = vex['FREQ'][freq[0]]['sample_rate'].split()
        self.sample_rate = float(value[0])
        if value[1] == 'Gs/sec':
            self.sample_rate *= 1e9
        elif value[1] == 'Ms/sec':
            self.sample_rate *= 1e6
            pass
        channels = vex['FREQ'][freq[0]].getall('chan_def')
        for chan_def in channels:
            value = chan_def[1].split()
            frequency = float(value[0])
            if value[1] == 'GHz':
                frequency *= 1e9
            elif value[1] == 'MHz':
                frequency *= 1e6
            elif value[1] == 'KHz':
                frequency *= 1e3
                pass
            if not frequency in self.frequencies:
                self.frequencies.append(frequency)
                pass
            if chan_def[2] == 'U':
                self.upper = True
            else:
                self.lower = True
                pass
            value = chan_def[3].split()
            bandwidth = float(value[0])
            if value[1] == 'GHz':
                bandwidth *= 1e9
            elif value[1] == 'MHz':
                bandwidth *= 1e6
            elif value[1] == 'KHz':
                bandwidth *= 1e3
                pass
            if bandwidth > self.max_bandwidth:
                self.max_bandwidth = bandwidth
                pass
            continue
        self.frequencies.sort()

        for chan_def in channels:
            value = chan_def[1].split()
            frequency = float(value[0])
            if value[1] == 'GHz':
                frequency *= 1e9
            elif value[1] == 'MHz':
                frequency *= 1e6
            elif value[1] == 'KHz':
                frequency *= 1e3
                pass
            frequency = self.frequencies.index(frequency)
            if chan_def[2] == 'U':
                sideband = 1
            else:
                sideband = 0
                pass
            bbcs = vex['BBC'][bbc[0]].getall('BBC_assign')
            for bbc_assign in bbcs:
                if chan_def[5] == bbc_assign[0]:
                    break
                continue
            ifs = vex['IF'][_if[0]].getall('if_def')
            for if_def in ifs:
                if bbc_assign[2] == if_def[0]:
                    break
                continue
            if if_def[2] == 'R':
                pol = 0
            else:
                pol = 1
                pass
            self.channels.append((frequency, sideband, 0))
            continue

        stations = []
        for station_def in vex['STATION']:
            stations.append(station_def)
            continue
        stations.sort()
        self.station = stations.index(station)

        source = vex['SCHED'][scan]['source']
        self.source = vex['SOURCE'][source]['source_name']
        self.start = vex2time(vex['SCHED'][scan]['start'])
        self.length = 0
        for transfer in vex['SCHED'][scan].getall('station'):
            if transfer[0] == station:
                self.length = int(transfer[2].split()[0])
                pass
            continue

        clock = vex['STATION'][station]['CLOCK']
        clock_early = vex['CLOCK'][clock]['clock_early']
        self.clock = np.zeros(2)
        self.clock_epoch = vex2time(clock_early[2])
        value = clock_early[1].split()
        self.clock[0] = float(value[0])
        if value[1] == 'usec':
            self.clock[0] *= 1e-6
            pass
        value = clock_early[3].split()
        # If clock rate unit is missing, assume usec/sec.
        if len(value) == 1:
            value.append('usec/sec')
            pass
        self.clock[1] = float(value[0])
        if value[1] == 'usec/sec':
            self.clock[1] *= 1e-6
            pass
        return

    def chan_name(self, frequency, sideband, pol):
        sidebands = ['L', 'U']
        polarisations = ['R', 'L']

        band = get_band(self.frequencies[frequency])
        sideband = sidebands[sideband]
        pol = polarisations[pol]

        return "%c%02d%c%c" % (band, frequency, sideband, pol)

    pass

bands = [
    [ 'P', 225e6, 390e6 ],
    [ 'L', 390e6, 1750e6 ],
    [ 'S', 1750e6, 3900e6 ],
    [ 'C', 3900e6, 6200e6 ],
    [ 'X', 6200e6, 10900e6 ],
    [ 'K', 10900e6, 36000e6 ],
    [ 'Q', 36000e6, 46000e6 ],
]

def get_band(freq):
    for band in bands:
        if band[1] <= freq and freq <= band[2]:
            return band[0]
        continue
    return 'B'

ident = "3s2s3x16s32s8x"

type3 = "3s2s3x"
type300 = "!c2s32sx12sfH2x"
type301 = "!H32s6x6d"
type302 = "!H32s6x6d"
type303 = "!H32s6x144x18d"
type309_header = "!IIdd"
type309_channel = "!8sd128i"

# Convert UNIX time to a Mk4 date
def time2date(secs):
    tm = time.gmtime(secs)
    return struct.pack("!HHHHf", tm.tm_year, tm.tm_yday, tm.tm_hour,
                       tm.tm_min, tm.tm_sec)

def vex2time(str):
    tupletime = time.strptime(str, "%Yy%jd%Hh%Mm%Ss");
    return time.mktime(tupletime)

delay_header_v0 = "=I2sx"
delay_header_v1 = "=II2sx"
delay_scan = "=80sx"
delay_source = "=80sxI"
delay_entry = "=7d"

def parse_model(info, delay_file):
    fp = open(delay_file, 'r')
    buf = fp.read(struct.calcsize(delay_header_v0))
    hdr = struct.unpack(delay_header_v0, buf)
    if hdr[0] > 3:
        fp.seek(0)
        buf = fp.read(struct.calcsize(delay_header_v1))
        hdr = struct.unpack(delay_header_v1, buf)
        version = hdr[1]
    else:
        version = 0
        pass

    scan = info.scan

    if version > 0:
        buf = fp.read(struct.calcsize(delay_scan))
        scan = struct.unpack(delay_scan, buf)[0]
        pass
    buf = fp.read(struct.calcsize(delay_source))
    src = struct.unpack(delay_source, buf)
    source = src[0].strip()
    mjd = src[1]
    start = None
    t = []
    d = []
    u = []; v = [];  w = []
    while fp:
        buf = fp.read(struct.calcsize(delay_entry))
        if not buf:
            break
        delay = struct.unpack(delay_entry, buf)
        if delay[0] == 0 and delay[4] == 0:
            if (source == info.source and
                scan == info.scan and
                start >= info.start and
                start < info.start + info.length):
                return (t, d, u, v, w)
            if version > 0:
                buf = fp.read(struct.calcsize(delay_scan))
                if not buf:
                    break
                scan = struct.unpack(delay_scan, buf)[0]
                pass
            buf = fp.read(struct.calcsize(delay_source))
            if not buf:
                break
            src = struct.unpack(delay_source, buf)
            source = src[0].strip()
            mjd = src[1]
            start = None
            t = []
            d = []
            u = []; v = [];  w = []
            continue
        if start == None:
            start = (mjd - 40587) * 86400 + delay[0]
            pass
        t.append(delay[0])
        u.append(delay[1])
        v.append(delay[2])
        w.append(delay[3])
        d.append(delay[4])
        continue
    return

def create_splines(interval, x, y):
    splines = []
    diff_max = 0.0
    x = np.array(x)
    y = np.array(y)
    akima = sp.interpolate.Akima1DInterpolator(x, y)
    while len(x) > 1:
        points = min(interval + 1, len(x))
        u = np.linspace(0, points - 1, 7)
        v = np.array(map(lambda x: akima(x), x[0] + u))
        z = np.polyfit(u, v, 5)
        splines.append(np.flip(z, 0))
        poly = np.poly1d(z)
        a = np.arange(0, points - 1, 0.125)
        b = np.array(map(lambda x: poly(x), a))
        q = np.arange(x[0], x[0] + points - 1, 0.125)
        r = np.array(map(lambda x: akima(x), q))
	diff = np.max(np.absolute(r - b))
        if diff > diff_max:
            diff_max = diff
            pass
        x = x[points - 1:]
        y = y[points - 1:]
        continue
    return (splines, diff_max)

def write_type300(info, interval, out_fp):
    # Write Type300 header
    buf = struct.pack(type3, '300', '00')
    out_fp.write(buf)
    nsplines = (info.length + interval -1) // interval
    buf = struct.pack(type300, info.station_code, info.station_id,
                      info.station_name, time2date(info.start),
                      interval, nsplines)
    out_fp.write(buf)
    return

def write_type301(info, index, frequency, sideband, pol, spline, out_fp):
    # Write Type301 header
    buf = struct.pack(type3, '301', '00')
    out_fp.write(buf)
    chan_name = info.chan_name(frequency, sideband, pol)
    coefficients = spline
    buf = struct.pack(type301, index, chan_name, *coefficients)
    out_fp.write(buf)
    return

def write_type302(info, index, frequency, sideband, pol, spline, out_fp):
    # Write Type302 header
    buf = struct.pack(type3, '302', '00')
    out_fp.write(buf)
    chan_name = info.chan_name(frequency, sideband, pol)
    coefficients = spline * info.frequencies[frequency]
    buf = struct.pack(type302, index, chan_name, *coefficients)
    out_fp.write(buf)
    return

def write_type303(info, index, frequency, sideband, pol, u, v, w, out_fp):
    # Write Type303 header
    buf = struct.pack(type3, '303', '00')
    out_fp.write(buf)
    chan_name = info.chan_name(frequency, sideband, pol)
    coefficients = np.concatenate((u, v, w))
    buf = struct.pack(type303, index, chan_name, *coefficients)
    out_fp.write(buf)
    return

def write_type309(info, in_file, out_fp):
    # Calculate offset of first phasecal tone
    offset = [
        math.floor(info.frequencies[0] / 1e6) * 1e6 - info.frequencies[0],
        math.ceil(info.frequencies[0] / 1e6) * 1e6 - info.frequencies[0]
    ]

    # Calculate (relative) tone frequencies
    tone_frequencies = []
    if info.upper:
        freq = offset[1]
        while freq < info.max_bandwidth:
            tone_frequencies.append(freq)
            freq += 1e6
            continue
        pass
    if info.lower:
        freq = offset[0]
        while freq > -info.max_bandwidth:
            tone_frequencies.append(freq)
            freq -= 1e6
            continue
        pass
    ntones = len(tone_frequencies)
    if info.upper:
        upper = 0
        lower = ntones / 2
    else:
        upper = 0
        lower = 0
    pass

    in_fp = open(in_file, 'r')
    header = in_fp.read(68)

    data = {}
    times = []
    while in_fp:
        header = in_fp.read(20)
        if len(header) == 0:
            break
        header = struct.unpack("BBBBIIII", header)
        station = header[0]
        frequency = header[1]
        sideband = header[2]
        polarisation = header[3]
        mjd = header[4]
        secs = (mjd - 40587) * 86400 + header[5]
        sign = -1 if sideband == 0 else 1
        acc_period = float(header[6])
        num_samples = header[7]
        buf = in_fp.read(num_samples * 4)
        dd = np.frombuffer(buf, dtype='int32').astype(float)

        if not secs in data:
            data[secs] = {}
            pass
        if not station in data[secs]:
            data[secs][station] = {}
            pass
        key = (frequency, sideband, polarisation)
        if not secs in times:
            times.append(secs)

        # Filter out vectors with no valid data
        ddx = dd.max()
        if ddx == 0:
            continue

        dmvec = np.exp((sign * 2 * math.pi * 1j * offset[sideband] /
            info.sample_rate) * np.array(range(num_samples)))

        sdd = sp.fftpack.ifft(dd)
        sddh = np.concatenate((sdd[0:num_samples / 2], np.zeros(num_samples / 2)))
        dda = sp.fftpack.fft(sddh)

        ddc = dda * dmvec

        pp = np.sum(ddc.reshape((100, -1)), axis=0)

        spp = sp.fftpack.ifft(pp)

        data[secs][station][key] = spp[0:ntones/2]
        continue

    times.sort()

    station = info.station
    for secs in times:
        if not station in data[secs]:
            continue

        # Blast from the past: calculate ROT
        tm = time.gmtime(secs)
        rot = 32e6 * ((((tm.tm_yday - 1) * 24 + tm.tm_hour) * 60 +
                       tm.tm_min) * 60 + tm.tm_sec)

        # Write Type309 header
        buf = struct.pack(type3, '309', '01')
        out_fp.write(buf)
        buf = struct.pack(type309_header, 0, ntones, rot, acc_period)
        out_fp.write(buf)

        # Start with first tone frequency again
        tone = 0

        for key in data[secs][station]:
            frequency = key[0]
            sideband = key[1]
            pol = key[2]

            chan_name = info.chan_name(frequency, sideband, pol)

            sidebands = ['L', 'U']
            sideband = sidebands[sideband]

            if tone < ntones:
                freq = tone_frequencies[tone]
            else:
                freq = 0
                pass

            acc = np.zeros((64, 2))
            idx = 0
            norm = 128.0
            for entry in data[secs][station][key]:
                if sideband == 'U':
                    acc[upper + idx][0] = norm * entry.real
                    acc[upper + idx][1] = norm * entry.imag
                else:
                    acc[lower + idx][0] = -norm * entry.real
                    acc[lower + idx][1] = -norm * entry.imag
                    pass
                idx += 1
                continue

            buf = struct.pack(type309_channel, chan_name, freq, *acc.flatten())
            out_fp.write(buf)
            tone += 1
            continue
        out_fp.write('\0' * (64 - tone) * struct.calcsize(type309_channel))
        continue
    return

def process_job(vex, json_input, rootid, basename="1234", interval=120):
    exper_name = json_input['exper_name']

    phasecal_uri = json_input['phasecal_file']
    phasecal_file = urlparse.urlparse(phasecal_uri).path

    delay_uri = json_input['delay_directory']
    delay_directory = urlparse.urlparse(delay_uri).path

    scan = json_input['scans'][0]

    for station in json_input['stations']:

        info = ScanInfo(vex, station, scan)

        delay_file = exper_name + '_' + station + '.del'
        delay_file = os.path.join(delay_directory, delay_file)
        (t, d, u, v, w) = parse_model(info, delay_file)
        (d, diff_max) = create_splines(interval, t, d)
        (u, dummy) = create_splines(interval, t, u)
        (v, dummy) = create_splines(interval, t, v)
        (w, dummy) = create_splines(interval, t, w)
        splines = zip(d, u, v, w)

        filename = basename + "/" + scan + "/" + info.station_code + \
                   ".." + rootid
        fp = open(filename, 'w')
        buf = struct.pack(ident, '000', '01', "2001001-123456", str(filename))
        fp.write(buf)

        index = 0
        start = info.start
        write_type300(info, interval, fp)
        for spline in splines:
            d = spline[0]
            (u, v, w) = spline[1:]
            # Apply clock
            d[0] += info.clock[0] + (start - info.clock_epoch) * info.clock[1]
            d[1] += info.clock[1]
            for channel in info.channels:
                write_type301(info, index, channel[0], channel[1], channel[2],
                              d, fp)
                write_type302(info, index, channel[0], channel[1], channel[2],
                              d, fp)
                write_type303(info, index, channel[0], channel[1], channel[2],
                              u, v, w, fp)
                continue
            index += 1
            start += interval
            continue
        write_type309(info, phasecal_file, fp)
        continue
    return

if __name__ == "__main__":
    vex = Vex(sys.argv[1])

    fp = open(sys.argv[2], 'r')
    json_input = json.load(fp)
    fp.close()

    rootid = sys.argv[3]

    process_job(vex, json_input, rootid)
    pass
