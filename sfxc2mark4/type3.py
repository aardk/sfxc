import json
import math
import struct
import sys
import time
import urlparse

import numpy as np
import scipy as sp
from scipy import fftpack, signal

from vex import Vex

class ScanInfo:
    def __init__(self, vex, station, scan):
        self.frequencies = []
        self.channels = []
        self.max_bandwidth = 0
        self.upper = False
        self.lower = False

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
type309_header = "!IIdd"
type309_channel = "!8sd128i"

# Convert UNIX time to a Mk4 date
def time2date(secs):
    tm = time.gmtime(secs)
    return struct.pack("!HHHHf", tm.tm_year, tm.tm_yday, tm.tm_hour,
                       tm.tm_min, tm.tm_sec)

def write_type300(info, start, interval, nsplines, out_fp):
    # Write Type300 header
    buf = struct.pack(type3, '300', '00')
    out_fp.write(buf)
    buf = struct.pack(type300, 'N', 'Ny', 'NY', time2date(start), interval,
                      nsplines)
    out_fp.write(buf)
    return

def write_type301(info, frequency, sideband, pol, out_fp):
    # Write Type301 header
    buf = struct.pack(type3, '301', '00')
    out_fp.write(buf)
    chan_name = info.chan_name(frequency, sideband, pol)
    coefficients = [0.0] * 6
    buf = struct.pack(type301, 0, chan_name, *coefficients)
    out_fp.write(buf)
    return

def write_type302(info, frequency, sideband, pol, out_fp):
    # Write Type301 header
    buf = struct.pack(type3, '301', '00')
    out_fp.write(buf)
    chan_name = info.chan_name(frequency, sideband, pol)
    coefficients = [0.0] * 6
    buf = struct.pack(type302, 0, chan_name, *coefficients)
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
            norm = math.sqrt(2) * 128.0
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

station = 'Ny'
scan = '191-0534'

fp = open(sys.argv[2], 'r')
json_input = json.load(fp)
fp.close()

vex = Vex(sys.argv[1])
info = ScanInfo(vex, station, scan)

phasecal_uri = json_input['phasecal_file']
phasecal_file = urlparse.urlparse(phasecal_uri).path

filename = "1234/191-0534/N..mtuwvo"
fp = open(filename, 'w')
buf = struct.pack(ident, '000', '01', "2001001-123456", filename)
fp.write(buf)

write_type300(info, 1404970440, 240, 1, fp)
for channel in info.channels:
    write_type301(info, channel[0], channel[1], channel[2], fp)
    write_type302(info, channel[0], channel[1], channel[2], fp)
    continue
write_type309(info, phasecal_file, fp)
