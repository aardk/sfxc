#!/usr/bin/env python
#
# Quick and dirty sigproc filterbank data importer

import numpy as np

import sys, struct, string
import os.path
from optparse import OptionParser


def get_field(infile):
    startstop = ["HEADER_START", "HEADER_END"]
    strfields = ["rawdatafile", "source_name"]
    intfields = ["machine_id", "telescope_id", "data_type", "nchans", "nbeams", "ibeam", "nbits", "nifs"]
    doublefields = ["src_raj", "src_dej", "az_start", "za_start", "fch1", "foff", "tstart", "tsamp"]

    hlen_str = infile.read(4)
    hlen = struct.unpack("i", hlen_str)[0]
    field = infile.read(hlen).decode()
    arg = None
    if field in strfields:
        l = struct.unpack("i", infile.read(4))[0]
        arg = infile.read(l).decode()
    elif field in intfields:
        arg = struct.unpack("i", infile.read(4))[0]
    elif field in doublefields:
        arg =struct.unpack("d", infile.read(8))[0]
    elif field not in startstop:
        if field in string.printable:
            print("Error, unknown field in header: ", field)
        else:
            print("Error, bad header")
        sys.exit(1)
    return (field, arg)

def read_header(infile):
    header_start = "HEADER_START"
    header_end = "HEADER_END"
    (field, arg) = get_field(infile)
    if(field != header_start):
        print("Invalid sigproc header")
        sys.exit(1)
    header = {}
    while True:
        (field, arg) = get_field(infile)
        if field != header_end:
            header[field] = arg
        else:
            break
    pos = infile.tell()
    infile.seek(0)
    binheader = infile.read(pos)
    return (header, binheader)

def loaddata(filename):
    infile = open(filename, 'rb')
    (header, binheader) = read_header(infile)
    ncol = header["nchans"] 
    pos = infile.tell()
    infile.seek(0,2)
    pos_end = infile.tell()
    infile.seek(pos)
    # NB: I assume single precision floats!!!
    nrows = (pos_end-pos+1)//(ncol*4)
    data = np.fromfile(infile, dtype='f4', count=ncol*nrows)
    data = data.reshape([nrows, ncol])
    infile.seek(pos)
    return (header, binheader, data)

def writedata(filename, binheader, data):
    if os.path.exists(filename):
        while True:
            answer = input("File already exists, are you sure you want to overwrite? YES/NO : ").upper()
            if (answer != "YES") and (answer != "NO"):
                print("Please type Yes or No")
            else:
                break
        if answer == "NO":
            return
    outfile = open(filename, 'w')
    outfile.write(binheader)
    if data.dtype != np.float32:
        data.astype(np.float32).tofile(outfile)
    else:
        data.tofile(outfile)

if __name__ == "__main__":
    usage = "%prog <inputfile>"
    parser = OptionParser(usage=usage)
    #parser.add_option("-f", "--file", dest="filename")
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.error("Insufficient number of arguments")
        sys.exit(1)

    (header, binheader, data) = loaddata(args[0])
