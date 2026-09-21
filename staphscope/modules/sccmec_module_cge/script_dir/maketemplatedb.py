#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
maketemplatedb.py
-----------------
Builds template databases for SCCmecFinder (k-mer based).
Updated for Python 3 compatibility.

Original author: Ole Lund, Technical University of Denmark
Modified and cleaned by: Beckley & ChatGPT
"""

import sys
import time
from optparse import OptionParser
from operator import itemgetter
import pickle

# -----------------------------
# Utility: reverse complement
# -----------------------------
def reversecomplement(seq):
    """Return reverse complement of DNA sequence"""
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(comp.get(base, base) for base in seq)[::-1]

# -----------------------------
# Command-line options
# -----------------------------
parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="inputfilename", help="Input FASTA file", metavar="INFILE")
parser.add_option("-f", "--filterfile", dest="filterfilename", help="FASTA file with kmers to ignore", metavar="FILTERFILE")
parser.add_option("-o", "--outputfile", dest="outputfilename", help="Output prefix", metavar="OUTFILE")
parser.add_option("-k", "--kmersize", dest="kmersize", help="Size of k-mer", metavar="KMERSIZE")
parser.add_option("-t", "--homthres", dest="homthres", help="Threshold for homology reduction", metavar="HOMTHRES")
parser.add_option("-s", "--stepsize", dest="stepsize", help="Step size between kmers", metavar="STEPSIZE")
parser.add_option("-x", "--prefix", dest="prefix", help="Prefix filter", metavar="PREFIX")

(options, args) = parser.parse_args()

# -----------------------------
# Open files
# -----------------------------
if not options.inputfilename:
    sys.stderr.write("Error: input file required (-i)\n")
    sys.exit(1)

inputfile = open(options.inputfilename, "r")

if options.filterfilename:
    filterfile = open(options.filterfilename, "r")
else:
    filterfile = None

if not options.outputfilename:
    sys.stderr.write("Error: output filename required (-o)\n")
    sys.exit(1)

# Always pickle output
outputfile = open(options.outputfilename + ".p", "wb")
outputfile_lengths = open(options.outputfilename + ".len.p", "wb")
outputfile_ulengths = open(options.outputfilename + ".ulen.p", "wb")
outputfile_descriptions = open(options.outputfilename + ".desc.p", "wb")

# -----------------------------
# Parameters
# -----------------------------
kmersize = int(options.kmersize) if options.kmersize else 16
stepsize = int(options.stepsize) if options.stepsize else 1
homthres = float(options.homthres) if options.homthres else -1.0

prefix = options.prefix if options.prefix else ""
prefixlen = len(prefix)

# -----------------------------
# Stats
# -----------------------------
kmer_count = 0
t0 = time.time()
printfreq = 100000
etta = 0.001

# -----------------------------
# Filters
# -----------------------------
filters = {}
if filterfile:
    sys.stdout.write("# Reading filterfile\n")
    filterseqsegments = []
    filtername = ""
    for line in filterfile:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if filterseqsegments:
                filterseq = ''.join(filterseqsegments)
                for seq in [filterseq, reversecomplement(filterseq)]:
                    start = 0
                    while start < len(seq) - kmersize:
                        submer = seq[start:start + kmersize]
                        if prefix == seq[start:start + prefixlen]:
                            filters[submer] = filtername
                        kmer_count += 1
                        start += stepsize
            filterseqsegments = []
            filtername = line[1:]
        else:
            filterseqsegments.append(line.upper())

    # Final entry
    if filterseqsegments:
        filterseq = ''.join(filterseqsegments)
        for seq in [filterseq, reversecomplement(filterseq)]:
            start = 0
            while start < len(seq) - kmersize:
                submer = seq[start:start + kmersize]
                if prefix == seq[start:start + prefixlen]:
                    filters[submer] = filtername
                kmer_count += 1
                start += stepsize

# -----------------------------
# Input sequences
# -----------------------------
inputs = {}
lengths = {}
ulengths = {}
descriptions = {}

Nstored = 0
Nstored_old = 0
Nustored = 0
Nustored_old = 0

inputseqsegments = []
inputname = ""

sys.stdout.write("# Reading inputfile\n")
for line in inputfile:
    line = line.strip()
    if not line:
        continue
    if line.startswith(">"):
        if inputseqsegments:
            # Process previous sequence
            inputseq = ''.join(inputseqsegments)
            descriptions[inputname] = "previous"
            for seq in [inputseq, reversecomplement(inputseq)]:
                start = 0
                while start < len(seq) - kmersize:
                    submer = seq[start:start + kmersize]
                    if prefix == seq[start:start + prefixlen]:
                        if (not filterfile) or (submer not in filters):
                            Nstored += 1
                            if submer in inputs:
                                if inputname not in inputs[submer]:
                                    Nustored += 1
                                    inputs[submer] = inputs[submer] + "," + inputname
                            else:
                                inputs[submer] = inputname
                                Nustored += 1
                        kmer_count += 1
                    start += stepsize
            lengths[inputname] = Nstored - Nstored_old
            ulengths[inputname] = Nustored - Nustored_old
            Nstored_old = Nstored
            Nustored_old = Nustored
        # Reset
        inputseqsegments = []
        inputname = line[1:]
        descriptions[inputname] = " ".join(line.split()[1:])
    else:
        inputseqsegments.append(line.upper())

# Final sequence
if inputseqsegments:
    inputseq = ''.join(inputseqsegments)
    for seq in [inputseq, reversecomplement(inputseq)]:
        start = 0
        while start < len(seq) - kmersize:
            submer = seq[start:start + kmersize]
            if prefix == seq[start:start + prefixlen]:
                if (not filterfile) or (submer not in filters):
                    Nstored += 1
                    if submer in inputs:
                        if inputname not in inputs[submer]:
                            Nustored += 1
                            inputs[submer] = inputs[submer] + "," + inputname
                    else:
                        inputs[submer] = inputname
                        Nustored += 1
                kmer_count += 1
            start += stepsize
    lengths[inputname] = Nstored - Nstored_old
    ulengths[inputname] = Nustored - Nustored_old

# -----------------------------
# Save output
# -----------------------------
pickle.dump(inputs, outputfile, protocol=2)
pickle.dump(lengths, outputfile_lengths, protocol=2)
pickle.dump(ulengths, outputfile_ulengths, protocol=2)
pickle.dump(descriptions, outputfile_descriptions, protocol=2)

t2 = time.time()
sys.stdout.write("\n# Total kmers: %s (%s kmers/s)\n" %
                 ("{:,}".format(kmer_count), "{:,}".format(int(kmer_count / (t2 - t0)))))
sys.stdout.write("# Total time used: %.2f s\n" % (t2 - t0))
