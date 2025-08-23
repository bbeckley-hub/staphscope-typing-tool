#!/usr/bin/env python

# Copyright (c) 2014, Ole Lund, Technical University of Denmark
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#
# Import libraries
#
import sys, time
import os
from optparse import OptionParser
from operator import itemgetter
import re
import cPickle as pickle
#
# Functions
#
# Construct the reverse complement from a sequence
#
def reversecomplement(seq):
    '''Reverse complement'''
    comp = ''
    for s in seq:
        if s == 'A': comp = comp + 'T'
        elif s == 'T': comp = comp + 'A'
        elif s == 'C': comp = comp + 'G'
        elif s == 'G': comp = comp + 'C'
        else: comp = comp + s
    return comp[::-1]
#
# Parse command line options
#
parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="inputfilename", help="read from INFILE", metavar="INFILE") 
parser.add_option("-f", "--filterfile", dest="filterfilename", help="filter (ignore) K-mers present in FILTERFILE", metavar="FILTERFILE") 
parser.add_option("-o", "--outputfile", dest="outputfilename", help="write to OUTFILE", metavar="OUTFILE") 
parser.add_option("-k", "--kmersize", dest="kmersize", help="Size of KMER", metavar="KMERSIZE") 
parser.add_option("-t", "--homthres", dest="homthres", help="Threshold for homology reduction", metavar="HOMTHRES") 
parser.add_option("-s", "--stepsize", dest="stepsize", help="Size of step between K-mers", metavar="STEPSIZE") 
parser.add_option("-x", "--prefix", dest="prefix", help="type of prefix", metavar="PREFIX") 
#parser.add_option("-p", "--pickleoutput", dest="pickleoutput",action="store_true", help="use pickle output") 
(options, args) = parser.parse_args() 
# 
# Open file for input sequence with kmers to save in database
# 
if options.inputfilename != None:
  if options.inputfilename == "--":
    inputfile = sys.stdin
  else:
    inputfile = open(options.inputfilename,"r")
#
# Open file to filter on (kmers not to save in database)
#
if options.filterfilename != None:
  filterfile = open(options.filterfilename,"r")
#
# Harcode this to be true so I do not need to use the -p option
#
options.pickleoutput = True
if options.outputfilename != None:
  if options.pickleoutput == True:
    outputfile = open(options.outputfilename+".p", "wb")
    outputfile_lengths = open(options.outputfilename+".len.p", "wb")
    outputfile_ulengths = open(options.outputfilename+".ulen.p", "wb")
    outputfile_descriptions = open(options.outputfilename+".desc.p", "wb")
  else:
    outputfile = open(options.outputfilename,"w")
else:
  outputfile = sys.stdout
#
# Size of K-mer
#
if options.kmersize != None:
  kmersize = int(options.kmersize)
else:
  kmersize = 16
oligolen=kmersize
#
# Size of step when looking for K-mers in the sequence
#
if options.stepsize != None:
  stepsize = int(options.stepsize)
else:
  stepsize = 1
#
# Homology threshold for when to include entry in database
#
if options.homthres != None:
  homthres = float(options.homthres)
else:
  homthres = -1
#
# Prefix to use fro filtering sequences
#
if options.prefix != None:
  #prefix = int(options.prefix)
  prefix = options.prefix
  prefixlist = [prefix]
  prefixlen = len(prefixlist[0])
else:
  prefix = ''
  prefixlist = [prefix]
  prefixlen = len(prefixlist[0])
#
# Initialize statistics
#
# # of kmers
kmer_count = 0
# Start time to keep track of progress
t0 = time.time()
# Print progress
printfreq = 100000
# frequenct to save sorted list to db
dbsavefreq = 30000000
#
#
#
etta = 0.001
#
# Read sequences from filterfile ans save K-mers
#
filterseqsegments = []
filters = {}
Nfilters=0
i=0
t1 = time.time()
if options.filterfilename != None:
  sys.stdout.write("%s\n" % ("# Reading filterfile"))
  for line in filterfile:
    line = line.rstrip('\n')
    fields=line.split()
    if len(line)>1:
      if fields[0][0] == ">":
        if (i>0):
	  #
	  # Fasta entry read
	  #
          filterseq = ''.join(filterseqsegments)
	  for seq in [filterseq,reversecomplement(filterseq)]:
            start=0
            while start < len(seq)-kmersize:
              submer = seq[start:start+kmersize]
	      if prefix == seq[start:start+prefixlen]:
                if not submer in filters:
                  filters[submer] = filtername
              kmer_count += 1
              if kmer_count % printfreq == 0:
                t1 = time.time()
                sys.stdout.write("\r%s kmers (%s kmers / s)" % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t1-t0))))
                sys.stdout.flush()
              start +=stepsize	
        del filterseqsegments
        filterseqsegments = []
        i=0
        filterseq = ""
        filtername = fields[0][1:] 
      else:
        filterseqsegments.append("")
        filterseqsegments[i] = fields[0]
        i+=1
  filterseq = ''.join(filterseqsegments)
  for seq in [filterseq,reversecomplement(filterseq)]:
    start=0
    while start < len(seq)-kmersize:
      submer = seq[start:start+kmersize]
      if prefix == seq[start:start+prefixlen]:
        if not submer in filters:
          filters[submer] = filtername
      kmer_count += 1
      if kmer_count % printfreq == 0:
        t1 = time.time()
        sys.stdout.write("\r%s kmers (%s kmers / s)" % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t1-t0))))
        sys.stdout.flush()
      start +=stepsize	
  #
  # Print final statistics for filterfile
  #
  t1 = time.time()
  sys.stdout.write("\r%s kmers (%s kmers / s)" % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t1-t0))))
  sys.stdout.flush()
  sys.stdout.write("\n")
del filterseqsegments
#
# Read input file and save K-mers not in filterfile (if -f option is used)
#
kmer_count = 0
Nstored = 0
Nstored_old = Nstored
# count only one per sequence
Nustored = 0
Nustored_old = Nustored
inputseqsegments = []
inputs = {}
lengths = {}
ulengths = {}
descriptions = {}
Ninputs=0
i=0
if options.inputfilename != None:
  sys.stdout.write("%s\n" % ("# Reading inputfile"))
  for line in inputfile:
    line = line.rstrip('\n')
    fields=line.split()
    if len(line)>1:
      if fields[0][0] == ">":
        if (i>0):
          inputseq = ''.join(inputseqsegments)
	  sys.stdout.write("%s %s\n" % ("# Entry read", inputname))
	  if options.homthres != None:
            #
	    # Check for homology
	    #
	    sys.stdout.write("%s\n" % ("# Checking for homology"))
	    #
            # Make list of unique k-mers in entry
            #
            queryindex = {}
            qtotlen=0
            querymers=0
            uquerymers=0
            seqlen = len(inputseq)
            qtotlen += seqlen;
            #for qseq in [queryseq[i],reversecomplement(queryseq[i])]:
            for qseq in [inputseq]:
              for j in range(0, seqlen-oligolen+1):
                submer = qseq[j:j+oligolen]
                if prefix == qseq[j:j+prefixlen]:
                  if submer in queryindex:
                    queryindex[submer] += 1
                    querymers += 1 
                  else:
                    queryindex[submer] = 1
                    querymers += 1 
                    uquerymers += 1
            # print number of querymers
	    sys.stdout.write("#querymers:\t" + str(querymers) + "\tuquerymers:\t" + str(uquerymers) + "\n")
            # Search for matches
            #	      
            mincoverage = 1
            Nhits=0
            templateentries = {}
            templateentries_tot = {}
            for submer in queryindex:
              if submer in inputs:
                if queryindex[submer] >= mincoverage:
                  matches = inputs[submer].split(",")
		  #changed by vanessa:
		  matches = list(set(matches))
                  for match in matches:
                    Nhits += 1
                    if match in templateentries:
                      templateentries[match] += 1
                    else:
                      templateentries[match] = 1
                  umatches = list(set(matches))
                  for match in umatches:
                    if match in templateentries_tot:
                      templateentries_tot[match] += queryindex[submer]
                    else:
                      templateentries_tot[match] = queryindex[submer]
            #
            #
            #
	    frac_q = 0.0
	    hitname = ""
	    score = 0
            sortedlist= sorted(templateentries.items(), key = itemgetter(1), reverse=True)
            for template,score in sortedlist:
              #frac_q = score/(float(querymers)+etta)
	      frac_q = score/(float(uquerymers)+etta)
	      hitname= template
	      break 
            sys.stdout.write("# Max frac_q similarity of %s to %s frac_q: %s Score: %s\n" % (inputname, hitname ,frac_q, score))
            del templateentries
            del templateentries_tot
            del queryindex 
	  if (options.homthres != None and frac_q >= homthres):
            sys.stdout.write("# Skipping entry: %s in databade due to similarity to %s frac_q: %s\n" % (inputname, hitname ,frac_q))	  
	  if (options.homthres == None or (options.homthres != None and frac_q < homthres)):
	    sys.stdout.write("%s %s\n" % ("# Including entry: ", inputname))
            #
	    # Start of database update
	    #
            for seq in [inputseq,reversecomplement(inputseq)]:
              start=0
              while start < len(seq)-kmersize:
                submer = seq[start:start+kmersize]
	        if prefix == seq[start:start+prefixlen]:
	          if (options.filterfilename != None and submer not in filters) or options.filterfilename == None:
		    Nstored += 1
                    if submer in inputs:
		      if (inputs[submer].find(inputname) == -1):
		        Nustored += 1
                      inputs[submer] = inputs[submer]+","+inputname
                    else:
                      inputs[submer] = inputname
		      Nustored += 1
                kmer_count += 1
		if options.homthres == None:
                  if kmer_count % printfreq == 0:
                    t1 = time.time()
                    sys.stdout.write("\r%s kmers (%s kmers / s)" % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t1-t0))))
                    sys.stdout.flush()
                start +=stepsize	
            lengths[inputname] = Nstored - Nstored_old
	    ulengths[inputname] = Nustored - Nustored_old
	    #print inputname,lengths[inputname], Nstored, Nstored_old,"i: ",i, len(inputseq)
	    Nstored_old = Nstored
            Nustored_old = Nustored
	    #
	    # End of database update
	    #
        del inputseqsegments	
        inputseqsegments = []
        i=0
        inputseq = ""
        inputname = fields[0][1:]
	descriptions[inputname] = ' '.join(fields[1:len(fields)])
	kmer_count_old = kmer_count
      else:
        inputseqsegments.append("")
        inputseqsegments[i] = fields[0].upper()
        i = i + 1
  lengths[inputname] = Nstored - Nstored_old
  ulengths[inputname] = Nustored - Nustored_old
  inputseq = ''.join(inputseqsegments)
  #
  # Update database with last entry in inputfile
  # 
  sys.stdout.write("%s %s\n" % ("# Entry read: ",inputname))
  if options.homthres != None:
    #
    # Check for homology
    #
    sys.stdout.write("%s\n" % ("# Checking for homology"))
    #
    # Make list of unique k-mers in entry
    #
    queryindex = {}
    qtotlen=0
    querymers=0
    uquerymers=0
    seqlen = len(inputseq)
    qtotlen += seqlen;
    #for qseq in [queryseq[i],reversecomplement(queryseq[i])]:
    for qseq in [inputseq]:
      for j in range(0, seqlen-oligolen+1):
        submer = qseq[j:j+oligolen]
        if prefix == qseq[j:j+prefixlen]:
          if submer in queryindex:
            queryindex[submer] += 1
            querymers += 1 
          else:
            queryindex[submer] = 1
            querymers += 1 
            uquerymers += 1
    #
    # Search for matches
    #	      
    mincoverage = 1
    Nhits=0
    templateentries = {}
    templateentries_tot = {}
    for submer in queryindex:
      if submer in inputs:
        if queryindex[submer] >= mincoverage:
          matches = inputs[submer].split(",")
	  # changed vanessa:
	  matches = list(set(matches))
          for match in matches:
            Nhits += 1
            if match in templateentries:
              templateentries[match] += 1
            else:
              templateentries[match] = 1
          umatches = list(set(matches))
          for match in umatches:
            if match in templateentries_tot:
              templateentries_tot[match] += queryindex[submer]
            else:
              templateentries_tot[match] = queryindex[submer]
    #
    #
    #
    frac_q = 0.0
    score = 0
    template = "" 
    sortedlist= sorted(templateentries.items(), key = itemgetter(1), reverse=True)
    for template,score in sortedlist:
      #frac_q = score/(float(querymers)+etta)
      frac_q = score/(float(uquerymers)+etta)
      hitname=template
      break
    sys.stdout.write("# Max frac_q similarity of %s to %s frac_q: %s Score: %s\n" % (inputname, hitname ,frac_q, score))
    del templateentries
    del templateentries_tot
    del queryindex 
  if (options.homthres != None and frac_q >= homthres):
    sys.stdout.write("# Skipping entry: %s in databade due to similarity to %s frac_q: %s\n" % (inputname, hitname ,frac_q))	  
  if (options.homthres == None or (options.homthres != None and frac_q < homthres)):
    #
    # Start of database update
    #
    sys.stdout.write("%s %s\n" % ("# Including entry: ", inputname))
    for seq in [inputseq,reversecomplement(inputseq)]:
      start=0
      while start < len(seq)-kmersize:
        submer = seq[start:start+kmersize]
        if prefix == seq[start:start+prefixlen]:
          if (options.filterfilename != None and submer not in filters) or options.filterfilename == None:
	    Nstored += 1
            if submer in inputs:
	      if (inputs[submer].find(inputname) == -1):
	        Nustored += 1
              inputs[submer] = inputs[submer]+","+inputname
            else:
              inputs[submer] = inputname
              Nustored += 1
        kmer_count += 1
	if options.homthres == None:
          if kmer_count % printfreq == 0:
            t1 = time.time()
            sys.stdout.write("\r%s kmers (%s kmers / s)" % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t1-t0))))
            sys.stdout.flush()
        start +=stepsize	
    lengths[inputname] = Nstored - Nstored_old
    ulengths[inputname] = Nustored - Nustored_old
    #print inputname,lengths[inputname], Nstored, Nstored_old,"i: ",i, len(inputseq)
    #Nstored_old = Nstored
    #
    # End of database update with hast entry
    #
del inputseqsegments
#
# Print database of nmers
#
if options.pickleoutput == True:
  pickle.dump(inputs, outputfile,2)
  pickle.dump(lengths, outputfile_lengths,2)
  pickle.dump(ulengths, outputfile_ulengths,2)
  pickle.dump(descriptions, outputfile_descriptions,2)
#
# Print final statistics for inputfile
#
t2 = time.time()
sys.stdout.write("\r%s kmers (%s kmers / s)" % ("{:,}".format(kmer_count), "{:,}".format(kmer_count / (t2-t1))))
sys.stdout.flush()
sys.stdout.write("\n")
sys.stdout.write("# Total time used: %s s\n" % (t2-t0))
#
# Done
#
