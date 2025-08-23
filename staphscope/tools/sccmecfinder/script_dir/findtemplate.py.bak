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

##########################################################################
# IMPORT LIBRARIES
##########################################################################

import sys
import time
import os
from math import sqrt, pow
from optparse import OptionParser
from operator import itemgetter
import re
import cPickle as pickle


##########################################################################
# FUNCTIONS
##########################################################################

#--------------------------------------
# reverse complement sequence:
#--------------------------------------
def reversecomplement(seq):
    '''Reverse complement'''
    comp = ''
    for s in seq:
        if s == 'A':
            comp = comp + 'T'
        elif s == 'T':
            comp = comp + 'A'
        elif s == 'C':
            comp = comp + 'G'
        elif s == 'G':
            comp = comp + 'C'
        else:
            comp = comp + s
    return comp[::-1]

#--------------------------------------
# save kmers of querysequence:
#--------------------------------------


def save_kmers(queryseq):

    global qtotlen, prefix, queryindex, querymers, uquerymers

    seqlen = len(queryseq)
    qtotlen += seqlen

    # store kmers in original and reverse complement sequence:
    for qseq in[queryseq, reversecomplement(queryseq)]:
        for j in range(0, seqlen - kmersize + 1):
            submer = qseq[j:j + kmersize]
            if prefix == qseq[j:j + prefixlen]:
                if submer in queryindex:
                    if submer in templates:
                        queryindex[submer] += 1
                    querymers += 1
                else:
                    if submer in templates:
                        queryindex[submer] = 1
                    querymers += 1
                    uquerymers += 1

#-------------------------------------
# search for matches:
#-------------------------------------


def find_matches():

    global queryindex, mincoverage, templates

    templateentries = {}
    templateentries_tot = {}
    Nhits = 0

    for submer in queryindex:
        if queryindex[submer] >= mincoverage:
            matches = templates[submer].split(",")

            # get list of unique matches:
            umatches = list(set(matches))

            # Nhits = sum of scores over all templates:
            Nhits = Nhits + len(umatches)

            for match in umatches:

                # get unique scores:
                if match in templateentries:
                    templateentries[match] += 1
                else:
                    templateentries[match] = 1

                # get total amount of kmers found in template (total score):
                if match in templateentries_tot:
                    templateentries_tot[match] += queryindex[submer]
                else:
                    templateentries_tot[match] = queryindex[submer]

    return(templateentries, templateentries_tot, Nhits)

#------------------------------------------------
# Conservative two sided p-value from z-score:
#------------------------------------------------


def z_from_two_samples(r1, n1, r2, n2):
    '''Comparison of two fractions, Statistical methods in medical research, Armitage et al. p. 125'''
    #
    # r1: positives in sample 1
    # n1: size of sample 1
    # r2: positives in sample 2
    # n2: size of sample 2

    p1 = float(r1) / (float(n1) + etta)
    p2 = float(r2) / (float(n2) + etta)
    q1 = 1 - p1
    q2 = 1 - p2
    p = (float(r1) + r2) / (n1 + n2 + etta)
    q = 1 - p
    z = (p1 - p2) / sqrt(p * q * (1 / (n1 + etta) + 1 / (n2 + etta)) + etta)

    return z

#-------------------------------------------------
# Conservative two sided p-value from z-score:
#-------------------------------------------------


def fastp(z):
    '''Conservative two sided p-value from z-score'''
    if z > 10.7016:
        p = 1e-26
    elif z > 10.4862:
        p = 1e-25
    elif z > 10.2663:
        p = 1e-24
    elif z > 10.0416:
        p = 1e-23
    elif z > 9.81197:
        p = 1e-22
    elif z > 9.5769:
        p = 1e-21
    elif z > 9.33604:
        p = 1e-20
    elif z > 9.08895:
        p = 1e-19
    elif z > 8.83511:
        p = 1e-18
    elif z > 8.57394:
        p = 1e-17
    elif z > 8.30479:
        p = 1e-16
    elif z > 8.02686:
        p = 1e-15
    elif z > 7.73926:
        p = 1e-14
    elif z > 7.4409:
        p = 1e-13
    elif z > 7.13051:
        p = 1e-12
    elif z > 6.8065:
        p = 1e-11
    elif z > 6.46695:
        p = 1e-10
    elif z > 6.10941:
        p = 1e-9
    elif z > 5.73073:
        p = 1e-8
    elif z > 5.32672:
        p = 1e-7
    elif z > 4.89164:
        p = 1e-6
    elif z > 4.41717:
        p = 1e-5
    elif z > 3.89059:
        p = 1e-4
    elif z > 3.29053:
        p = 1e-3
    elif z > 2.57583:
        p = 0.01
    elif z > 1.95996:
        p = 0.05
    elif z > 1.64485:
        p = 0.1
    else:
        p = 1.0
    return p


##########################################################################
#	DEFINE GLOBAL VARIABLES
##########################################################################

global templates, queryindex, prefix, qtotlen, mincoverage, querymers, uquerymers


##########################################################################
# PARSE COMMAND LINE OPTIONS
##########################################################################

parser = OptionParser()
parser.add_option("-i", "--inputfile", dest="inputfilename",help="read from INFILE", metavar="INFILE")
parser.add_option("-t", "--templatefile", dest="templatefilename",help="read from TEMFILE", metavar="TEMFILE")
parser.add_option("-o", "--outputfile", dest="outputfilename",help="write to OUTFILE", metavar="OUTFILE")
parser.add_option("-k", "--kmersize", dest="kmersize",help="Size of k-mer, default 16", metavar="KMERSIZE")
parser.add_option("-x", "--prefix", dest="prefix",help="prefix, e.g. ATGAC, default none", metavar="_id")
parser.add_option("-a", "--printall", dest="printall", action="store_true",help="Print matches to all templates in templatefile unsorted")
parser.add_option("-w", "--winnertakesitall", dest="wta", action="store_true",help="kmer hits are only assigned to most similar template")
parser.add_option("-e", "--evalue", dest="evalue", help="Maximum E-value", metavar="EVALUE")
(options, args) = parser.parse_args()

# set up prefix filtering:
if options.prefix != None:
    prefix = options.prefix
else:
    prefix = ''
prefixlen = len(prefix)


# get e-value:
if options.evalue != None:
    evalue = float(options.evalue)
else:
    evalue = float(0.05)


# Open queryfile:
t0 = time.time()
if options.inputfilename != None:
    if options.inputfilename == "--":
        inputfile = sys.stdin
    else:
        inputfile = open(options.inputfilename, "r")

# open templatefile:
if options.templatefilename != None:
    templatefile = open(options.templatefilename + ".p", "rb")
    templatefile_lengths = open(options.templatefilename + ".len.p", "rb")
    try:
        templatefile_ulengths = open(
            options.templatefilename + ".ulen.p", "rb")
    except:
        # do nothing
        two = 2
    templatefile_descriptions = open(
        options.templatefilename + ".desc.p", "rb")
else:
    sys.exit("No template file specified")

# open outputfile:
if options.outputfilename != None:
    outputfile = open(options.outputfilename, "w")
else:  # If no output filename choose the same as the input filename
    outputfilename = os.path.splitext(options.inputfilename)[0]
    outputfile = open(outputfilename, "w")

# get kmer size:
if options.kmersize != None:
    kmersize = int(options.kmersize)
else:
    kmersize = 16


##########################################################################
# READ DATABASE OF TEMPLATES
##########################################################################

templates = {}
Ntemplates = 0


# Read Template file:
sys.stdout.write("%s\n" % ("# Reading database of templates"))
templates = pickle.load(templatefile)
templates_lengths = pickle.load(templatefile_lengths)
try:
    templates_ulengths = pickle.load(templatefile_ulengths)
except:
    sys.stderr.write('No ulen.p file found for database')
    SystemExit()
templates_descriptions = pickle.load(templatefile_descriptions)

# Count number of k-mers, and sum of unique k-mers over all templates:
template_tot_len = 0
template_tot_ulen = 0
Ntemplates = 0
# length added
for name in templates_lengths:
    template_tot_len += templates_lengths[name]
    template_tot_ulen += templates_ulengths[name]
    Ntemplates += 1


##########################################################################
# READ INPUTFILE
##########################################################################

queryseq = ""
queryseqsegments = []
Nquerys = 0
queryindex = {}
qtotlen = 0
querymers = 0
uquerymers = 0
i = 0

if options.inputfilename != None:
    sys.stdout.write("%s\n" % ("# Reading inputfile"))
    for line in inputfile:
        fields = line.split()
        if len(fields) >= 1:
            # FASTA file:
            if fields[0][0] == ">":
                if (i > 0):
                    queryseq = ''.join(queryseqsegments)

                    # Update dictionary of kmers
                    save_kmers(queryseq)

                del queryseqsegments
                queryseqsegments = []
                i = 0

            # Fastq file:
            elif fields[0][0] == "@":
                # Fastq file
                if (i > 0):
                    queryseq = ''.join(queryseqsegments)

                    # Update dictionary of kmers:
                    save_kmers(queryseq)

                del queryseqsegments
                queryseqsegments = []
                i = 0
                queryseq = ""

                try:
                    line = inputfile.next()
                    fields = line.split()
                    if len(fields) == 0:
                        continue
                    queryseqsegments.append("")
                    queryseqsegments[i] = fields[0]
                    i += 1
                    line = inputfile.next()
                    line = inputfile.next()
                except:
                    break
            else:
                queryseqsegments.append("")
                queryseqsegments[i] = fields[0].upper()
                i += 1

    queryseq = ''.join(queryseqsegments)

    # Update dictionary of K-mers:
    save_kmers(queryseq)

del queryseqsegments


##########################################################################
# SEARCH FOR MATCHES
##########################################################################

sys.stdout.write("%s\n" % ("# Searching for matches of input in template"))
mincoverage = 1
templateentries = {}
templateentries_tot = {}
Nhits = 0

(templateentries, templateentries_tot, Nhits) = find_matches()


##########################################################################
#	DO STATISTICS
##########################################################################


minscore = 0
etta = 1.0e-8  # etta is a small number to avoid division by zero

# report search statistics:
sys.stdout.write("%s\n" % ("# Search statistics"))
sys.stdout.write("%s\n" % ("# Total number of hits: %s") % (Nhits))
sys.stdout.write("%s\n" % ("# Total number of kmers in templates : %s") % (template_tot_len))
sys.stdout.write("%s\n" % ("# Minimum number of k-mer hits to report template: %s") % (minscore))
sys.stdout.write("%s\n" % ("# Maximum multiple testing corrected E-value to report match : %s") % (evalue))
sys.stdout.write("%s\n" % ("# Printing best matches"))


# print heading of outputfile:
if options.wta != True:
    outputfile.write("#Template\tScore\tExpected\tz\tp_value\tquery coverage [%]\ttemplate coverage [%]\tdepth\tKmers in Template\tDescription\n")
elif options.wta == True:
    outputfile.write("#Template\tScore\tExpected\tz\tp_value\tquery coverage [%]\ttemplate coverage [%]\tdepth\ttotal query coverage [%]\ttotal template coverage [%]\ttotal depth\tKmers in Template\tDescription\n")


##########################################################################
#	STANDARD SCORING SCHEME
##########################################################################


if not options.wta == True:

    sortedlist = sorted(
        templateentries.items(), key=itemgetter(1), reverse=True)
    for template, score in sortedlist:
        if score > minscore:
            expected = float(
                Nhits) * float(templates_ulengths[template]) / float(template_tot_ulen)
            #z = (score - expected)/sqrt(score + expected+etta)
            #p  = fastp(z)
            #
            # If expected < 1 the above poisson approximation is a poor model
            # Use instead: probabilyty of seing X hits is p**X if probability
            # of seing one hit is p (like tossing a dice X times)
            #
            # if expected <1:
            #  p = expected**score
            #
            # Comparison of two fractions, Statistical methods in medical
            # research, Armitage et al. p. 125:
            z = z_from_two_samples(
                score, templates_ulengths[template], Nhits, template_tot_ulen)
            p = fastp(z)

            # Correction for multiple testing:
            p_corr = p * Ntemplates
            frac_q = ( score / (float(uquerymers) + etta)) * 100
            frac_d = (score / (templates_ulengths[template] + etta)) * 100
            coverage = templateentries_tot[
                template] / float(templates_lengths[template])

            # print str(p) + " " + str(p_corr) + " " + str(evalue)
            # p_corr=0
            if p_corr <= evalue:
                outputfile.write("%-12s\t%8d\t%8d\t%8.2f\t%4.1e\t%8.2f\t%8.2f\t%8.2f\t%8d\t%s\n" %
                                 (template, score, int(round(expected)), round(z, 1), p_corr, frac_q, frac_d, coverage, templates_ulengths[template], templates_descriptions[template].strip()))


##########################################################################
#	WINNER TAKES IT ALL
##########################################################################

if options.wta == True:

    w_templateentries = templateentries.copy()
    w_templateentries_tot = templateentries_tot.copy()
    w_Nhits = Nhits

    maxhits = 100
    hitcounter = 1
    stop = False
    while not stop:
        hitcounter += 1
        if hitcounter > maxhits:
            stop = True
        sortedlist = sorted(
            w_templateentries.items(), key=itemgetter(1), reverse=True)
        for template, score in sortedlist:
            if score > minscore:
                expected = float(
                    w_Nhits) * float(templates_ulengths[template]) / float(template_tot_ulen)
                #z = (score - expected)/sqrt(score + expected+etta)
                #p  = fastp(z)
                #
                # If expected < 1 the above poisson approximation is a poor model
                #
                # if expected <1:
                #  p = expected**score
                #
                # Comparison of two fractions, Statistical methods in medical
                # research, Armitage et al. p. 125:
                z = z_from_two_samples(
                    score, templates_ulengths[template], w_Nhits, template_tot_ulen)
                p = fastp(z)

                # correction for multiple testing:
                p_corr = p * Ntemplates
                # print score,float(uquerymers),etta
                frac_q = (score / (float(uquerymers) + etta)) * 100
                frac_d = (score / (templates_ulengths[template] + etta)) * 100
                coverage = w_templateentries_tot[
                    template] / float(templates_lengths[template])

                # calculate total values:
                tot_frac_q = (templateentries[template] / (float(uquerymers) + etta)) * 100
                tot_frac_d = (templateentries[template] / (templates_ulengths[template] + etta)) * 100
                tot_coverage = templateentries_tot[template] / float(templates_lengths[template])

                # print results to outputfile:
                if p_corr <= evalue:
                    outputfile.write("%-12s\t%8d\t%8d\t%8.1f\t%4.2e\t%8.2f\t%8.2f\t%4.2f\t%8.2f\t%8.2f\t%4.2f\t%8d\t%s\n" %
                                     (template, score, int(round(expected)), round(z, 1), p_corr, frac_q, frac_d, coverage, tot_frac_q, tot_frac_d, tot_coverage, templates_ulengths[template], templates_descriptions[template].strip()))

                    # remove all kmers in best hit from queryindex
                    for submer in queryindex:
                        matches = templates[submer].split(",")
                        if template in matches:
                            queryindex[submer] = 0

                    # find best hit like before:
                    del w_templateentries
                    del w_templateentries_tot
                    w_templateentries = {}
                    w_templateentries_tot = {}
                    w_Nhits = 0

                    (w_templateentries, w_templateentries_tot,
                     w_Nhits) = find_matches()

                else:
                    stop = True
            break

##########################################################################
# CLOSE FILES
##########################################################################

t1 = time.time()

sys.stdout.write("\r# %s kmers (%s kmers / s). Total time used: %s sec" %("{:,}".format(querymers), "{:,}".format(querymers / (t1 - t0)), int(t1 - t0)))
sys.stdout.flush()
sys.stdout.write("\n")
sys.stdout.write("# Closing files\n")
