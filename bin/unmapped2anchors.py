#!/usr/bin/env python
import pysam
import numpy
import os,sys,re
from optparse import OptionParser
from logging import error

COMPLEMENT = {
    'a' : 't',
    't' : 'a',
    'c' : 'g',
    'g' : 'c',
    'k' : 'm',
    'm' : 'k',
    'r' : 'y',
    'y' : 'r',
    's' : 's',
    'w' : 'w',
    'b' : 'v',
    'v' : 'b',
    'h' : 'd',
    'd' : 'h',
    'n' : 'n',
    'A' : 'T',
    'T' : 'A',
    'C' : 'G',
    'G' : 'C',
    'K' : 'M',
    'M' : 'K',
    'R' : 'Y',
    'Y' : 'R',
    'S' : 'S',
    'W' : 'W',
    'B' : 'V',
    'V' : 'B',
    'H' : 'D',
    'D' : 'H',
    'N' : 'N',
}

def complement(s):
    return "".join([COMPLEMENT[x] for x in s])

def rev_comp(seq):
    return complement(seq)[::-1]

usage = """

  %prog <alignments.bam> > unmapped_anchors.qfa

Extract anchor sequences from unmapped reads. Optionally permute.
"""

parser = OptionParser(usage=usage)
parser.add_option("-a","--anchor",dest="asize",type=int,default=20,help="anchor size")
parser.add_option("-q","--minqual",dest="minqual",type=int,default=5,help="min avg. qual along both anchors (default=5)")
parser.add_option("-r","--rev",dest="rev",type="choice",choices=["A","B","R","N","C","P"],default="N",help="P-ermute read parts or reverse A,B,R,C,N for control purposes.")

options,args = parser.parse_args()

import random
perm_A = []
perm_I = []
perm_B = []

perm_burn_in = []
N_perm = 100

def randomchoice(l):
    return l.pop(random.randint(0,len(l)-1))

def passthru(x):
    return x

def reverse(x):
    return x[::-1]

funcs = {
    'A' : (passthru,reverse,passthru),
    'B' : (passthru,passthru,reverse),
    'R' : (reverse,passthru,passthru),
    'N' : (passthru,passthru,passthru),
    'P' : (passthru,passthru,passthru),
    'C' : (passthru,passthru,passthru),
}

read_f,A_f,B_f = funcs[options.rev]

def handle_read(read):
    if not read.is_unmapped:
        return

    seq,qual = read_f(read.seq),read_f(read.qual)
    
    # minimal quality scores
    nq = numpy.fromstring(qual,dtype=numpy.uint8) - 35
    if nq[:options.asize].mean() < options.minqual or nq[-options.asize:].mean() < options.minqual:
        # read is junk
        #print "qual.fail",nq[:options.asize].mean(),nq[-options.asize:].mean()
        return
    
    if options.rev == "P":
        perm_A.append((seq[:options.asize],qual[:options.asize]))
        perm_B.append((seq[-options.asize:],qual[-options.asize:]))
        perm_I.append((seq[options.asize:-options.asize:],qual[options.asize:-options.asize:]))
        
        if len(perm_burn_in) < N_perm:
            # collect some reads for permutation control first.
            perm_burn_in.append(read)
            return

        A_seq,A_qual = randomchoice(perm_A)
        B_seq,B_qual = randomchoice(perm_B)
        I_seq,I_qual = randomchoice(perm_I)
        
        seq,qual = A_seq+I_seq+B_seq,A_qual+I_qual+B_qual
        
    if options.rev == "C":
        seq = rev_comp(seq)
        qual = reverse(qual)
    
    print "@%s_A__%s" % (read.qname,seq)
    print A_f(seq[:options.asize])
    print "+"
    print A_f(qual[:options.asize])
        
    print "@%s_B" % read.qname
    print B_f(seq[-options.asize:])
    print "+"
    print B_f(qual[-options.asize:])

for read in pysam.Samfile(args[0],'rb'):
    handle_read(read)

for read in perm_burn_in:
    handle_read(read)
