#!/usr/bin/env python
import sys,os
from collections import defaultdict


def read_to_hash(fname,ds=0,de=0,flank=0,cover=False):
    #print "loading",fname
    pos = {}
    for line in file(fname):
        if line.startswith("#"):
            continue
        line = line.strip()
        chrom,start,end,name,score,sense = line.split('\t')[:6]
        start,end = int(start)+ds,int(end)+de

        #print (chrom,start,end,sense)
        pos[(chrom,start,end,sense)] = line
        
        if flank:
            for x in xrange(flank):
                pos[(chrom,start-x,end,sense)] = line
                pos[(chrom,start+x,end,sense)] = line
                pos[(chrom,start,end-x,sense)] = line
                pos[(chrom,start,end+x,sense)] = line
        
        #if cover:
            #for x in xrange
    return pos

N = defaultdict(int)
anna = read_to_hash(sys.argv[1],flank=0)
N['unique_input1'] = len(anna)
#print len(anna.keys())

marv = read_to_hash(sys.argv[2])
N['unique_input2'] = len(marv)
#print len(marv.keys())

for circ,line in marv.items():
    if circ in anna:
        if len(sys.argv) > 3:
            print "%s\t%s" % (anna[circ].split('\t')[3],line.split('\t')[3])
        else:
            print anna[circ]
        #print "M",line
        N['overlap'] += 1        
        del anna[circ]
    else:
        N['input2_not_in_input1'] += 1
    #print len(anna.keys())
        
for k,l in anna.items():
    #if "HEK" in l:
        print "MISSING\t%s" % l
        N['input1_not_in_input2'] += 1

for k in sorted(N.keys()):
    sys.stderr.write("%s\t%d\n" % (k,N[k]))
        
found = N['overlap']
detected = N['unique_input2']
total = N['unique_input1']
fp = N['input2_not_in_input1']

#print "#sensitivity %d/%d = %.2f %%" % (found,total,float(found)/total*100)
#print "#FDR %d/%d = %.2f %%" % (fp,detected,float(fp)/detected*100)