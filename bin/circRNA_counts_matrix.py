#!/usr/bin/env python

import sys, glob
from collections import defaultdict

"""
Script to parse circRNA quantification bed files, creating a count matrix suitable for statistical analyses
requires python 2.7
Author: Declan Bennett @declan93
License: MIT
"""


par = sys.argv
hld = defaultdict(list)
samps = defaultdict(list)

path = "./"
files = [f for f in glob.glob("*.bed")]

for fi in files:
    name = fi.strip("./").split(".bed")[0]
    with open(fi) as IN:
        for li in IN:
            fds = li.strip().split()
            k = "_".join(fds[0:4])
            if k in samps:
                samps[k].append((name, fds[4]))
            else:
                samps[k] = [(name, fds[4])]


def Diff(list1, list2):
    return list(list(set(list1) - set(list2)) + list(set(list2) - set(list1)))


tmp = [x.replace(".bed", "") for x in files]
print("Chr\tStart\tStop\tStrand\t" + "\t".join([x for x in tmp]))

for k, v in samps.items():
    disj = Diff(tmp, [a_tuple[0] for a_tuple in v])
    for val in disj:
        v.append((val, "0"))
    li = []
    for h in tmp:
        for x in v:
            if x[0] == h:
                li.append(x[1])
    print("\t".join(k.split("_")) + "\t" + "\t".join(li))
sys.exit()
