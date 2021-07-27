#!/usr/bin/env python
import sys
cols = [int(col) for col in sys.argv[1].split(",")]
flags = [(col >= 0)*2 -1 for col in cols]
cols = [abs(c)-1 for c in cols]
threshs = [float(th) for th in sys.argv[2].split(",")]*len(cols)
#print flags,threshs
for line in sys.stdin:
    if line.startswith("#"):
        print line,
        continue

    data = line.split("\t") + [0,]
    values = [float(data[col]) for col in cols]
    # default logic: threshold met on one of the columns is sufficient,
    # if threshold is required on each column you can chain scorethresh.py

    valid = False
    for value,flag,thresh in zip(values,flags,threshs):
        if value * flag >= thresh * flag:
            valid = True

    if valid:
        print line,

