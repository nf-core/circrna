#!/usr/bin/env python
import sys
cols = [int(col) for col in sys.argv[1].split(",")]
flags = [(col >= 0)*2 -1 for col in cols]
cols = [abs(c)-1 for c in cols]


name = "sum"
if len(sys.argv) > 2:
    name = sys.argv[2]
    
for line in sys.stdin:
    if line.startswith("#"):
        print line.rstrip()+'\t%s' % name
        continue
    data = line.rstrip().split("\t")
    values = [float(data[col]) for col in cols]

    res = 0.
    for value,flag in zip(values,flags):
        res += value * flag

    data.append(str(res))
    print "\t".join(data)

