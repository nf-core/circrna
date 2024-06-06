#!/usr/bin/env python3
import os
meta1 = "$meta"
depth = "$depth"
bed = "$bed"

with open("output.txt", "w") as file:
    # Write the variables to the file
    depth = os.path.abspath(depth)
    bed = os.path.abspath(bed)
    file.write(meta1)
    file.write(depth)
    file.write(bed)