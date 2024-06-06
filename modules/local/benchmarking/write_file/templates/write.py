#!/usr/bin/env python3
import os
m = "$meta"
v = "$value"


with open("output.txt", "w") as file:
    # Write the variables to the file
    v = os.path.abspath(v)
    file.write(m + '\t')
    file.write(v)