#!/usr/bin/env python3

# remember to mount the working directory so python within the container can access it
import sys
from scipy.signal import argrelmin
import numpy as np

hist_dict = {}

with open(sys.argv[1]) as i:
    for line in i:
        line = line.strip()
        # Keys are sequencing coverage values, values are counts of that coverage value
        hist_dict[line.split(',')[0]] = int(line.split(',')[1])

"""
In an unphased genome assembly, the haplotigs will exhibit around
half of the coverage of the diploid co-assembled contigs. This can
be exploited to separate haplotigs: see the tutorial in

https://bitbucket.org/mroachawri/purge_haplotigs/wiki/Tutorial

Under the assumption that the input assembly exhibits the expected
pattern of coverage, the low, midpoint, and high cutoffs required by
phase_haplotigs can be estimated as the local minima of the coverage
graph when the order of the coverage vs. count equals 3.

This assumption cannot be true in all situations. Check the graphical
output of purge_haplotigs readhist.
"""
local_minima = argrelmin(np.array(list(hist_dict.values())), order = 3)
local_minima_x = [list(hist_dict.keys())[X] for X in local_minima[0]]
local_minima_y = [list(hist_dict.values())[X] for X in local_minima[0]]

#with open("critical_values.csv", "w") as cv_out:
#    cv_out.write("Critical point, X, Y\nLow, %s, %s\nMidpoint, %s, %s\nHigh, %s, %s\n" % (local_minima_x[0], local_minima_y[0], local_minima_x[1], local_minima_y[1], local_minima_x[2], local_minima_y[2]))

#with open("low_mid_high.csv", "w") as lmh_out:
#    lmh_out.write("%s,%s,%s\n" % (local_minima_x[0], local_minima_x[1], local_minima_x[2]))
print("%s\t%s\t%s" % (local_minima_x[0], local_minima_x[1], local_minima_x[2]))
