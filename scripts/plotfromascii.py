#! /usr/bin/env python

import sys
import os
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

commandname, inputfile, outputfile = sys.argv

'''
Takes path to an ASCII file containing PSR data in the usual (time, real, imag)
format. Returns a plot, saved to the path specified by the second argument.
'''

# load data
input_data = np.loadtxt(inputfile)
time = input_data[:, 0]
data = input_data[:, 1] + 1j * input_data[:, 2]

# plot
plt.figure()
plt.plot(time, data.real, 'b', label='Re')
plt.plot(time, data.imag, 'r', label='Im')
plt.ylabel('h (strain)')
plt.xlabel('Time (GPS)')
plt.xlim(time[0], time[-1])
plt.legend(numpoints=1)

plt.savefig(outputfile, bbox_inches='tight')
plt.close()
