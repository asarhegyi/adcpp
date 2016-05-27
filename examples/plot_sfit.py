#!/bin/python

import numpy as np
import scipy.optimize as optimize
import scipy.fftpack as fftpack
import matplotlib.pyplot as plt
import sys

time = [] #time 
data = [] #measured data
fitted = [] #fitted data
residuals = []
#filename = sys.argv[1]

#filename = './input_data.dat'
filename = sys.argv[1]
print("Loading "+filename)
FILE = open(filename, "r")
index = 0
for lineItems in FILE:     #Read in file content
    #lineItems = lineItems.rstrip().split(",")
    lineItems = lineItems.rstrip().split()
    if (len(lineItems)>1):
		time.append(float(lineItems[0]))
    else:
		time.append(index)
    index=index+1
    data.append(float(lineItems[-1]))
FILE.close() #EOF, close the FILE!


data = np.array(data)
N = len(data)
#time = np.arange(0,N,1)


filename = './fitted.dat'
print("Loading "+filename)
FILE = open(filename, "r")
for lineItems in FILE:     #Read in file content
    lineItems = lineItems.rstrip().split(",")
    fitted.append(float(lineItems[-1]))
FILE.close() #EOF, close the FILE!


filename = './residuals.dat'
print("Loading "+filename)
FILE = open(filename, "r")
for lineItems in FILE:     #Read in file content
    lineItems = lineItems.rstrip().split(",")
    residuals.append(float(lineItems[-1]))
FILE.close() #EOF, close the FILE!



# plot the real data
plt.figure(figsize = (15, 10))
plt.subplot(2, 1, 1)
plt.plot(time, data, 'r.', label = 'Real Values')
plt.plot(time, fitted, 'g', label = 'Fitted Sinewave')
plt.title('Sine fitting results')
plt.xlabel('Samples')
plt.ylabel('Sampled data')
plt.legend(shadow = True, fancybox = True)
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(time, residuals, 'r.-', label = 'Residual')
plt.xlabel('Samples')
plt.ylabel('Residuals')
plt.legend(shadow = True, fancybox = True)
plt.grid(True)
plt.show()
