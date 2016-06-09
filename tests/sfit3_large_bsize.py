#!/usr/bin/env python

"""
    This script characterizes the execution
    time of the sfit3_large algorithm against blockSize.
"""

__author__      = "Attila Sarhegyi"
__copyright__   = "Copyright 2009-2015, Planet Earth"


import sys
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from collections import namedtuple

import testbench

args = sys.argv[1:]

#############################################################
##################### Link source files #####################
#############################################################
print("\nLinking source files:\n")

linkList = ["../examples/sfit3_large.cpp"]
testbench.symlink_srcfiles(linkList)


#############################################################
################# Compile C++ source files ##################
#############################################################
print("\nCompiling source files:\n")

ccList = ["sfit3_large"]
testbench.cppCompile(ccList)


#############################################################
# Run sfit3_large with diffrent block sizes, capture the
# execution time and the coefficients estimated by the LS method.
#############################################################
execTime = []   # Store average execution time for each block size
N = 100;        # Repeate measurement N times for each block size to get a distribution

# 
frequency = "191655.29098"
fSample   = "2500000000"
#inputFile = sys.argv[2]
inputFile = "../charData/vref_xtalk_short.txt"

# sfit3 parameters
harmonics = "0"
blockArray = ["524288", "262144", "131072", "65536", "32768", "16384", "4096", "1024", "256", "64", "16", "4"]

# Create a new tuple subclass named Params
Params = namedtuple("Params", "frequency fs")

# instantiate the subclass as Stimulus
Stimulus = Params(frequency, fSample)

# blockSize loop
for blockSize in blockArray:
    population = [] # repeat the same conversion N times and calculate the average

    # Save coefficients from each iteration for verification
    fname_coeffs = "coeffs_b"+blockSize+".est"
    testbench.add_loghead(fname_coeffs, "w")

    # repeat N times to get a distribution
    for x in range(N):
        # Run sfit3_large
        stdout_value = testbench.run_sfitcpp(3, "large", Stimulus, harmonics, inputFile, blockSize)

        # Parse stdout and append results to the end of the log file
        results = testbench.parse_stdout(stdout_value, fname_coeffs, "a")

        # Save execution time in an array 
        population.append(results.runtime) 

    execTime.append(np.sum(population)/N)
    del population[:]


#############################################################
# Save execution time results into a file
#############################################################
fname_time = "avg_exec_time_sfit3_large.txt"
FILEO = open(fname_time, "w")

FILEO.write("## Block Size\tExecution Time\n")
for index, item in enumerate(execTime):
    #FILEO.write("\t%s\t\t%f\n" %(blockArray[index], item))
    FILEO.write('{0:13d} {1:> 10.6f}\n'.format(int(blockArray[index]), item))
FILEO.close() #EOF, close the FILE!

# flush stdout
sys.stdout.flush()


#############################################################
# Plot execution time 
#############################################################
plt.plot(blockArray, execTime, 'r.-', label = 'Execution Time')
plt.xscale('log')
#plt.semilogx(blockArray, execTime, 'r.-', label = 'Execution Time')
plt.xlabel('Array Size')
plt.ylabel('Execution Time [s]')
plt.legend(shadow = True, fancybox = True)
plt.grid(True)
plt.show()


# delete symlinks from the current directory        
testbench.cleanup(linkList)
# delete compiled files        
testbench.cleanup(ccList)
