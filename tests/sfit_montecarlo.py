#!/usr/bin/env python

"""
    This script runs a Monte Carlo simulation on sfit3_linear,
    sfit3_large, sfit4_linear, and sfit4_large algorithms
    from the ADC++ code base and compares the results with
    the MATLAB version of the sfit3/4 algorithms
"""

__author__      = "Attila Sarhegyi"
__copyright__   = "Copyright 2009-2015, Planet Earth"


import sys
import subprocess
import numpy as np
import os, glob, shutil
import matplotlib.pyplot as plt
from collections import namedtuple

import testbench

args = sys.argv[1:]

#############################################################
##################### Link source files #####################
#############################################################
print("\nLinking source files:\n")

linkList = ["../src/stimulus.h",
            "../src/sfit_class.h",
            "../examples/sfit3_linear.cpp",
            "../examples/sfit3_large.cpp",
            "../examples/sfit4_linear.cpp",
            "../examples/sfit4_large.cpp",
            "../matlab/sfit3.m",
            "../matlab/sfit4.m",
            "../matlab/sfit4imp.m",
            "../matlab/wrapper_sfit.m"]

if ('debug_resolution' in args):
    linkList.append("../examples/debug_resolution.cpp")
else:
    linkList.append("../examples/gen_rand_sine.cpp")

testbench.symlink_srcfiles(linkList)


#############################################################
################# Compile C++ source files ##################
#############################################################
print("\nCompiling source files:\n")

ccList = [ "sfit3_linear",
            "sfit3_large",
            "sfit4_linear",
            "sfit4_large"]

if ('debug_resolution' in args):
    ccList.append("debug_resolution")
else:
    ccList.append("gen_rand_sine")

testbench.cppCompile(ccList)


#############################################################
################ Run Monte Carlo simulations ################
#############################################################

sfit3_linlog = "sfit3_linear_cpp.est"
testbench.add_loghead(sfit3_linlog, "a");

sfit3_lrglog = "sfit3_large_cpp.est"
testbench.add_loghead(sfit3_lrglog, "a");

sfit3_matlog = "sfit3_matlab.est"
testbench.add_loghead(sfit3_matlog, "a");

sfit4_linlog = "sfit4_linear_cpp.est"
testbench.add_loghead(sfit4_linlog, "a");

sfit4_lrglog = "sfit4_large_cpp.est"
testbench.add_loghead(sfit4_lrglog, "a");

sfit4_matlog = "sfit4_matlab.est"
testbench.add_loghead(sfit4_matlog, "a");


# Number of Monte Carlo runs
if ('debug_resolution' in args):
    mcruns = 101
elif ('single' in args):
    mcruns = 1
else:
    mcruns = 1000

# Start matlab session
mlab = testbench.start_matlab()

# Run testbench
for mcindex in range(mcruns):

    if ('debug_resolution' in args):
        vmax = 1.0+mcindex*0.1;     # 1:0.1:11
        #vmax = 2.0+mcindex*0.02;    # 2:0.02:4


    #############################################################
    ################ Generate sine-wave stimulus ################
    #############################################################
    print("\nGenerating sine-wave stimulus:\n")

    if ('debug_resolution' in args):
        cmd="./debug_resolution -v"+'{0:<4.2f}'.format(vmax)
    else:
        cmd="./gen_rand_sine"

    subprocess.check_call(["echo", "-e", cmd+"\n"])
    subprocess.check_call(cmd, shell=True)

    #############################################################
    ################### Parse parameter file ####################
    #############################################################
    #Read in the parameters from file
    filename = './parameters.dat'
    print("Loading "+filename)
    FILE = open(filename, "r")
    for lineItems in FILE:
        lineItems = lineItems.rstrip().split(" ")
    FILE.close()

    # Print parameters into a summary file
    fname_stimulus = "stimulus.dat"
    FSTOUT = open(fname_stimulus, "a")
    for indx, param in enumerate(lineItems):
        if (indx == 0) | (indx == 2) | (indx == 3) | (indx == 7):
            FSTOUT.write('{0:> 24.20f}'.format(float(param)))
        elif (indx == 1) | (indx == 4):
            FSTOUT.write('{0:> 27.20f}'.format(float(param)))
        elif (indx == 5):
            FSTOUT.write('{0:> 8d}'.format(int(param)))
        elif (indx == 6):
            FSTOUT.write('{0:> 4d}'.format(int(param)))
        elif (indx == 8) | (indx == 9):
            FSTOUT.write('{0:> 10.6f}'.format(float(param)))
    FSTOUT.write('\n')
    FSTOUT.close()

    # Create a new tuple subclass named Params
    Params = namedtuple("Params", "amplitude frequency phase dc fs nofSamples qbits Q Vmax Vmin") 

    # instantiate the subclass as Stimulus
    Stimulus = Params(lineItems[0], lineItems[1], lineItems[2], lineItems[3], lineItems[4], lineItems[5], lineItems[6], lineItems[7], lineItems[8], lineItems[9])

    harmonics = 0
    inputFile = "./quantizedData.dat"

    #############################################################
    ###################### Run sfit3_linear #####################
    #############################################################
    stdout_value = testbench.run_sfitcpp(3, "linear", Stimulus, harmonics, inputFile)

    # parse stdout and append results to the end of the log file
    testbench.parse_stdout(stdout_value, sfit3_linlog, "a")

    #############################################################
    ###################### Run sfit3_large ######################
    #############################################################
    stdout_value = testbench.run_sfitcpp(3, "large", Stimulus, harmonics, inputFile)

    # parse stdout and append results to the end of the log file
    testbench.parse_stdout(stdout_value, sfit3_lrglog, "a")

    #############################################################
    #################### Run MATLAB : sfit3.m ###################
    #############################################################
    matlab_out = testbench.run_sfitmat(mlab, 3, Stimulus, inputFile)

    # Parse the matlab results and append it to the end of the log file
    testbench.parse_matout(matlab_out, sfit3_matlog, "a")

    #############################################################
    ###################### Run sfit4_linear #####################
    #############################################################
    stdout_value = testbench.run_sfitcpp(4, "linear", Stimulus, harmonics, inputFile)

    # parse stdout and append results to the end of the log file
    testbench.parse_stdout(stdout_value, sfit4_linlog, "a")

    #############################################################
    ###################### Run sfit4_large ######################
    #############################################################
    stdout_value = testbench.run_sfitcpp(4, "large", Stimulus, harmonics, inputFile)

    # parse stdout and append results to the end of the log file
    testbench.parse_stdout(stdout_value, sfit4_lrglog, "a")

    #############################################################
    #################### Run MATLAB : sfit4.m ###################
    #############################################################
    matlab_out = testbench.run_sfitmat(mlab, 4, Stimulus, inputFile)

    # Parse the matlab results and append it to the end of the log file
    testbench.parse_matout(matlab_out, sfit4_matlog, "a")


# Stop matlab session
mlab.stop()


# Clean up the current working directory
if not 'debug' in args:
    # delete symlinks        
    testbench.cleanup(linkList)

    # delete compiled files        
    testbench.cleanup(ccList)


# Move estimation results into a subdirectory for post-processing
os.mkdir("results")
for file in glob.glob('*.est'):
    shutil.move(file, "results/")

shutil.move("stimulus.dat", "results/")

