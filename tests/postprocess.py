#!/usr/bin/env python

"""
    Plot the simulation results of the Monte Carlo regression
"""

__author__      = "Attila Sarhegyi"
__copyright__   = "Copyright 2009-2015, Planet Earth"

import os
import sys
import glob
import shutil

import testbench

args = sys.argv[1:]

#############################################################
############### find stimulus.dat and set cwd ###############
#############################################################
#path = testbench.find("stimulus.dat", "./")
#print("Found stimulus.dat in "+path+"\n")

path = sys.argv[1]

os.chdir(path)
print("Current working directory: "+os.getcwd()+"\n")


#############################################################
##################### Link source files #####################
#############################################################
print("\nLinking source files:\n")

linkList = ["../../matlab/plot_mchist.m",
            "../../matlab/plot_deltahist.m"]

testbench.symlink_srcfiles(linkList)


#############################################################
################### Plot the histograms #####################
#############################################################
# Post-Process Monte Carlo simulation results and plot histograms
# for all of the sfit algorithms
testbench.plot_mchist("sfit3_matlab.est", "subplot")
testbench.plot_mchist("sfit3_linear_cpp.est", "subplot")
testbench.plot_mchist("sfit3_large_cpp.est", "subplot")

testbench.plot_mchist("sfit4_matlab.est", "subplot")
testbench.plot_mchist("sfit4_linear_cpp.est", "subplot")
testbench.plot_mchist("sfit4_large_cpp.est", "subplot")


# Analyze the delta between the algorithms by plotting
# the histograms of the differences.
testbench.plot_deltahist("sfit3_linear_cpp.est", "sfit3_matlab.est")
testbench.plot_deltahist("sfit3_large_cpp.est", "sfit3_matlab.est")
testbench.plot_deltahist("sfit3_linear_cpp.est", "sfit3_large_cpp.est")

testbench.plot_deltahist("sfit4_linear_cpp.est", "sfit4_matlab.est")
testbench.plot_deltahist("sfit4_large_cpp.est", "sfit4_matlab.est")
testbench.plot_deltahist("sfit4_linear_cpp.est", "sfit4_large_cpp.est")


# Move figures into subdirectories
os.mkdir("EPS")
for file in glob.glob('*.eps'):
    shutil.move(file, "EPS/")

os.mkdir("PNG")
for file in glob.glob('*.png'):
    shutil.move(file, "PNG/")

# Clean up the current working directory
if not 'debug' in args:
    # delete symlinks        
    testbench.cleanup(linkList)
