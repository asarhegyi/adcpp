#!/usr/bin/env python

"""
    This script intends to verify and visualize the 
    waveforms generated with the stimulus class.
"""

__author__      = "Attila Sarhegyi"
__copyright__   = "Copyright 2009-2015, Planet Earth"


import sys
import subprocess
import os, glob, shutil
from pymatbridge import Matlab

import testbench

args = sys.argv[1:]

nofBits=sys.argv[1]

#############################################################
##################### Link source files #####################
#############################################################
print("\nLinking source files:\n")

linkList = ["../matlab/plot_stimulus.m",
            "../matlab/quantizer.m"]

testbench.symlink_srcfiles(linkList)

#############################################################
################# Compile C++ source files ##################
#############################################################
print("\nCompiling source files:\n")

ccList = [ "gen_stimulus"]

testbench.cppCompile(ccList)


#############################################################
##################### Run test_stimulus #####################
#############################################################
print("\nRunning test_stimulus:\n")

cmd="./gen_stimulus -n"+'{0:<d}'.format(int(nofBits))

subprocess.check_call(["echo", "-e", cmd+"\n"])
subprocess.check_call(cmd, shell=True)


#############################################################
######################### Run MATLAB ########################
#############################################################
# Start matlab session
mlab = testbench.start_matlab()

# Call the gen_stimulus_plots MATLAB function to plot the waveforms.
result1 = mlab.run_func('plot_stimulus.m', {'arg1': int(nofBits)})

#Call the quantizer MATLAB function to generate quantizedData and ideat Tl in MATLAB
result2 = mlab.run_func('quantizer.m', {'arg1': int(nofBits)})

# Stop matlab session
mlab.stop()



# Clean up the current working directory
if not 'debug' in args:
    # delete symlinks        
    testbench.cleanup(linkList)

    # delete compiled files        
    testbench.cleanup(ccList)


# Move results into a subdirectory
os.mkdir("results")
for file in glob.glob('*.dat'):
    shutil.move(file, "results/")

for file in glob.glob('*.eps'):
    shutil.move(file, "results/")
os.chdir("./results")


#############################################################
######################## Show results #######################
#############################################################

cmd="ghostview noisy_signal_path.eps &"
subprocess.check_call(["echo", "-e", "\n"+cmd+"\n"])
subprocess.check_call(cmd, shell=True)

cmd="ghostview noisy_ctl_path.eps &"
subprocess.check_call(["echo", "-e", "\n"+cmd+"\n"])
subprocess.check_call(cmd, shell=True)

cmd="tkdiff Tl.dat Tl.mat.dat &"
subprocess.check_call(["echo", "-e", "\n"+cmd+"\n"])
subprocess.check_call(cmd, shell=True)

cmd="tkdiff quantizedData.dat quantizedData.mat.dat &"
subprocess.check_call(["echo", "-e", "\n"+cmd+"\n"])
subprocess.check_call(cmd, shell=True)
