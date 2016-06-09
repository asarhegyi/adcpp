#!/usr/bin/env python

"""
testbech.py
    This module defines functions to executes the 3- and
    4-parameter sine wave fitting algorithms written in C++
    and/or MATLAB.
    There are currently 4 sine fit algorithms implemented in
    C++ such as

        sfit3_linear -- 3-parameter sine fit algorithm for a
        small set of data ( < 32768 number of samples).

        sfit3_large -- a recursive version of the sfit3_linear
        algorithm to be able to process large data sets. It
        can process ten million samples with ease.

        sfit4_linear --  4-parameter sine fit algorithm for a
        small set of data ( < 32768 number of samples).

        sfit4_large -- a recursive version of the sfit4_linear
        algorithm to be able to process large data sets. It
        still can process ten million samples with ease but a
        little bit slower sfit3_large.

    The purpose if this module is to develop a standardized
    platform to simulate and compare the performance of the
    various algorithms. One way to do it is to run Monte Carlo
    simulations and do a statistical analysis.

    TODO: 
        -- recode plot_deltahist.m & plot_mchist.m ROI?
"""

__author__      = "Attila Sarhegyi"
__copyright__   = "Copyright 2009-2015, Planet Earth"


import re
import sys
import os.path
import subprocess
from pymatbridge import Matlab
from collections import namedtuple
from distutils.spawn import find_executable


#------------------------------------------------------------------------------
"""
run_sfitcpp
    This function runs the C++ version of the 3- and 4-parameter
    sine wave fitting algorithms such as  sfit3_linear,
    sfit3_large, sfit4_linear, and sfit4_large.

Parameters:
    numParams -- selects the 3- or 4-parameter fit
    sfitMode -- it can be "linear" or "large". See more above.
    Stimulus -- stores the parameters of the stimulus. Only the
        sampling and the fundamental frequencies are used
        in this function.
    numHarmonics -- number of harmonics to estimate besides the
        fundamental. If it is zero, only the fundamental component
        is estimated.
    inputFile -- the file that stores the observations
"""
def run_sfitcpp(numParams, sfitMode, Stimulus, numHarmonics, inputFile, *varags):

    if not varags:  # no blockSize define
        cmd1="./sfit"+str(numParams)+"_"+sfitMode+" -f"+Stimulus.frequency+" -s"+Stimulus.fs+" -r"+str(numHarmonics)+" "+inputFile
    else:       # blockSize defined explicitly
        cmd1="./sfit"+str(numParams)+"_"+sfitMode+" -f"+Stimulus.frequency+" -s"+Stimulus.fs+" -r"+str(numHarmonics)+" -b"+str(varags[-1])+" "+inputFile

    print("\n"+cmd1+"\n")
    sys.stdout.flush()

    p1=subprocess.Popen(cmd1,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True)

    p1.wait()    # wait to finish the subprocess
    p1_stdout = p1.communicate()[0]
    print(p1_stdout)
    sys.stdout.flush()

    return p1_stdout


#------------------------------------------------------------------------------
"""
start_matlab
    Invoke matlab, connect it to python and return with a
    matlab session class instance
"""
def start_matlab():

    prog = find_executable('matlab')
    if not prog:
        sys.exit("Could not find MATLAB on the PATH. Exiting...")
    else:
        print("Found MATLAB in "+prog) 

    mlab = Matlab(executable=prog)
    mlab.start()

    return mlab


#------------------------------------------------------------------------------
"""
run_sfitmat
    This function calls sfit3.m or sfit4.m though their wrapper
    scripts (wrapper_sfit3.m, wrapper_sfit4.m), prints results
    to stdout and returns with a results object.

Parameters:
    mlab -- the handler the open matlab session
    numParams -- selects the 3- or 4-parameter fit
    Stimulus -- stores the parameters of the stimulus. Only the
        sampling and the fundamental frequencies are used
        in this function.
    inputFile -- the file that stores the observations
"""
def run_sfitmat(mlab, numParams, Stimulus, inputFile):

    results = mlab.run_func('wrapper_sfit.m', {\
            'arg1': str(inputFile),\
            'arg2': int(numParams),\
            'arg3': float(Stimulus.fs),\
            'arg4': float(Stimulus.frequency)\
            })

    print('\nMATLAB results')
    print('Parameters of the sine-wave:')
    print("\tAmp0  : %.20e" %float(results['result']['A']))
    print("\tPhase0: %.20e" %float(results['result']['phi']))
    print("\tFreq  : %.20e" %float(results['result']['f']))
    print("\tdc    : %.20e" %float(results['result']['dc']))

    print("\nfs: %.20f" %float(results['result']['sample_rate']))
    print("Execution time: %.6f" %float(results['result']['exec_time']))

    print("\n\terms: %.20e\n" %float(results['result']['erms']))

    return results


#------------------------------------------------------------------------------
"""
add_logheads
    This function adds the parameter header line to the log files
    to make it more user friendly.

Parameters:
    filename -- name of the log file
    mode -- describes the way in which the file will be opened.
        'w' open for writing only; 'a' open for appending; etc.
"""
def add_loghead(filename, mode):

    FILEOUT = open(filename, mode)
    FILEOUT.write('%')
    for param in ['Amplitude', 'Frequency', 'Phase', 'dc', 'erms']:
        FILEOUT.write('{0:^27s}'.format(param))
    FILEOUT.write('{0:^10s}\n'.format('rtime'))
    FILEOUT.close()


#------------------------------------------------------------------------------
"""
parse_stdout
    This function parses the standard output and looks for the
    estimated parameters of the sine fit algorithms. It is used
    only for parsing the output of the C++ code. Results are
    stored in a log file.

Parameters:
    stdout_str -- data stream from the standard output
    filename -- name of the log file
    mode -- describes the way in which the file will be opened.
        'w' open for writing only; 'a' open for appending; etc.
"""
def parse_stdout(stdout_str, filename, mode):

    Amp = []
    Phase = []
    Freq = []
    dc = []
    erms = []
    fs = []
    rtime = []

    match = re.search(r'Amp0\s*:\s*([+-]?\d*.\d+e[+-]?\d+)', stdout_str)
    Amp.append(float(match.group(1)))

    match = re.search(r'Phase0\s*:\s*([+-]?\d*.\d+e[+-]?\d+)', stdout_str)
    Phase.append(float(match.group(1)))

    match = re.search(r'Freq\s*:\s*([+-]?\d*.\d+e[+-]?\d+)', stdout_str)
    Freq.append(float(match.group(1)))

    match = re.search(r'dc\s*:\s*([+-]?\d*.\d+e[+-]?\d+)', stdout_str)
    dc.append(float(match.group(1)))

    match = re.search(r'erms\s*:\s*([+-]?\d*.\d+e[+-]?\d+)', stdout_str)
    erms.append(float(match.group(1)))

    match = re.search(r'fs\s*:\s*(\d*.\d+)', stdout_str)
    fs.append(float(match.group(1)))

    match = re.search(r'Execution time\s*:\s*(\d*.\d+) seconds', stdout_str)
    rtime.append(float(match.group(1)))

    # Print coeffs into a log file
    FILEOUT = open(filename, mode)
    for param in [Amp[-1], Freq[-1], Phase[-1], dc[-1], erms[-1]]:
        FILEOUT.write('{0:> 27.20f}'.format(param))
    FILEOUT.write('{0:> 10.6f}\n'.format(rtime[-1]))
    FILEOUT.close()

    # Create a new tuple subclass named Results
    Results = namedtuple("Results", "amplitude phase frequency dc erms fs runtime")

    # Instantiate the subclass as cppResults
    cppResults = Results(Amp[-1], Phase[-1], Freq[-1], dc[-1], erms[-1], fs[-1], rtime[-1])

    return cppResults


#------------------------------------------------------------------------------
"""
parse_matout
    This function parses the returned data structure from the
    matlab function calls (sfit3.m, sfit4.m). Results are
    stored in a log file.

Parameters:
    matlab_out -- data structure returned from the matlab call
    filename -- name of the log file
    mode -- describes the way in which the file will be opened.
        'w' open for writing only; 'a' open for appending; etc.
"""
def parse_matout(matlab_out, filename, mode):

    Amp = []
    Phase = []
    Freq = []
    dc = []
    erms = []
    fs = []
    rtime = []

    Amp.append(float(matlab_out['result']['A']))

    Phase.append(float(matlab_out['result']['phi']))

    Freq.append(float(matlab_out['result']['f']))

    dc.append(float(matlab_out['result']['dc']))

    erms.append(float(matlab_out['result']['erms']))

    fs.append(float(matlab_out['result']['sample_rate']))

    rtime.append(float(matlab_out['result']['exec_time']))

    FILEOUT = open(filename, mode)
    for param in [Amp[-1], Freq[-1], Phase[-1], dc[-1], erms[-1]]:
        FILEOUT.write('{0:> 27.20f}'.format(param))
    FILEOUT.write('{0:> 10.6f}\n'.format(rtime[-1]))
    FILEOUT.close()

    # Create a new tuple subclass named Results
    Results = namedtuple("Results", "amplitude phase frequency dc erms fs runtime")

    # Instantiate the subclass as matlabResults
    matlabResults = Results(Amp[-1], Phase[-1], Freq[-1], dc[-1], erms[-1], fs[-1], rtime[-1])

    return matlabResults


#------------------------------------------------------------------------------
"""
symlink_srcfiles
    Symlink all source files defined on the linkList

Parameters:
    linkList -- list of the source files with their relative path
"""
def symlink_srcfiles(linkList):

    for linkItem in linkList:
        splitItem = linkItem.rstrip().split("/")
        if os.path.lexists(splitItem[-1]):
            print("Symlink to "+splitItem[-1]+" exists")
        else:
            cmd="ln -s "+linkItem
            subprocess.check_call(["echo", cmd])
            subprocess.check_call(cmd, shell=True)


#------------------------------------------------------------------------------
"""
cppCompile
    Compile C++ source files

Parameters:
    compileList -- list of the source files to compile
"""
def cppCompile(compileList):

    for compItem in compileList:
        cmd="make -f ../Makefile TARGET="+compItem
        subprocess.check_call(["echo", cmd])
        subprocess.check_call(cmd, shell=True)


#------------------------------------------------------------------------------
"""
cleanup
    Delete files and symlinks in the current directory

Parameters:
    deleteList -- list of the files and symlinks to delete
"""
def cleanup(deleteList):

    for deleteItem in deleteList:
        splitItem = deleteItem.rstrip().split("/")
        if os.path.lexists(splitItem[-1]):
            print("Symlink "+splitItem[-1]+" exists in the current directory")
            cmd="rm "+splitItem[-1]
        elif os.path.exists(splitItem[-1]):
            print("File "+splitItem[-1]+" exists in the current directory")
            cmd="rm "+splitItem[-1]
        else:
            print(splitItem[-1]+" does not exist in the current directory") 
            cmd="echo"

        subprocess.check_call(["echo", cmd])
        subprocess.check_call(cmd, shell=True)


#------------------------------------------------------------------------------
"""
plot_mchist
    Plot histograms of the Monte Carlo simulation results
"""
def plot_mchist(fileName, plotFormat):

    prog = find_executable('matlab')
    if not prog:
        sys.exit("Could not find MATLAB on the PATH. Exiting...")
    else:
        print("Found MATLAB in "+prog) 

    mlab = Matlab(executable=prog)
    mlab.start()
    results = mlab.run_func('./plot_mchist.m', {'arg1': fileName, 'arg2': plotFormat})
    mlab.stop()


#------------------------------------------------------------------------------
"""
plot_deltahist
    Plot the delta of the Monte Carlo simulation results
"""
def plot_deltahist(fileName1, fileName2):

    prog = find_executable('matlab')
    if not prog:
        sys.exit("Could not find MATLAB on the PATH. Exiting...")
    else:
        print("Found MATLAB in "+prog) 

    mlab = Matlab(executable=prog)
    mlab.start()
    results = mlab.run_func('./plot_deltahist.m', {'arg1': fileName1, 'arg2': fileName2})
    mlab.stop()


#------------------------------------------------------------------------------
"""
find
    Find file on the path
"""
def find(filename, path):
    for root, dirs, files in os.walk(path):
        if filename in files:
            #return os.path.join(root, name)
            return root


#------------------------------------------------------------------------------

if __name__ == "__main__":
    sys.exit("main has not defined yet")
