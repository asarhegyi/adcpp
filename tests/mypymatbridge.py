#!/usr/bin/env python

'''
    This script tests the python-matlab-bridge setup.
    First it searches for a matlab executable. If matlab
    is on the path it creates a matlab session and connects
    the python interpreter. After that it calls jk.m mcruns
    times which emulates a counter. Once it is done with
    counting it shuts down the matlab server.   
'''

from pymatbridge import Matlab
from distutils.spawn import find_executable


mcruns = 800
res = {'result': 0}

# Search for matlab	
prog = find_executable('matlab')
if not prog:
    sys.exit("Could not find MATLAB on PATH. Exiting...")
else:
    print("Found MATLAB in "+prog) 

# Create matlab session
mlab = Matlab(executable=prog)
mlab.start()

# Call jk.m in a for loop to emulate a counter
for mcindex in range(mcruns):
    res = mlab.run_func('../matlab/adder.m', {'arg1': res['result'], 'arg2': 1})
    print(res['result'])

# Shut down the matlab server
mlab.stop()

