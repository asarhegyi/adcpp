#!/bin/bash

# convert GSL matrix dump into a matrix array

sed -i ':a;N;$!ba;s/\n/\t/g' m_X.dat

if [ $1 -eq 3 ]; then
    sed -i -r 's/(-?[-+0-9.e]+)\t(-?[-+0-9.e]+)\t(-?[-+0-9.e]+)\t/\1\t\2\t\3\n/g' m_X.dat
fi 

if [ $1 -eq 4 ]; then
    sed -i -r 's/(-?[-+0-9.e]+)\t(-?[-+0-9.e]+)\t(-?[-+0-9.e]+)\t(-?[-+0-9.e]+)\t/\1\t\2\t\3\t\4\n/g' m_X.dat
fi
