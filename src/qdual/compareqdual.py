#!/usr/bin/python

import os
import sys

offFile1 = 'out1.off'
offFile2 = 'out2.off'
numTest = 5
axisSize = 20
maxVal = 255
isovalue = 127.5
seed = 12345

# Function runqdual
def runqdual(qdual1, qdual2):

    ijkgenscalar_command_line = 'ijkgenscalar -s -dim 3 -asize ' + str(axisSize) + ' -field randomint -bzero -seed ' + str(seed) + ' -maxval ' + str(maxVal) + ' random.nrrd';

    print 'Executing: ', ijkgenscalar_command_line;
    os.system(ijkgenscalar_command_line);

    qdual_options = '-s -trimesh -move_vertex -collapseC -del_isolate ' + str(isovalue) + ' random.nrrd';
    qdual_command_line = qdual1 + ' -o ' + offFile1 + ' ' + qdual_options;
    print 'Executing: ', qdual_command_line;
    os.system(qdual_command_line);

    qdual_command_line = qdual2 + ' -o ' + offFile2 + ' ' + qdual_options;
    print 'Executing: ', qdual_command_line;
    os.system(qdual_command_line);

    diff_command_line = 'diff -q ' + offFile1 + ' ' + offFile2
    os.system(diff_command_line)

    print ''
    
    return;

# main routine
if (len(sys.argv) != 3) :
    print 'Usage error.  Usage: ', sys.argv[0], " <qdual1> <qdual2>"
    exit(10)


for i in range(numTest):
    seed = seed+10*i;
    runqdual(sys.argv[1], sys.argv[2]);


