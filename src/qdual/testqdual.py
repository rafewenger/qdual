#!/usr/bin/python

import os
import sys

offFile = 'out.off'
numTest = 10
axisSize = 50
maxVal = 5
isovalue = 2.1
seed = 12345

# Function runqdual
def runqdual():

    ijkgenscalar_command_line = 'ijkgenscalar -s -dim 3 -asize ' + str(axisSize) + ' -field randomint -bzero -seed ' + str(seed) + ' -maxval ' + str(maxVal) + ' random.nrrd';

    print 'Executing: ', ijkgenscalar_command_line;
    os.system(ijkgenscalar_command_line);

    qdual_command_line = 'qdual -s -trimesh -move_vertex -o ' + offFile + ' ' + str(isovalue) + ' random.nrrd';
    
    print 'Executing: ', qdual_command_line;
    os.system(qdual_command_line);

    ijkmeshinfo_command_line = 'ijkmeshinfo -manifold -terse out.off';
    os.system(ijkmeshinfo_command_line);

    ijkmeshinfo_command_line = 'ijkmeshinfo -out_min_angle out.off';
    os.system(ijkmeshinfo_command_line);

    return;

# main routine
if (len(sys.argv) > 1) :
    numTest = sys.argv[1]

for i in range(numTest):
    seed = seed+10*i;
    runqdual();


