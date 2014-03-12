#!/usr/bin/python

import os
import sys

offFile = 'out.off'
numTest = 10
axisSize = 20
maxVal = 5
isovalue = 2.1
seed = 12345
epsilon = 0.3


# Function runqdual
def runqdual():

    ijkgenscalar_command_line = 'ijkgenscalar -s -dim 3 -asize ' + str(axisSize) + ' -field randomint -bzero -seed ' + str(seed) + ' -maxval ' + str(maxVal) + ' random.nrrd';

    print 'Executing: ', ijkgenscalar_command_line;
    sys.stdout.flush()
    os.system(ijkgenscalar_command_line);

    qdual_command_line = 'qdual -s -trimesh -move_vertex ' + \
                         '-epsilon ' + str(epsilon) + \
                         ' -o ' + offFile + ' ' + str(isovalue) + \
                         ' random.nrrd';
    
    print 'Executing: ', qdual_command_line;
    sys.stdout.flush()
    ret_code = os.system(qdual_command_line);
    if (ret_code != 0):
        print 'qdual return code: ', ret_code
        print 'Exiting ', sys.argv[0]
        exit(ret_code)
    
    ijkmeshinfo_command_line = 'ijkmeshinfo -manifold -terse out.off';
    os.system(ijkmeshinfo_command_line);

    ijkmeshinfo_command_line = 'ijkmeshinfo -out_min_angle out.off';
    os.system(ijkmeshinfo_command_line);

    return;

# main routine
if (len(sys.argv) > 1) :
    if (not((sys.argv[1]).isdigit())) :
        print 'Usage: ', sys.argv[0], ' [ <numTest> [ <seed> [ <maxval> ]]]';
        exit(10);
        
if (len(sys.argv) > 1) :
    numTest = int(sys.argv[1]);

if (len(sys.argv) > 2) :
    seed = int(sys.argv[2]);

if (len(sys.argv) > 3) :
    maxVal = int(sys.argv[3]);
    isovalue = maxVal/2.0 + 0.1;

for i in range(numTest):
    runqdual();
    seed = seed+10;


