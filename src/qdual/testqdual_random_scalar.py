#!/usr/bin/python

import os
import sys

offFile = 'out_tqrs.off';
numTest = 10;
axisSize = 10;
maxVal = 5;
isovalue = 2.1;
seed0 = 12345;
epsilon = 0.3;

# Function get_argument
def get_argument(option, option_list):
    flag = False;
    argument = None;
    for i in range(len(option_list)-1) :
        if (option_list[i] == option):
            argument = option_list[i+1];
            flag = True;
    return(flag, argument);

# Remove option and arguments from command_line.
def strip_option(option, num_arg, option_list):
    new_option_list = [];
    i = 0;
    while (i < len(option_list)):
        if (option_list[i] != option):
            new_option_list.append(option_list[i]);
        else:
            i = i + num_arg;
        i = i+1;

    return(new_option_list);


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
    
    ijkmeshinfo_command_line = 'ijkmeshinfo -manifold -terse -report_deep ' + \
                               offFile;
    os.system(ijkmeshinfo_command_line);

    ijkmeshinfo_command_line = 'ijkmeshinfo -out_min_angle ' + offFile;
    os.system(ijkmeshinfo_command_line);

    return;

# main routine
command_options = sys.argv[1:];

flagFound, argument = get_argument('-numTest', command_options);
if (flagFound):
    numTest = int(argument);
    command_options = strip_option('-numTest', 1, command_options);

flagFound, argument = get_argument('-seed', command_options);
if (flagFound):
    seed0 = int(argument);
    command_options = strip_option('-seed', 1, command_options);

flagFound, argument = get_argument('-maxval', command_options);
if (flagFound):
    maxVal = int(argument);
    isovalue = maxVal/2.0 + 0.1;
    command_options = strip_option('-maxval', 1, command_options);

flagFound, argument = get_argument('-asize', command_options);
if (flagFound):
    axisSize = int(argument);
    command_options = strip_option('-asize', 1, command_options);


seed = seed0;
for i in range(numTest):
    runqdual();
    seed = seed+10;


