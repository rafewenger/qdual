#!/usr/bin/python

import os
import sys

offFile = 'out_tqrs.off';
numTest = 10;
axisSize = 10;
maxVal = 5;
isovalue = 2.1;
seed0 = 12345;
flag_collapseB = False;
flag_collapseC = False;
flag_trimesh = True;
epsilon = 0.3;

# Function get_flag
def get_flag(option, option_list):
    flag = False;
    argument = None;
    for i in range(len(option_list)) :
        if (option_list[i] == option):
            flag = True;
    return(flag);

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

    qdual_command_line = 'qdual -s -move_vertex -epsilon ' + str(epsilon);
    
    if (flag_trimesh) :
        qdual_command_line = qdual_command_line + ' -trimesh';

    if (flag_collapseB):
        qdual_command_line = qdual_command_line + ' -collapseB';

    if (flag_collapseC):
        qdual_command_line = qdual_command_line + ' -collapseC';

    qdual_command_line = qdual_command_line + \
                         ' -o ' + offFile + ' ' + str(isovalue) + \
                         ' random.nrrd';
    
    print 'Executing: ', qdual_command_line;
    sys.stdout.flush()
    ret_code = os.system(qdual_command_line);
    if (ret_code != 0):
        print 'qdual return code: ', ret_code
        print 'Exiting ', sys.argv[0]
        exit(ret_code)

    if (flag_trimesh):
        ijkmeshinfo_command_line = 'ijkmeshinfo -mesh_dim 2 -manifold -terse -report_deep ' + offFile;
        os.system(ijkmeshinfo_command_line);

    ijkmeshinfo_command_line = 'ijkmeshinfo -mesh_dim 2 -internal -out_min_angle ' + offFile;
    os.system(ijkmeshinfo_command_line);

    ijkmeshinfo_command_line = 'ijkmeshinfo -mesh_dim 2 -out_min_angle ' \
                               + offFile;
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

flagFound, argument = get_argument('-epsilon', command_options);
if (flagFound):
    epsilon = float(argument);
    command_options = strip_option('-epsilon', 1, command_options);

flag_collapseB = get_flag('-collapseB', command_options);
if (flag_collapseB):
    command_options = strip_option('-collapseB', 0, command_options);

flag_collapseC = get_flag('-collapseC', command_options);
if (flag_collapseC):
    command_options = strip_option('-collapseC', 0, command_options);

flagFound = get_flag('-quad', command_options);
if (flagFound):
    flag_trimesh = False;
    command_options = strip_option('-quad', 0, command_options);


seed = seed0;
for i in range(numTest):
    runqdual();
    seed = seed+10;
