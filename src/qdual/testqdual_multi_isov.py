#!/usr/bin/python

import os
import sys

offFile = 'out_tqmi.off';
numTest = 10;
isovalue0 = 0;
increment = 0.1;
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
def runqdual(qdual_options):

    qdual_command_line = 'qdual -s -trimesh -move_vertex ' + \
                         '-epsilon ' + str(epsilon) + \
                         ' ' + ' '.join(qdual_options) + \
                         ' -o ' + offFile + ' ' + str(isovalue) + ' ' + \
                         nrrd_file;
    
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
if (len(sys.argv) < 3) :
    print 'Usage: ', sys.argv[0], ' [-num_test <N>] [-increment <I>] [qdual options] <isovalue0> {nrrd_file}';
    exit(10);

command_options = sys.argv[1:-2];

nrrd_file = sys.argv[-1];
isovalue0 = float(sys.argv[-2]);

flagFound, argument = get_argument('-numTest', command_options);
if (flagFound):
    numTest = int(argument);
    command_options = strip_option('-numTest', 1, command_options);

flagFound, argument = get_argument('-increment', command_options);
if (flagFound):
    increment = float(argument);
    command_options = strip_option('-increment', 1, command_options);
    

isovalue = isovalue0;
for i in range(numTest):
    runqdual(command_options);
    isovalue = isovalue + increment;
