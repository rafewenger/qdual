#!/usr/bin/python
# test qdual file using random positions

import os
import sys

offFile = 'out.off';
numTest = 10;
seed = 7654;
flagNoCollapse = False;
flagMoveVertex = True;

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

    qdual_command_line = 'qdual -position random -seed ' + str(seed) + \
                         ' ' + ' '.join(qdual_options)
    
    print 'Executing: ', qdual_command_line;
    sys.stdout.flush()
    os.system(qdual_command_line);
    ret_code = os.system(qdual_command_line);
    if (ret_code != 0):
        print 'qdual return code: ', ret_code
        print 'Exiting ', sys.argv[0]
        exit(ret_code)

    ijkmeshinfo_command_line = 'ijkmeshinfo -manifold -report_deep -terse out.off';
    os.system(ijkmeshinfo_command_line);

    ijkmeshinfo_command_line = 'ijkmeshinfo -out_min_angle out.off';
    os.system(ijkmeshinfo_command_line);

    return;

# main routine

if (len(sys.argv) <= 2) :
    print 'Usage error.  Usage: ', sys.argv[0], " [qdual options] <isovalue> {nrrdFile}"
    exit(10)

command_options = sys.argv[1:];
flagFound, argument = get_argument('-seed', command_options);
if (flagFound):
    seed = int(argument);
    command_options = strip_option('-seed', 1, command_options);

flagFound, argument = get_argument('-numTest', command_options);
if (flagFound):
    numTest = int(argument);
    command_options = strip_option('-numTest', 1, command_options);

if ('-no_collapse' in command_options):
    flagNoCollapse = True;
    flagMoveVertex = False;

qdual_options = ['-s', '-trimesh', '-o', offFile ];
if (flagMoveVertex):
    qdual_options.append('-move_vertex');
qdual_options = qdual_options + command_options;

for i in range(numTest):
    runqdual(qdual_options);
    seed = seed+10;
