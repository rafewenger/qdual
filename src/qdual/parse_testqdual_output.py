#!/usr/bin/python

import re
import os
import sys

# global variables
separator = ' ';
flag_genscalar = False;
genscalar_header = [ 'fileName', 'axisSize', 'isovalue', 'epsilon',
                     'genscalarSeed', 'maxVal',
                     'minAngle', 'minInternalAngle', 'Manifold', 'Boundary'];
random_pos_header = [ 'fileName', 'isovalue', 'epsilon',
                      'qdualSeed', 'randomNumI',
                      'minAngle', 'minInternalAngle', 'Manifold', 'Boundary'];

# get argument
def getArgument(option, line, init_value):
    flag = False;
    words = line.split();
    argument = init_value;
    for i in range(len(words)-1) :
        if (words[i] == option):
            argument = words[i+1];
            flag = True;

    return(argument);

# search for string s in file
def searchFile(filename, s):

    infile = open(filename, 'r')

    flag = False;
    for line in infile:
        if "ijkgenscalar" in line:
            flag = True;

    infile.close();
    return(flag);


# Convert bool to character
def bool2char(flag, char_false, char_true):

    c = char_false;
    if (flag):
        c = char_true;

    return(c);


# main routine
if (len(sys.argv) != 3) :
    print 'Usage error.  Usage: ', sys.argv[0], " <infile> <outfile>"
    exit(10)

inFilename = sys.argv[1]
outFilename = sys.argv[2]

flag_genscalar = searchFile(inFilename, 'ijkgenscalar');

infile = open(inFilename, 'r')
outfile = open(outFilename, 'w')

# defaults
ijkgenscalar_seed = 0;
qdual_seed = 0;
axis_size = 0;
maxval = 0;
epsilon = '***';
flag_passed_manifold = False;
flag_passed_boundary = False;
min_angle = 180;
min_internal_angle = 180;

if (flag_genscalar) :
    header = separator.join(genscalar_header);
else:
    header = separator.join(random_pos_header);

header = header + '\n';
outfile.write(header);


for line in infile:
    matchObj = re.match('Executing...qdual', line);
    if (matchObj):
        qdual_seed = getArgument('-seed', line, 0);
        epsilon = getArgument('-epsilon', line, '***');
        random_num_intervals = getArgument('-random_num_intervals', line, 1);
        words = line.split();
        qdualFile = words[-1];
        isovalue = words[-2];
        min_angle = 180;
        min_internal_angle = 180;

    matchObj = re.match('Executing...ijkgenscalar', line);
    if (matchObj) :
        axis_size = getArgument('-asize', line, 0);
        ijkgenscalar_seed = getArgument('-seed', line, 0);
        maxval = getArgument('-maxval', line, 0);
        words = line.split();

    matchObj = re.match('Passed all', line);
    if (matchObj) :
        flag_passed_manifold = True;
        flag_passed_boundary = True;

    matchObj = re.match('Failed manifold', line);
    if (matchObj) :
        flag_passed_manifold = False;

    matchObj = re.match('Failed boundary', line);
    if (matchObj) :
        flag_passed_boundary = False;

    matchObj = re.match('Min internal polygon angle', line)
    if (matchObj) :
        words = line.split()
        min_internal_angle = words[-1]
        
    matchObj = re.match('Min polygon angle', line)
    if (matchObj) :
        words = line.split()
        min_angle = words[-1]
        flagM = bool2char(flag_passed_manifold, 'F', 'P');
        flagB = bool2char(flag_passed_boundary, 'F', 'P');
        if (flag_genscalar) :
            data_list = [ qdualFile, str(axis_size), str(isovalue), \
                          str(epsilon), \
                          str(ijkgenscalar_seed), str(maxval), \
                          min_angle, min_internal_angle,
                          flagM, flagB ];
        else:
            data_list = [ qdualFile, str(isovalue), \
                          str(epsilon), \
                          str(qdual_seed), str(random_num_intervals), \
                          min_angle, min_internal_angle,
                          flagM, flagB ];

        s = separator.join(data_list);
        s = s + '\n';

        outfile.write(s)

    

