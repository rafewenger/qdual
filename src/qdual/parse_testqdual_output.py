#!/usr/bin/python

import re
import os
import sys

# global variables
separator = ' '

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


# main routine
if (len(sys.argv) != 3) :
    print 'Usage error.  Usage: ', sys.argv[0], " <infile> <outfile>"
    exit(10)

inFilename = sys.argv[1]
outFilename = sys.argv[2]
infile = open(inFilename, 'r')
outfile = open(outFilename, 'w')

# defaults
ijkgenscalar_seed = 0;
qdual_seed = 0;

header = 'fileName' + separator + 'isovalue' + \
         separator + 'genscalarSeed' + \
         separator + 'qdualSeed' + \
         separator + 'randomNumI' + \
         separator + 'minAngle' + \
         separator + 'Manifold' + \
         separator + 'Boundary' + '\n'
    
outfile.write(header);


for line in infile:
    matchObj = re.match('Executing...qdual', line);
    if (matchObj):
        qdual_seed = getArgument('-seed', line, 0);
        random_num_intervals = getArgument('-random_num_intervals', line, 1);
        words = line.split();
        qdualFile = words[-1];
        isovalue = words[-2];

    matchObj = re.match('Executing...ijkgenscalar', line);
    if (matchObj) :
        ijkgenscalar_seed = getArgument('-seed', line, 0);
        words = line.split();

    flag_passed_manifold = True;
    flag_passed_boundary = True;
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

    flagM = 'F';
    if (flag_passed_manifold) :
        flagM = 'P';

    flagB = 'F';
    if (flag_passed_boundary) :
        flagB = 'P';
    
    matchObj = re.match('Min polygon angle', line)
    if (matchObj) :
        words = line.split()
        min_angle = words[-1]
        s = str(qdualFile) + separator + str(isovalue) + \
            separator + str(ijkgenscalar_seed) + \
            separator + str(qdual_seed) + \
            separator + str(random_num_intervals) + \
            separator + str(min_angle) + \
            separator + flagM + separator + flagB + '\n'
        outfile.write(s)
    

    

