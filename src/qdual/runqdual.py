#!/usr/bin/python

import os
import sys

offFile = 'out.off'

if (len(sys.argv) <= 1) :
    print 'Usage error.  Usage: ', sys.argv[0], " [qdual options] <isovalue> {nrrdFile}"
    exit(10)

qdual_command_line = ' '.join(sys.argv[1:])
qdual_command_line = 'qdual -s -trimesh -move_vertex -o ' + offFile + ' ' + qdual_command_line

print 'Executing: ', qdual_command_line

os.system(qdual_command_line)

ijkmeshinfo_command_line = 'ijkmeshinfo -manifold -terse out.off'
os.system(ijkmeshinfo_command_line)

ijkmeshinfo_command_line = 'ijkmeshinfo -out_min_angle out.off'
os.system(ijkmeshinfo_command_line)




