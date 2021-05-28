#!/usr/bin/env python
#-*- coding: utf-8 -*-
#===============================================================================
import sys
import re
import os
import glob

def process(infile, outfile, normals):
    out = open(outfile, 'w')
    with open(infile) as f:
        header = f.readline()
        out.write(header)
        for line in f:
            l = line.rstrip().split('\t')
            key = '{}\t{}\t{}\t{}'.format(l[0], l[1], l[2], l[5])
            if not key in normals:
                out.write(line)
    out.close()

normals = set([])
with open('../normals/normal.stats.txt') as f:
    f.readline()
    for line in f:
        l = line.rstrip().split('\t')
        if int(l[-1]) > 3:
            normals.add('{}'.format('\t'.join(l[:-1])))

for i in glob.iglob('../tumors/*.exitron'):
    name = os.path.basename(i).split('.')[0]
    process(i, f'{name}.exitron', normals)
