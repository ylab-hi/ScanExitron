#!/usr/bin/env python
#-*- coding: utf-8 -*-
#===============================================================================
import sys
import re
import os
import glob
from collections import defaultdict

count = defaultdict(int)

for i in glob.iglob('*.exitron'):
    barcode = os.path.basename(i).split('.')[0]
    with open(i) as f:
        f.readline()
        for line in f:
            l = line.rstrip().split('\t')
            key = f'{l[0]}\t{l[1]}\t{l[2]}\t{l[5]}'
            count[key] += 1

out = open('normal.stats.txt','w')
out.write('chrom\tstart\tend\tstrand\n')
for i in count:
    out.write(f'{i}\t{count[i]}\n')
out.close()
