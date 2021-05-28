#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys
import re
import os
import glob

hla = {}
with open('PRAD.info.with.HLA.txt') as f:
    f.readline()
    for line in f:
        l = line.rstrip().split('\t')
        barcode = l[0]
        hla_type = l[-1]
        hla[barcode] = hla_type
        
out = open('run_ScanNeo_neoantigen.sh', 'w')
for i in glob.iglob('../VCFs/*.vep.vcf'):
    name = os.path.basename(i).split('.')[0]
    out.write(f'ScanNeo.py hla -t 16 -i {i} --alleles {hla[name]} --af PSO -e 9 -p /home/usr/IDEB/ -o {name}.tsv\n')
out.close()
