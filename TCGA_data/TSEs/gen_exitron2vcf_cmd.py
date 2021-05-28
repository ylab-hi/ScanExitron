#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
import glob
out = open('exitron2vcf.sh', 'w')
for i in glob.iglob('*.exitron'):
    name = os.path.basename(i).split('.')[0]
    out.write(f'exitron2vcf.py -i {i} -o ./VCFs/{name}.vcf\n')
out.close()
