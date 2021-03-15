#!/usr/bin/env python
#-*- coding: utf-8 -*-
import os
import glob


for i in glob.iglob('*.bam'):
  print('ScanExitron.py -i {0} -r hg38'.format(i))
