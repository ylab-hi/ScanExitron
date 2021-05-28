#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    output = open('run_ScanExitron.sh', 'w')
    with open('../PRAD.info.txt') as f:
        f.readline()
        for line in f:
            l = line.rstrip().split('\t')
            barcode = l[0]
            output.write(f'ScanExitron.py -i {barcode}.bam -m 50\n')
    output.close()
