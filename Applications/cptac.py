#!/usr/bin/env python
#-*- coding: utf-8 -*-
import sys
import re
import os
import argparse
import pandas as pd
import numpy as np
import glob
import random
import datetime
import subprocess
import shutil

def status_message(msg):
    print(msg)
    sys.stdout.flush()

def run_cmd(cmd, msg=None):
    status_message(cmd)
    if ',' in msg:
        begin, finish = msg.split(',')
        status_message(begin)
    else:
        finish = msg
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'Error happend!: {}\n{}'.format(err, err.output)
    else:
        error_msg = ''
    if not error_msg:
        status_message(finish)
        return True
    else:
        status_message(error_msg)
        return False

def get_spectra_list(folder):
    rand_num = random.getrandbits(50)
    tmpList = []
    for i in glob.iglob(os.path.abspath(folder) + '/*.mzML.gz'):
        tmpList.append('{}'.format(i))
    return tmpList

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")
    parser.add_argument('-i', '--input', action='store', dest='spectra', help="spcetra folder including (*.mzML.gz)", required=True)
    parser.add_argument('-n', '--name', action='store', dest='name', help="sample name", required=True)
    parser.add_argument('-f', '--fasta', action='store', dest='fasta', help="Input fasta (decoy sequence included)", required=True)
    parser.add_argument('-c', '--config', action='store', dest='config', help="output type (default: %(default)s)", default='./config/MSGF_iTRAQ.ini')
    parser.add_argument('-m', '--msgfplus', action='store', dest='msgfplus', help='MSGF+ executable (default: %(default)s)', default='/home/tywang/software/anaconda2/share/msgf_plus-2017.07.21-2/MSGFPlus.jar')
    parser.add_argument('-t', '--tmp', action='store', dest='tmp',help='temp directory (default: %(default)s)', default='/tmp/')
    parser.add_argument('-o', '--output', action='store', dest='output', help="Output directory (default: %(default)s)", default='output')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    memory=20000
    threads=16

    tmp_dir = args.tmp
    output_dir = args.output
    sample_name = args.name
    fasta = args.fasta
    msgfplus = args.msgfplus
    config = args.config

    src_spectra_list = get_spectra_list(args.spectra)
    if not os.path.isdir(output_dir):
        os.mkdir('{}'.format(output_dir))

    status_message('Processing {} begins'.format(sample_name))
    status_message('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    start_execute_timestamp = datetime.datetime.now()
    for spectra in src_spectra_list:
        rand = random.getrandbits(25)
        os.mkdir('{}'.format(os.path.join(tmp_dir, str(rand))))
        status_message('Processing {} begins'.format(spectra))
        status_message('##########################################')
        label = '{}.{}'.format(sample_name, os.path.basename(spectra).replace('.mzML.gz',''))

        idXML = '{}/{}.idXML'.format(output_dir, label)
        pi_idXML = '{}/{}_pi.idXML'.format(output_dir, label)
        # pi_idXML is newer than idXML
        if os.path.isfile(pi_idXML) and os.path.getsize(pi_idXML) > 0 and os.path.getmtime(pi_idXML) > os.path.getmtime(idXML):
            print('File {} exists and seems up-to-date. Skipping MSGFPlusAdapter and PeptideIndexer ...'.format(pi_idXML))
            continue
        if not os.path.isfile(idXML):
            spectra_in_tmp = os.path.join(tmp_dir, str(rand), label)
            unzip = 'zcat {} > {}.mzML'.format(spectra, spectra_in_tmp)
            run_cmd(unzip,'unzip {} accomplished!'.format(spectra))
            msgfplus_cmd = 'MSGFPlusAdapter -ini {} -in {}.mzML -out {} -database {} -executable {} -java_memory {} -threads {}'.format(config, spectra_in_tmp, idXML, fasta, msgfplus, memory, threads)
            run_cmd(msgfplus_cmd,'MSGFPlus {} accomplished!'.format(spectra_in_tmp))
        else:
            print('File {} exists. Skipping MSGFPlusAdapter...'.format(idXML))

        if os.path.isfile(idXML) and os.path.getsize(idXML) > 0:
            if not os.path.isfile(pi_idXML) or os.path.getmtime(idXML) > os.path.getmtime(pi_idXML):
                peptideindex_cmd = "PeptideIndexer -in {} -fasta {} -out {} -allow_unmatched -enzyme:specificity 'semi'".format(idXML, fasta, pi_idXML)
                run_cmd(peptideindex_cmd,'PeptideIndexer {} accomplished!'.format(idXML))
            else:
                print('File {} exists and seems up-to-date. Skipping PeptideIndexer ...'.format(pi_idXML))
        else:
            print('Error no idXML file generated!')
            sys.exit(1)

        status_message('###############################################')
        shutil.rmtree('{}'.format(os.path.join(tmp_dir, str(rand))))

    status_message('MSGFPlus and PeptideIndexer accomplished!')
    status_message('Merging and filtering begins!')

    if not os.path.isdir(os.path.join(output_dir,'results')):
        os.mkdir('{}'.format(os.path.join(output_dir,'results')))
    results_dir = os.path.join(output_dir, 'results')
    merged_file = '{}/merged.idXML'.format(results_dir)
    fdr_file = '{}/fdr.idXML'.format(results_dir)
    filtered_file = '{}/filtered.idXML'.format(results_dir)
    file_info = '{}/info.txt'.format(results_dir)
    filtered_csv = '{}/filtered.csv'.format(results_dir)

    pi_idXMLs = []
    for j in glob.iglob('{}/*_pi.idXML'.format(output_dir)):
        pi_idXMLs.append(j)

    id_merge = 'IDMerger -in {} -out {}'.format(' '.join(pi_idXMLs), merged_file)
    fdr_cal = 'FalseDiscoveryRate -in {} -out {}'.format(merged_file, fdr_file)
    id_filter = 'IDFilter -in {} -out {} -score:pep 0.05'.format(fdr_file, filtered_file)
    file_info_gen = 'FileInfo -in {} | tee {}'.format(filtered_file, file_info)
    to_csv = 'TextExporter -in {} -out {}'.format(filtered_file, filtered_csv)

    if run_cmd(id_merge, '{} generated!'.format(merged_file)):
        if run_cmd(fdr_cal, '{} generated!'.format(fdr_file)):
            if run_cmd(id_filter, '{} generated!'.format(filtered_file)):
                if run_cmd(file_info_gen, '{} generated!'.format(file_info)):
                    if run_cmd(to_csv, '{} generated!'.format(filtered_csv)):
                        status_message('IDMerger, FalseDiscoveryRate, IDFilter and TextExporter accomplished!')
                        end_execute_timestamp = datetime.datetime.now()
                        elapsed_time = ( end_execute_timestamp - start_execute_timestamp ).total_seconds()
                        status_message('Processing {} completed, consumed {} seconds'.format(sample_name, elapsed_time))
                        status_message('<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<')
