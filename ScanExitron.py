#!/usr/bin/env python3
#-*- coding: utf-8 -*-
#===============================================================================
#===============================================================================
import sys
import re
import os
import argparse
import numpy as np
import glob
from pyfaidx import Fasta
import subprocess
import random
import shutil
import pysam
from collections import OrderedDict
from io import BytesIO ## for Python 3
try:
    import configparser
except ImportError:
    import ConfigParser as configparser
    
    
def get_value_by_key(src, delimiter='='):
    out_dict = {}
    src = src.strip().strip(';')
    tmpList = re.split(r';\s{0,2}', src)■
    for i in tmpList:
        if i:
            k = i.replace('"', '')
            m, n = k.split(delimiter)
            out_dict[m] = n
    return out_dict

   

def status_message(msg):
    print(msg)
    sys.stdout.flush()


def config_getter():
    this_dir = os.path.dirname(os.path.realpath(__file__))
    config_default = os.path.join(this_dir, 'config.ini')
    config = configparser.SafeConfigParser(os.environ)
    config.read(config_default)
    hg38_ref   = config.get("fasta", "hg38")
    hg19_ref   = config.get("fasta", "hg19")
    hg38_anno  = config.get("annotation", "hg38")
    hg19_anno  = config.get("annotation", "hg19")
    hg38_cds = config.get("cds", "hg38")
    hg19_cds = config.get("cds", "hg19")
    return {'hg38_ref': hg38_ref, 'hg19_ref': hg19_ref,'hg38_anno':hg38_anno, 'hg19_anno':hg19_anno, 'hg38_cds':hg38_cds, 'hg19_cds':hg19_cds}

chrms = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5',
    'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11','chr12', 'chr13', 'chr14', 'chr15',
    'chr16','chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY', 'chrM'}

non_mito_chrms = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5',
    'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
    'chr11','chr12', 'chr13', 'chr14', 'chr15',
    'chr16','chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX', 'chrY'}

chrms_dict = {'1':'chr1', '2':'chr2', '3':'chr3', '4':'chr4', '5':'chr5',
        '6':'chr6', '7':'chr7', '8':'chr8', '9':'chr9', '10':'chr10',
        '11':'chr11','12':'chr12', '13':'chr13', '14':'chr14', '15':'chr15',
        '16':'chr16','17':'chr17', '18':'chr18', '19':'chr19', '20':'chr20',
        '21':'chr21', '22':'chr22', 'X':'chrX', 'Y':'chrY', 'MT':'chrM'}

reverse_chrms_dict = dict((chrms_dict[i], i) for i in chrms_dict)

def BED_handler(inbed):
    ''' keep only canonical chromosomes and convert b37/b38 to hg19/38
    '''
    x_ = int(random.random()*10000)
    tmp_file = 'tmp.{}{}.txt'.format(os.getpid(), x_)
    tmp = open(tmp_file, 'w')
    with open(inbed, 'r') as f:
        for line in f:
            l = line.rstrip('\n').split('\t')
            if l[0] in chrms:
                tmp.write('\t'.join(l)+'\n')
            elif l[0] in chrms_dict:
                l[0] = chrms_dict[l[0]]
                tmp.write('\t'.join(l)+'\n')
    tmp.close()
    shutil.move(tmp_file, inbed)
    return os.path.abspath(inbed)

def junction_caller(bam_file, ref='hg38', out_name=None, config=config_getter()):
    '''
    Call junctions using regtools
    output: out_name.janno
    '''
    if ref == 'hg19':
        fasta == config['hg19_ref']
        gtf = config['hg19_anno']
    elif ref == 'hg38':
        fasta == config['hg38_ref']
        gtf = config['hg38_anno']

    prefix = os.path.splitext(os.path.basename(bam_file))[0] 
    cmd = 'regtools junctions extract -i 5 -I 10000000 {} -o {}.bed'.format(bam_file, prefix) 
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error! {}'.format(err)
    else:
        error_msg = ''
    if not error_msg:
        junction_flag = True
        print('Calling junctions finished!')

    if junction_flag:
        bed = BED_handler('{}.bed'.format(prefix))

    if not out_name:
        out_name = prefix
    cmd = 'regtools junctions annotate {0} {1} {2} -o {3}.janno'.format(bed, fasta, gtf, out_name)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error! {}'.format(err.output)
    else:
        error_msg = ''
    if not error_msg:
        print('{}.janno generated!'.format(out_name))
        os.remove('{}'.format(bed))
        return '{}.janno'.format(out_name)


def junction_overlap_CDS_to_position_BED(janno, ref='hg38', config=config_getter()):
    '''
        intersect junctions with annotated CDS to search exitrons 
    '''
    if ref == 'hg19':
        cds = config['hg19_cds']
    elif ref == 'hg38':
        cds = config['hg38_cds']

    genome_seq = seq_dict(ref=ref)

    print('Reading {}'.format(janno))

    # write all the novel junctions with canonical splicing sites to file (junction.bed)
    junction_bed = '{}.junction.bed'.format(os.getpid())
    total_junctions = 0
    out = open(junction_bed, 'w')
    with open(janno) as f:
        f.readline()
        for line in f:
            l = line.rstrip().split('\t')
            total_junctions += int(l[4])
            chrm = l[0]
            start = int(l[1])
            end = int(l[2])
            stats = l[10]
            strand = l[5]
            spliced_site = l[6].upper()
            #if stats == 'N' and strand != '?' and spliced_site in {'GT-AG','GC-AG','AT-AC'}:
            if stats == 'N' and strand != '?':
                if strand == '+':
                    left_site = genome_seq[chrm][start:start+2].seq
                    right_site = genome_seq[chrm][end-3:end-1].seq
                elif strand== '-':
                    left_site = genome_seq[chrm][end-3:end-1].reverse.complement.seq
                    right_site = genome_seq[chrm][start:start+2].reverse.complement.seq
                l[6] = '{}-{}'.format(left_site, right_site)
                if l[6] in {'GT-AG','GC-AG','AT-AC'}:
                    out.write('{}\n'.format('\t'.join(l[:7])))
    out.close()

    overlap_file = '{}.overlap.bed'.format(os.getpid())             
    cmd = 'bedtools intersect -s -wo -a {} -b {} > {}'.format(junction_bed, cds, overlap_file)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'error!'
    else:
        error_msg = ''
    if not error_msg:
        pass
    # no overlap in CDS and junctions file
    if os.path.isfile(overlap_file) and os.path.getsize(overlap_file) == 0:
        os.remove(overlap_file)            
        os.remove(junction_bed)            
        print('No overlaps found in {} and gencode CDS'.format(janno))
        return False 
    # overlaps found in CDS and junctions file
    elif os.path.isfile(overlap_file) and os.path.getsize(overlap_file) > 0:
        tmp_dict =  OrderedDict()
        with open(overlap_file) as f:
            for line in f:
                l = line.rstrip().split('\t')
                chrm = l[0]
                length = int(l[-1])
                junc_start = int(l[1])
                junc_end = int(l[2])
                junc_id = l[3]
                junc_read_no = l[4]
                strand = l[5]
                splice_site = l[6]
                ref_start = int(l[8])
                ref_end = int(l[9])
                attr = get_value_by_key(l[16], delimiter=' ')
                #exon_number = attr['exon_number']
                gene_name = attr['gene_name']
                gene_id = attr['gene_id']
                pos_key = '{}:{}-{}'.format(chrm, junc_start, junc_end)
                if length == junc_end - junc_start  and junc_start > ref_start and junc_end < ref_end and chrm in non_mito_chrms and int(junc_read_no)>=3:
                    if not pos_key in tmp_dict:
                        info = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(chrm, junc_start, junc_end, junc_id, junc_read_no, strand, gene_name, length-1, splice_site, gene_id, total_junctions)
                        tmp_dict[pos_key] = info
        os.remove(overlap_file)
        os.remove(junction_bed)
        # exitrons found
        if len(tmp_dict) > 0:
            #x_ = int(random.random()*10000)
            position_bed_file = os.path.splitext(os.path.basename(janno))[0] + '.position.bed'
            out = open(position_bed_file, 'w')

            src_exitron_file = os.path.splitext(os.path.basename(janno))[0] + '.src'
            output = open(src_exitron_file, 'w')

            position_set = set([])
            for i in tmp_dict:
                chrm, junc_start, junc_end, junc_id, junc_read_no, strand, gene_name, junc_len, splice_site, gene_id, total_junctions = tmp_dict[i].split('\t')

                #chrm = reverse_chrms_dict[chrm]
                output.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(chrm, junc_start, junc_end, junc_id, junc_read_no, strand, gene_name, junc_len, splice_site, gene_id, total_junctions))
                start = int(junc_start)
                end = int(junc_end)
                #ao = int(junc_read_no)
                middle_point = int(np.median([start, end]))

                if not '{}\t{}'.format(chrm, start) in position_set:
                    out.write('{}\t{}\t{}\n'.format(chrm, start-1, start))
                    position_set.add('{}\t{}'.format(chrm, start))

                if not '{}\t{}'.format(chrm, end) in position_set:
                    out.write('{}\t{}\t{}\n'.format(chrm, end-1, end))
                    position_set.add('{}\t{}'.format(chrm, end))

                if not '{}\t{}'.format(chrm, middle_point) in position_set:
                    out.write('{}\t{}\t{}\n'.format(chrm, middle_point-1, middle_point))
                    position_set.add('{}\t{}'.format(chrm, middle_point))

            output.close()
            out.close()
        else:
            print('No exitron found in {}'.format(janno))
            return False
    return src_exitron_file, position_bed_file


def percent_spliced_out(bam_file, src_exitron_file, position_bed_file, mapq):
     
    print('Reading BAM file: {}'.format(bam_file))
    depth_dict = {}
    cmd = 'samtools bedcov {0} {1} -Q {2}'.format(position_bed_file, bam_file, mapq)
    try:
        result = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT,)
    except subprocess.CalledProcessError as err:
        error_msg = 'Error happend!: {}\n{}'.format(err, err.output)
    else:
        error_msg = ''
    if not error_msg:
        result_file = BytesIO(result)
        result_string = result_file.getvalue().decode("utf-8")
        for line in result_string.split('\n'):
            if line:
                chrm, _, pos, depth = line.rstrip().split()
                depth_dict['{}\t{}'.format(chrm, pos)] = int(depth)
    else:
        status_message(error_msg)
    outfile = os.path.basename(src_exitron_file).split('.')[0] + '.exitron'
    
    out = open(outfile, 'w')
    
    with open(src_exitron_file) as f:
        for line in f:
            l = line.rstrip('\n').split('\t')
            chrm = l[0]
            start = int(l[1])
            end = int(l[2])
            ao = int(l[4])
            strand = l[5]
            middle_point = int(np.median([start, end]))

            if strand == '+':
                five_prime_reads = depth_dict['{}\t{}'.format(chrm, start)] - ao
                three_prime_reads = depth_dict['{}\t{}'.format(chrm, end)] - ao
                middle_reads = depth_dict['{}\t{}'.format(chrm, middle_point)] - ao
            elif strand == '-':
                five_prime_reads = depth_dict['{}\t{}'.format(chrm, end)] - ao
                three_prime_reads = depth_dict['{}\t{}'.format(chrm, start)] - ao
                middle_reads = depth_dict['{}\t{}'.format(chrm, middle_point)] - ao
            ave_dp = (five_prime_reads+three_prime_reads+middle_reads)/3.0
            if five_prime_reads < 0 or three_prime_reads < 0:
                continue
            try:
                pso = float(ao) / (ave_dp + ao)
            except ZeroDivisionError:
                print('Error in {} {} {}'.format(chrm, junc_start, junc_end))
                pso = 0
            psi = 1.0 - pso
            dp = int(ao/pso)
            out.write('{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(l[:-1]), pso, psi, dp, l[-1]))
    os.remove(src_exitron_file)
    os.remove(position_bed_file)
    out.close()
    print('Finished reading BAM file: {}'.format(bam_file))
    return outfile

def external_tool_checking():
    """checking dependencies are installed"""
    software = ['regtools', 'bedtools', 'samtools']
    cmd = "which"
    for each in software:
        try:
            path = subprocess.check_output([cmd, each], stderr=subprocess.STDOUT)
            path = str(path, 'utf-8')
        except subprocess.CalledProcessError:
            print("Checking for '" + each + "': ERROR - could not find '" + each + "'", file=sys.stderr)
            print("Exiting.", file=sys.stderr)
            sys.exit(0)
        print("Checking for '" + each + "': found " + path)


def seq_dict(ref='hg38', config=config_getter()):
    if ref =='hg19':
        fasta = config['hg19_ref']
    elif ref == 'hg38':
        fasta = config['hg38_ref']
    genome_dict = Fasta(fasta, sequence_always_upper=True)
    return genome_dict



def parse_args():
    parser = argparse.ArgumentParser(description = "%(prog)s -i input_rna_seq_bam_file -r [hg38/hg19] -m mapping_quality", epilog="ScanExitron: detecting exitron splicing using RNA-Seq data")
    parser.add_argument('-i', '--input', action='store', dest='input', help="Input BAM/CRAM file along with BAI/CRAI file", required=True)
    parser.add_argument('-m', '--mapq', action='store', dest='mapq', type=int, help="mapq (default: %(default)s)", default=0)
    parser.add_argument('-r', '--ref', action='store', dest='ref', help="reference (default: %(default)s)", choices=['hg19', 'hg38'], default='hg38')
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()
    return args

def main():
    external_tool_checking()
    args = parse_args()

    janno_file = junction_caller(bam_file=args.input, ref=args.ref)
    src_exitron_file, position_bed_file = junction_overlap_CDS_to_position_BED(janno_file, ref=args.ref)

    percent_spliced_out(bam_file=args.input, src_exitron_file=src_exitron_file, position_bed_file=position_bed_file, mapq=args.mapq)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
