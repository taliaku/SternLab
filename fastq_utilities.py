#! /usr/local/python_anaconda/bin/python3.4

import re
import os
from file_utilities import check_dirname, make_dir

def filter_fastq_by_read_id(fastq_path, filtered_fastq_path, read_id_list):
    '''
    Filters fastq file, keeps only reads in the given list of read ids.
    :param fastq_path: input fastq path.
    :param filtered_dastq_path: path to create filtered fastq in.
    :param read_id_list: list of read ids to keep.
    '''
    new_fastq = ''
    with open(fastq_path) as f:
        old_fastq_str = f.read()
    for read_id in read_id_list:
        p = re.compile('(@' + read_id + '\n.*?\n)@M', re.DOTALL)
        read_data = p.findall(old_fastq_str)[0]
        new_fastq += read_data
    with open(filtered_fastq_path, 'w') as f:
        f.write(new_fastq)
    return

def filter_out_fastq_by_read_id(fastq_path, filtered_fastq_path, read_id_list):
    '''
    Filters fastq file, keeps only reads *not* in the given list of read ids.
    :param fastq_path: input fastq path.
    :param filtered_dastq_path: path to create filtered fastq in.
    :param read_id_list: list of read ids to filter out.
    '''
    with open(fastq_path) as f:
        fastq_str = f.read()
    for read_id in read_id_list:
        p = re.compile('(@' + read_id + '\n.*?\n)@M', re.DOTALL)
        read_data = p.findall(fastq_str)[0]
        fastq_str = fastq_str.replace(read_data, '')
    with open(filtered_fastq_path, 'w') as f:
        f.write(fastq_str)
    return

def unite_nextseq_lanes(input_dir, output_dir):
    '''
    Gets a directory with fastq files from nextseq, and unites the different lanes
    into one file for each sample.
    :param input_dir: directory with fastq files from the run.
    :param output_dir: directory to create and write new fastq files to.
    '''
    input_dir = check_dirname(input_dir)
    make_dir(check_dirname(output_dir, Truedir=False))
    files = [input_dir + f for f in os.listdir(input_dir)]
    samples = list(set([f.split('/')[-1].split('_L001')[0] for f in files if '_L001' in f]))
    for sample in samples:
        for r in ['R1', 'R2']:
            os.system('cat %s%s_L00[1234]_%s_001.fastq > %s%s_%s.fastq' % (input_dir, sample, r, output_dir, sample, r))
    return 

def merge_every_sample(input_dir, output_dir):
    '''
    Gets a directory with fastq files, and submits to que merge on R1 and R2 for every sample.
    Merged files are submitted to output_dir.
    :param input_dir: directory with R1 and R2 fastq files 
    :param output_dir: directory to create and write merged fastq files to.
    '''
    input_dir = check_dirname(input_dir)
    make_dir(check_dirname(output_dir, Truedir=False))
    pass
    