#! /usr/local/python_anaconda/bin/python3.4

import re
import os
from file_utilities import check_dirname, make_dir
from pbs_runners import merge_runner, pipeline_runner

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
    files = [input_dir + '/' + f for f in os.listdir(input_dir)]
    samples = list(set([f.split('/')[-1].split('_L001')[0] for f in files if '_L001' in f]))
    for sample in samples:
        for r in ['R1', 'R2']:
            os.system('cat %s%s_L00[1234]_%s_001.fastq > %s%s_%s.fastq' % (input_dir, sample, r, output_dir, sample, r))
    return 

def merge_every_sample(input_dir, output_dir, individual_directories=True):
    '''
    Gets a directory with fastq files, and submits to que merge on R1 and R2 for every sample.
    Merged files are submitted to output_dir.
    :param input_dir: directory with R1 and R2 fastq files 
    :param output_dir: directory to create and write merged fastq files to.
    :param individual_directories: True or False, should every merge file be in an
            individual directory inside input_dir, or files directly in input_dir.
    '''
    input_dir = check_dirname(input_dir)
    make_dir(check_dirname(output_dir, Truedir=False))
    files = [input_dir +  '/' + f for f in os.listdir(input_dir)]
    job_ids = []
    for f in files:
        if 'R1.fastq' in f.split('/')[-1]:
            if individual_directories == False:
                job_id = merge_runner(f, 
                                      f.replace('R1.fastq', 'R2.fastq'), 
                                      output_dir + f.split('/')[-1].replace('_R1', ''))
            else:
                individual_dir = output_dir + '/' + f.split('/')[-1].split('_R1')[0] + '/'
                make_dir(check_dirname(individual_dir, Truedir=False))
                job_id = merge_runner(f, 
                                      f.replace('R1.fastq', 'R2.fastq'), 
                                      individual_dir + f.split('/')[-1].replace('_R1', ''))
            job_ids.append(job_id)
    return job_ids

def pipeline_every_sample(input_dir, output_dir, ref_file, NGS_or_Cirseq, TYPE_OF_INPUT_FILE=None, start=None, end=None, gaps=None,
                    qscore=None, blast=None, rep=None, t=None, alias="pipeline"):
    '''
    Run the pipeline with the same parameters for many samples.
    Gets a directory with directories contatining fastq (or gz) files, 
    treats every inner directory as a sample and runs the pipeline on it. 
    Creates inner directory for very sample with the results inside output_dir.
    Other then input and output directories, gets all the varaibles that 
    pipeline_runner from pbs_runners gets.
    '''
    input_dir = check_dirname(input_dir)
    make_dir(check_dirname(output_dir, Truedir=False))
    dirs = [input_dir +  '/' + f for f in os.listdir(input_dir)]
    job_ids = []
    for d in dirs:
        out_d = output_dir + '/' + d.split('/')[-1]
        make_dir(check_dirname(out_d, Truedir=False))
        job_id = pipeline_runner(d, out_d, ref_file, NGS_or_Cirseq, 
                                 TYPE_OF_INPUT_FILE=TYPE_OF_INPUT_FILE, start=start, end=end, 
                                 gaps=gaps, qscore=qscore, blast=blast, rep=rep, 
                                 t=t, alias=alias)
        job_ids.append(job_id)
    return job_ids

    