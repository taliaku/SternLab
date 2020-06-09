#! /usr/local/python_anaconda/bin/python3.4

import os
from os import path
from file_utilities import set_filenames_for_pbs_runs, check_filename, check_dirname
import re

def get_sample_info_from_bowtie_command_line(command_line):
        input_file = command_line.split("-1")[1].split("-2")[0].strip().split("/")[-1]
        patient = input_file.split("_")[0]
        gen_mat = input_file.split("_")[1]
        lane = int(input_file.split("_")[3].split("L00")[1].split(".")[0])
        return(patient, gen_mat, lane)