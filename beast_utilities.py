#! /usr/local/python_anaconda/bin/python3.4

import os
from os import path
from file_utilities import set_filenames_for_pbs_runs, check_filename, check_dirname
import re

def configure_log_path(xmlpath, logpath, out=None):
    """
    self configuration of the log file in the xml
    :param xmlpath: a path to an xml file created by BEAUTi
    :param logpath: the path were the log file will be saved 
    :param out: a path to save the new xml file
    return: saves a new xml input file with a configured log
    """
    xmlpath = check_filename(xmlpath)

    if out == None:
        out=xmlpath

    with open(xmlpath, 'r') as base_file, open(out, 'w') as out_file:
        text = base_file.read()
        text_2_write = re.sub(r'"([\w+/+]*.log)"',logpath, text)
        cnt = out_file.write(text_2_write)
    return cnt
