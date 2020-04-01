#! /usr/local/python_anaconda/bin/python3.4

import os
from os import path
from file_utilities import set_filenames_for_pbs_runs, check_filename, check_dirname
import re

def configure_log_path(xml_in, logpath, xml_out=None):
    """
    self configuration of the log file in the xml
    :param xml_in: a path to an xml file created by BEAUTi
    :param logpath: the path were the log file will be saved 
    :param xml_out: a path to save the new xml file
    return: saves a new xml input file with a configured log
    """
    xml_in = check_filename(xml_in)

    if xml_out == None:
        xml_out=xml_in

    with open(xml_in, 'r') as base_file:
        text = base_file.read()
        text_2_write = re.sub(r'([\w+/+]*.log)',logpath, text)
    with open(xml_out, 'w') as out_file:
        cnt = out_file.write(text_2_write)
    return cnt


def replace_beast_configs(xml_in, template, xml_out=None):
    xml_in = check_filename(xml_in)

    if xml_out == None:
        xml_out = xml_in

    with open(template, "r") as t:
        template_configs = t.read()
    try:
        with open(xml_in, "r") as fp:
            content = fp.read()
            loc = re.search("	<patterns id.*", content)
            if loc:
                content = content[:loc.start()]
        with open(xml_out, "w") as fp:
            print(xml_out)
            fp.write(content + template_configs)

    except Exception as e:
        print(e, "error", xml_in, xml_out)