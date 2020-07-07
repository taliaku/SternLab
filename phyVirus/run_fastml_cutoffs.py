#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from phyVirus_utilities import get_basename, get_baltimore_classifiaction
import fastml_utilities


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="fasta", help="fasta file")

    (options, args) = parser.parse_args()

    fasta = options.fasta
    fasta = check_filename(fasta)

    base = get_basename(fasta).split(".")[0]
    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)

    joint_map_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/fastml_midpoint_tree/%s*joint*map" % base)[0]

    branch_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/branch_length_info/%s.branch_lengths_info.csv" % base)[0]
    print(joint_map_file, branch_file)
    fastml_utilities.different_cutoffs(joint_map_file, branch_file)






if __name__ == "__main__":
    main()

