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
    #family = base.split("_")[0]
    #baltimore = get_baltimore_classifiaction(family)
    #print('get_mutation_joint_df_codon')
    #map_file = fastml_utilities.get_mutation_joint_df_codon(base, infoDictionary={"base":base, "family":family, "baltimore":baltimore}, overwrite=False)
    #fastml_utilities.summmerize_mutations(map_file, overwrite=True)
    fastml_utilities.get_mutation_context_apobec(base)





if __name__ == "__main__":
    main()

