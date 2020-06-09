#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')

import os
from optparse import OptionParser
from file_utilities import check_filename, check_dirname
from PAML_utilities import write_ctl_file
import pandas as pd
import glob
from selecton_utilities import extract_selecton_final_params_single_file
from pbs_runners import script_runner
import os
from general_utilities import frange
from dirSel_utilities import run_several_context_analysis_on_file, analyze_context_dirSel, alignment_from_mutation_map_info, analyze_dirSel_mutation_map_to_nucs_several_cutoffs



def main():
    
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="file", help="file")


    (options, args) = parser.parse_args()
    file = options.file
    file = check_filename(file)
    basename = file.split("/")[-1].split("dirSel.results.mutation.map")[0]
    #analyze_context_dirSel(file, letter1="G", letter2="A", before_or_after="after", overwrite=True)
    #run_several_context_analysis_on_file(file, overwrite=True, mutations=["GA", "AG", "CT", "TC"], aln_type="from_mut_map")
    #alignment_from_mutation_map_info(file)
    output =  "/sternadi/home/volume3/taliakustin/phyVirus_analysis/" + "dirSel_analysis/%sdirSel.ratios_to_nucs_cutoffs" % basename


    analyze_dirSel_mutation_map_to_nucs_several_cutoffs(file, ratios_output=output, overwrite=True)


if __name__ == "__main__":
    main()

