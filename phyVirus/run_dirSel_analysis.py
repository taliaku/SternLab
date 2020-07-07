#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import *



def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-b", "--basename", dest="basename", help="basename")
    parser.add_option("-o", "--overwrite", action="store_true", default=False, dest="overwrite", help="overwrite files")
    parser.add_option("-m", "--midpoint", action="store_true", default=False, dest="midpoint", help="is this midpoint rooting dirSel?")





    (options, args) = parser.parse_args()
    basename = options.basename
    if basename[-1] != ".":
        basename += "."
    overwrite = options.overwrite
    midpoint = options.midpoint
    dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"

    dirSel_analysis

    if midpoint:
        mutation_map_file = dir + "phyVirus_analysis_midpoint/%sdirSel_midpoint.results.mutation.map" % basename
        dirSel_analysis = dir + "dirSel_analysis_midpoint/"
    else:
        mutation_map_file = dir + "dirSel/%sdirSel.results.mutation.map" % basename
        dirSel_analysis = dir + "dirSel_analysis/"
    mutation_map_file = check_filename(mutation_map_file)

    aln_output_file = dirSel_analysis" /%sdirSel.aln" % basename
    aln_output_file = check_filename(aln_output_file, Truefile=False)

    if not os.path.isfile(aln_output_file):
        alignment_from_mutation_map_info(mutation_map_file, aln_output_file)



    several_context_output_file = dirSel_analysis + "%sdirSel.context" % basename
    several_context_output_file = check_filename(several_context_output_file, Truefile=False)
    run_several_context_analysis_on_file(mutation_map_file, output=several_context_output_file, aln_type="original", overwrite=overwrite)
    print(several_context_output_file)

    several_context_aln_mut_map_output_file = dirSel_analysis + "%sdirSel.context_aln_mut_map" % basename
    several_context_aln_mut_map_output_file = check_filename(several_context_aln_mut_map_output_file, Truefile=False)
    run_several_context_analysis_on_file(mutation_map_file, output=several_context_aln_mut_map_output_file,
                                         aln_type="from_mut_map", aln_file=aln_output_file, overwrite=overwrite)

    analyze_dirSel_mutation_map_output_file = dirSel_analysis + "%sdirSel.ratios" % basename
    analyze_dirSel_mutation_map_output_file = check_filename(analyze_dirSel_mutation_map_output_file, Truefile=False)
    analyze_dirSel_mutation_map(mutation_map_file, ratios_output=analyze_dirSel_mutation_map_output_file, overwrite=overwrite)

    analyze_dirSel_mutation_map_to_nucs_output_file = dirSel_analysis + "%sdirSel.ratios_to_nucs" % basename
    analyze_dirSel_mutation_map_to_nucs_output_file = check_filename(analyze_dirSel_mutation_map_to_nucs_output_file, Truefile=False)
    analyze_dirSel_mutation_map_to_nucs(mutation_map_file, ratios_output=analyze_dirSel_mutation_map_to_nucs_output_file, overwrite=overwrite)

    analyze_dirSel_mutation_map_to_nucs_several_cutoffs_output_file =  dirSel_analysis + "%sdirSel.ratios_to_nucs_cutoffs" % basename
    analyze_dirSel_mutation_map_to_nucs_several_cutoffs_output_file = check_filename(analyze_dirSel_mutation_map_to_nucs_several_cutoffs_output_file, Truefile=False)
    analyze_dirSel_mutation_map_to_nucs_several_cutoffs(mutation_map_file, ratios_output=analyze_dirSel_mutation_map_to_nucs_several_cutoffs_output_file, overwrite=overwrite)

    analyze_dirSel_mutation_map_to_nucs_updated_output_file = dirSel_analysis + "%sdirSel.ratios_to_nucs_updated" % basename
    analyze_dirSel_mutation_map_to_nucs_updated_output_file = check_filename(analyze_dirSel_mutation_map_to_nucs_updated_output_file, Truefile=False)
    analyze_dirSel_mutation_map_to_nucs_updated(mutation_map_file, ratios_output=analyze_dirSel_mutation_map_to_nucs_updated_output_file, overwrite=overwrite)




if __name__ == "__main__":
    main()

