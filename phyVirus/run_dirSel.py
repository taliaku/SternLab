#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file

DIRSEL_PATH = "/sternadi/home/volume1/taliakustin/software/phylogenyCode/programs/directionalSelection/directionalSelection"

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-a", "--aln", dest="aln", help="alignment file")
    parser.add_option("-t", "--tree", dest="tree", help="tree file")
    parser.add_option("-m", "--model", dest="model", help="model - fixed_beta or alternative")
    parser.add_option("-g", "--gtr", dest="gtr", help="gtr output")
    parser.add_option("-d", "--dirSel_path", dest="dirSel_path", help="dirSel path",
                      default="/sternadi/home/volume1/taliakustin/software/phylogenyCode/programs/directionalSelection/directionalSelection")

    (options, args) = parser.parse_args()

    aln = options.aln
    aln = check_filename(aln)
    tree = options.tree
    tree = check_filename(tree)
    gtr = options.gtr
    gtr = check_filename(gtr)
    model = options.model
    dirSel_path = options.dirSel_path

    if not os.path.isfile(dirSel_path):
        raise("dirSel path is not correct")
    print(model)
    #if model not in ["fixed_beta", "alternative", "extInt", "dirSel_extIntNull", "dirSel_A", "dirSel_A_null"]:
     #   raise("Model must be fixed_beta or alternative or extInt or dirSel_extIntNull or dirSel_A or dirSel_A_null")
    init_beta = 0
    fixed_tau = 1
    fixed_kappa = 1
    bblOpt = 0
    if model == "toG_alt":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toG_alt/"
        fixed_beta = 0
        fixed_s = 0
        fixed_probS = 0
        init_probS = 0.01
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toG_with_posterior/programs/directionalSelection/directionalSelection"
    elif model == "toU_alt":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toU_alt/"
        fixed_beta = 0
        fixed_s = 0
        fixed_probS = 0
        init_probS = 0.01
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toU_with_posterior/programs/directionalSelection/directionalSelection"

    elif model == "toG_null":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toG_null"
        fixed_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toG_with_posterior/programs/directionalSelection/directionalSelection"
    elif model == "toU_null":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toU_null"
        fixed_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toU_with_posterior/programs/directionalSelection/directionalSelection"

    """
    if model == "beta_null":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_null_fixed_tau"
        fixed_beta = 1
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root/programs/directionalSelection/directionalSelection"
    elif model == "beta_alt":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_fixed_tau"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root/programs/directionalSelection/directionalSelection"
    elif model =="beta_0.05":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_0.05_fixed_tau"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_cutoff_0.05/programs/directionalSelection/directionalSelection"
    elif model =="beta_0.01":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_0.01_fixed_tau"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_cutoff_0.01/programs/directionalSelection/directionalSelection"
    elif model =="beta_0.005":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_0.005_fixed_tau"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_cutoff_0.005/programs/directionalSelection/directionalSelection"
    elif model =="beta_0.001":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_0.001_fixed_tau"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_cutoff_0.001/programs/directionalSelection/directionalSelection"
    elif model =="beta_0.1":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_0.1_fixed_tau"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_cutoff_0.1/programs/directionalSelection/directionalSelection"
    elif model=="beta_some_categories":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_some_categories"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_beta_some_categories/programs/directionalSelection/directionalSelection"
    elif model == "beta_null_some_categories":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_null_some_categories"
        fixed_beta = 1
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 1

    elif model == "beta_null_optimized":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_null_optimized"
        fixed_beta = 1
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 0
        fixed_kappa = 0
        bblOpt = 0
        optimizeLineSearch = 0


    elif model == "beta_alt_optimized":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_optimized"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 0
        fixed_kappa = 0
        bblOpt = 0
        optimizeLineSearch = 0
    elif model == "beta_alt_optimized_some":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/beta_alt_optimized_some_categories"
        fixed_beta = 0
        init_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 0
        fixed_kappa = 0
        bblOpt = 0
        optimizeLineSearch = 0

    """

    """
    if model == "toA_null":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toA_null_pos"
        fixed_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_with_posterior/programs/directionalSelection/directionalSelection"
    elif model == "toA_alt":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toA_alt_pos/"
        fixed_beta = 0
        fixed_s = 0
        fixed_probS = 0
        init_probS = 0.01
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_with_posterior/programs/directionalSelection/directionalSelection"
    elif model == "dirSel_withS":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_withS_noRoot/"
        fixed_beta = 0
        fixed_s = 0
        fixed_probS = 0
        init_probS = 0.01
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root/programs/directionalSelection/directionalSelection"
    elif model == "dirSel_null":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_null_noRoot/"
        fixed_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root/programs/directionalSelection/directionalSelection"
    elif model=="Hepe":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_tmp/"
        fixed_beta = 0
        fixed_s = 0
        fixed_probS = 1
        init_probS = 0.1
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root/programs/directionalSelection/directionalSelection"


    elif model == "toC_null":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toC_null"
        fixed_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toC_with_posterior/programs/directionalSelection/directionalSelection"
    elif model == "toC_alt":
        inter = "dirSel"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toC_alt/"
        fixed_beta = 0
        fixed_s = 0
        fixed_probS = 0
        init_probS = 0.01
        path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toC_with_posterior/programs/directionalSelection/directionalSelection"
    """
    """
        if model == "hky":
            inter = "dirSel"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/hky_midpoint"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0
            path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA/programs/directionalSelection/directionalSelection"
        elif model == "toA":
            inter = "dirSel"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toA_midpoint/"
            fixed_beta = 0
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
            path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA/programs/directionalSelection/directionalSelection"
        elif model == "toA_internal":
            inter = "dirSel"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/toA_internal_midpoint/"
            fixed_beta = 0
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
            path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_internal/programs/directionalSelection/directionalSelection"
        elif model == "hky_internal":
            inter = "dirSel"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/hky_midpoint_internal/"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0
            path = "/sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_internal/programs/directionalSelection/directionalSelection"
    
    """
    """
        if model == "midpoint":
            inter = "dirSel_midpoint"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_midpoint//"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0.01
        elif model == "midpoint_null":
            inter = "dirSel_midpoint_null"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_midpoint_null/"
            fixed_beta = 1
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0.01
        elif model == "toA_internal":
            inter = "dirSel_selection_toA_internal"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_selection_toA_internal/"
            fixed_beta = 0
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
        elif model ==  "toA_internal_null":
            inter = "dirSel_selection_toA_internal_null"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_selection_toA_internal_null/"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0.01
        elif model == "toA_external":
            inter = "dirSel_selection_toA_external"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_selection_toA_external/"
            fixed_beta = 0
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
        elif model == "toA_external_null":
            inter = "dirSel_selection_toA_external_null"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_selection_toA_external_null/"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0.01
        """
    """
        if model == "alternative":
            inter = "dirSel"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel/"
            fixed_beta = 0
            fixed_s = 1
            fixed_ProbS = 1
            init_probS = 0.01
        elif model == "fixed_beta":
            inter = "dirSel_beta0"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_beta0/"
            fixed_beta = 1
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0.01
        elif model == "extInt":
            inter = "dirSel_extInt"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_extInt/"
            fixed_beta = 0
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
        elif model == "dirSel_extIntNull":
            inter = "dirSel_extIntNull"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_extInt_null/"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0
        elif model == "dirSel_A_internal":
            inter = "dirSel_A_internal"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_A_internal/"
            fixed_beta = 0
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
        elif model == "dirSel_A_internal_null":
            inter = "dirSel_A_internal_null"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_A_internal_null/"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0
        elif model == "dirSel_A_external":
            inter = "dirSel_A_external"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_A_external/"
            fixed_beta = 0
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
        elif model == "dirSel_A_external_null":
            inter = "dirSel_A_external_null"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_A_external/"
            fixed_beta = 0
            fixed_s = 1
            fixed_probS = 1
            init_probS = 0
        elif model == "dirSel_A_internal_fixed_beta":
            inter = "dirSel_A_internal_fixed_beta"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_A_internal_fixed_beta/"
            fixed_beta = 1
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
        elif model == "dirSel_A_external_fixed_beta":
            inter = "dirSel_A_external_fixed_beta"
            output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/dirSel_A_external_fixed_beta/"
            fixed_beta = 1
            fixed_s = 0
            fixed_probS = 0
            init_probS = 0.01
        """

    basename = aln.split("/")[-1].split("aln.best.fas")[0]

    param_file = output_dir + "/%s%s.params" % (basename, inter)

    output_res = output_dir + "/%s%s.results" % (basename, inter)
    output_log = output_dir + "/%s%s.log" % (basename, inter)
    output_tree = output_dir + "/%s%s.tree" % (basename, inter)

    print(param_file, output_res, output_log, output_tree)

    write_params_file(param_file, aln, tree,output_res, output_log, output_tree, gtr_output=gtr,
                      fixed_beta=fixed_beta, fixed_s=fixed_s, fixed_ProbS = fixed_probS, init_probS=init_probS, init_beta=init_beta, fixed_tau=fixed_tau,
                      fixed_kappa=fixed_kappa, bblOpt=bblOpt)
    #os.system("%s %s" % (path, param_file))
if __name__ == "__main__":
    main()

