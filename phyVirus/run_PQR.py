#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file_PQR

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

    print(model)

    if model == "alternative":
        inter = "PQR_alt"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/PQR_alt"
        fixed_beta = 0
        fixed_s = 1
        fixed_probS = 1
        init_probS = 0
        fixed_tau = 0
        path = "/sternadi/home/volume1/taliakustin/software/PQR_no_root/programs/PQRmodel/PQRmodel"
    elif model == "null":
        inter = "PQR_null"
        output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/PQR_null/"
        fixed_beta = 0
        fixed_s = 0
        fixed_probS = 0
        init_probS = 0.01
        fixed_tau = 0
        path = "/sternadi/home/volume1/taliakustin/software/PQR_no_root/programs/PQRmodel/PQRmodel"

    basename = aln.split("/")[-1].split("aln.best.fas")[0]

    param_file = output_dir + "/%s%s.params" % (basename, inter)

    output_res = output_dir + "/%s%s.results" % (basename, inter)
    output_log = output_dir + "/%s%s.log" % (basename, inter)
    output_tree = output_dir + "/%s%s.tree" % (basename, inter)


    print(param_file, output_res, output_log, output_tree)

    write_params_file_PQR(param_file, aln, tree,output_res, output_log, output_tree, gtr_output=gtr,
                      fixed_beta=fixed_beta, fixed_s=fixed_s, fixed_ProbS = fixed_probS, init_probS=init_probS, fixed_tau=fixed_tau)
    os.system("%s %s" % (path, param_file))

if __name__ == "__main__":
    main()

