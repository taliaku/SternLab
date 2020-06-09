#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
import os
import itertools
from selecton_utilities import extract_selecton_final_params_single_file
import glob


SELECTON_WEB_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton/programs/selecton/selecton"
SELECTON_WEB_WRT_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton-wrt/programs/selecton-wrt/selecton"
SELECTON_WEB_WRTM_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton-wrtm/programs/selecton-wrtm/selecton"


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("--idx", dest="idx", help="index number")
    parser.add_option("-a", "--aln", dest="aln", help="alignment file")
    parser.add_option("-b", "--best_model", dest="best_model", help="best previous model file")
    parser.add_option("-o", "--output", dest="output", help="output directory")
    parser.add_option("-p", "--protein", dest="protein", help="protein name")
    (options, args) = parser.parse_args()


    idx = int(options.idx)
    aln = options.aln
    best_model_file = options.best_model
    output = options.output
    protein = options.protein



    parameters = extract_selecton_final_params_single_file(best_model_file)

    parameter = [0.1, 0.3, 0.5, 0.7, 0.9, 1, 1.1, 1.3, 1.5, 2, 3, 4]


    tree = glob.glob(os.path.dirname(best_model_file) + "/%s*tree*" % protein)[0]

    a = float(parameters["alpha"])
    x = float(parameters["beta"])
    k = float(parameters["kappa"])
    w = float(parameters["additional_omega_category"])
    p = float(parameters["prob(additional_omega_category)"])

    if idx <= 12:
        print("parameter checked is theta")
        print("mu is fixed to 1")
        t = parameter[idx - 1]
        output_GA = output + "wrtmModel_fixed1_mr_theta%s/" % str(t)
        if not os.path.exists(output_GA):
            os.system("mkdir %s" % output_GA)
        output_GA_protein = output_GA + protein
        output_path = "%s_output.txt" % output_GA_protein
        print(output_GA, output_path)
        command = "%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt " \
                  "-c %s_color.txt -t %s_tree.txt -v %s -y 1 -w %f -p %f -a %f -x %f -k %f -j 40 -z 3 -g 1 -fm -fr" % \
                                                    (SELECTON_WEB_WRTM_PATH, aln, tree, output_GA_protein, output_GA_protein, output_GA_protein
                                                       , output_GA_protein, output_GA_protein, output_GA_protein, t, w, 1 - p, a, x, k)


    elif idx > 12:
        idx = idx - 12
        print("parameter checked is mu")
        print("theta is fixed to 1")
        g = parameter[idx - 1]
        output_GA = output + "wrtmModel_fixed1_tr_mu%s/" % str(g)
        if not os.path.exists(output_GA):
            os.system("mkdir %s" % output_GA)
        output_GA_protein = output_GA + protein
        output_path = "%s_output.txt" % output_GA_protein
        command = "%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt " \
                  "-c %s_color.txt -t %s_tree.txt -v 1 -y 1 -w %f -p %f -a %f -x %f -k %f -j 40 -z 3 -g %f -ft -fr" % \
                                                    (SELECTON_WEB_WRTM_PATH, aln, tree, output_GA_protein, output_GA_protein, output_GA_protein
                                                       , output_GA_protein, output_GA_protein, output_GA_protein, w, 1 - p, a, x, k, g)


    if (os.path.exists(output_path) and os.stat(output_path).st_size == 0) or not os.path.exists(output_path):
        print(command)
        os.system(command)

    else:
        print("output file already exists - didn't run selecton")



if __name__ == "__main__":
    main()

