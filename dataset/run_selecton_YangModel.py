#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename, check_dirname
from PAML_utilities import write_ctl_file
import pandas as pd
import glob
from selecton_utilities import extract_selecton_final_params_single_file
from pbs_runners import script_runner
import os
from general_utilities import frange



SELECTON_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton/programs/selecton/selecton"

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--dir", dest="output_dir", help="output dir path (/sternadi/nobackup/volume1/phyVirus/HIV/selecton)")
    parser.add_option("-n", "--name", dest="output_dir_name", help="output dir name in the /sternadi/nobackup/volume1/phyVirus/HIV/selecton directory")
    parser.add_option("-a", "--aln", dest="alignment_dir", help="directory with alignemnts")
    parser.add_option("-t", "--tree", dest="trees_dir", help="directory with trees")
    parser.add_option("--action", dest="action_for_script", help="action for script:"
                        "f - YangModel run without any determination of parameters"
                        "s - YangModel run with parameter initialization")

    possible_actions = ["f", "s"]



    (options, args) = parser.parse_args()
    dir = options.output_dir
    name = options.output_dir_name
    alignment_dir = options.alignment_dir
    trees_dir = options.trees_dir
    action = options.action_for_script


    dir = check_dirname(dir)
    alignment_dir = check_dirname(alignment_dir)
    trees_dir = check_dirname(trees_dir)
    output_dir = dir + "/" +  name
    output_dir = check_dirname(output_dir, Truedir=False)

    if action not in possible_actions:
        raise Exception("actions must be one of: %s" % " ".join(possible_actions))



    alns = glob.glob(alignment_dir + '/*.aln')
    trees = [glob.glob(trees_dir + "/" + os.path.basename(a).split(".aln")[0].split("_no_stop")[0] + "*tree.txt")[0] for a in alns]

    new_dir = output_dir
    if not os.path.exists(new_dir):
        os.system("mkdir %s" % new_dir)

    if action =="f":


        for i in range(len(alns)):
            aln = alns[i]
            tree = trees[i]
            protein = os.path.basename(aln).split("_")[3]

            out = new_dir + "/" + protein
            script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                          "-c %s_color.txt -t %s_tree.txt -j 40"
                          % (SELECTON_PATH, aln, tree, out, out, out, out, out, out), alias="YangModelsmalls")

    elif action == "s":
        count = 0
        kappa = [1.8, 3.3, 4.4]
        alpha = [0.2, 0.7, 1.4]
        beta = [0.8, 1.7, 2.7]
        additional_omega = [1.1, 1.5, 2.3, 3.1]
        prob =[0.01, 0.1, 0.4]
        for k in kappa:
            for a in alpha:
                for b in beta:
                    for w in additional_omega:
                        for p in prob:
                            new_dir = output_dir + '_k%s_a%s_b%s_w%s_p%s' % (k, a, b, w, p)
                            if not os.path.exists(new_dir):
                                os.system("mkdir %s" % new_dir)

                            for i in range(len(alns)):
                                count += 1

                                aln = alns[i]
                                tree = trees[i]
                                protein = os.path.basename(aln).split("_")[3]

                                out = new_dir + "/" + protein
                                script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                              "-c %s_color.txt -t %s_tree.txt -w %f -p %f -a %f -x %f -k %f -j 40"
                                              % (SELECTON_PATH, aln, tree, out, out, out, out, out, out,
                                                 w, 1 - p, a, b, k), alias="YangModelSmallIterations")






if __name__ == "__main__":
    main()

