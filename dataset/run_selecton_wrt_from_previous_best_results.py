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

SELECTON_PATH_T_AND_THAN_R = "/sternadi/home/volume1/taliakustin/phylogenyCode/t_and_than_r/programs/selecton/selecton"
SELECTON_PATH_R_AND_THAN_T = "/sternadi/home/volume1/taliakustin/phylogenyCode/r_and_than_t/programs/selecton/selecton"
A_TO_G_SELECTON = "/sternadi/home/volume1/taliakustin/phylogenyCode/AG_selecton/programs/selecton/selecton"
SELECTON_WEB_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton/programs/selecton/selecton"


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-r", "--results", dest="results", help="results file")
    parser.add_option("-i", "--input_dir", dest="input_dir", help="dir with input files of best models")
    parser.add_option("-d", "--dir", dest="output_dir", help="output dir path (/sternadi/nobackup/volume1/phyVirus/HIV/selecton)")
    parser.add_option("-n", "--name", dest="output_dir_name", help="output dir name in the /sternadi/nobackup/volume1/phyVirus/HIV/selecton directory")
    parser.add_option("-a", "--action", dest="action_for_script", help="action for script:"
                                                                       "r - run selecton for each best parameter file with all r and t combinations,"
                                                                       "t - run selecton for each best parameter file with all best parameters,"
                                                                       "f - run selecton with 1 parameter fixed (t or r)" 
                                                                       "g - run the A>G selecton" 
                                                                       "e - run selecton web from non-selecton-web best results"
                                                                       "a - run web selecton for each best parameter file with all r and t combinations" 
                                                                       "w - write")

    possible_actions = ["r", "w", "t", "f", "g", "e"]
    (options, args) = parser.parse_args()
    results_file = options.results
    dir = options.output_dir
    input_dir = options.input_dir
    name = options.output_dir_name
    results_file = check_filename(results_file)
    dir = check_dirname(dir)
    input_dir = check_dirname(input_dir)
    output_dir = dir + "/" +  name
    output_dir = check_dirname(output_dir, Truedir=False)
    action = options.action_for_script
    if action not in possible_actions:
        raise Exception("actions must be one of: %s" % " ".join(possible_actions))

    print(output_dir)
    print("results file: %s" % results_file)

    results = pd.read_csv(results_file, index_col=None)
    best_lnL = results.groupby("Protein")["log-likelihood"].transform(max) == results["log-likelihood"]
    best_results = results[best_lnL].groupby("Protein", as_index=False).first()

    print(best_results)
    if action == "w":
        print(best_results)
        best_results.to_csv("/sternadi/home/volume1/taliakustin/best_results.csv")




    if action == "r" or action == "g":
        rt_dic = [0.1] + list(frange(0.5, 5.5, 0.5))
        count = 0
        for index, row in best_results.iterrows():
            #print(row["model"])
            best_model_path = input_dir + "/" + row["model"]
            protein = row["Protein"]

            best_tree = glob.glob(best_model_path + "/%s*tree.txt" % protein)[0]
            #print(dir + "/*" + protein + "*.aln")
            aln = glob.glob(input_dir + "/codon_aln/*" + protein + "*.aln")[0]
            #print (aln)
            best_result_file = glob.glob(best_model_path + "/*%s*output.txt" %protein)[0]
            print(best_result_file)

            (likelihood, kappa, theta, rho, alpha, beta, add_omega, prob) = extract_selecton_final_params_single_file(best_result_file)

            print(best_result_file, aln, best_tree)

            for i in rt_dic:
                for j in rt_dic:
                    print(i, j)
                    count += 1

                    if action == "r":
                        new_dir_first_t =output_dir + "_t%s_r%s_first_t" % (i, j)
                        new_dir_first_r = output_dir + "_t%s_r%s_first_r" % (i, j)
                        if not os.path.exists(new_dir_first_t):
                            os.system("mkdir %s" % new_dir_first_t)
                        if not os.path.exists(new_dir_first_r):
                            os.system("mkdir %s" % new_dir_first_r)

                        out_first_t = new_dir_first_t + "/" + protein
                        out_first_r = new_dir_first_r + "/" + protein

                        script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                      "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40"
                                      % (SELECTON_PATH_T_AND_THAN_R, aln, best_tree, out_first_t, out_first_t, out_first_t
                                         , out_first_t, out_first_t, out_first_t, i, j, add_omega, 1-prob, alpha, beta, kappa), alias="second_t")
                        script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                      "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40"
                                      % (SELECTON_PATH_R_AND_THAN_T, aln, best_tree, out_first_r, out_first_r, out_first_r
                                         , out_first_r, out_first_r, out_first_r, i, j, add_omega, 1 - prob, alpha, beta, kappa), alias="second_r")
                    if action == "g":
                        new_dir = output_dir + "_t%s_r%s" % (i, j)
                        if not os.path.exists(new_dir):
                            os.system("mkdir %s" % new_dir)


                        out = new_dir + "/" + protein

                        script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                      "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40"
                                      % (
                                      A_TO_G_SELECTON, aln, best_tree, out, out, out
                                      , out, out, out, i, j, add_omega, 1 - prob, alpha, beta,
                                      kappa), alias="A_to_G")



    if action == "f":
        rt_dic = [0.1] + list(frange(0.5, 5.5, 0.5))
        count = 0
        for index, row in best_results.iterrows():
            best_model_path = dir + "/" + row["model"]
            protein = row["Protein"]
            best_tree = glob.glob(best_model_path + "/%s*tree.txt" % protein)[0]

            # print(dir + "/*" + protein + "*.aln")
            aln = glob.glob(dir + "/codon_aln/*" + protein + "*.aln")[0]
            # print (aln)
            best_result_file = glob.glob(best_model_path + "/*%s*output.txt" % protein)[0]
            print(best_result_file)

            (likelihood, kappa, theta, rho, alpha, beta, add_omega,
             prob) = extract_selecton_final_params_single_file(best_result_file)

            print(best_result_file, aln, best_tree)

            for i in rt_dic:
                for j in rt_dic:
                    print(i, j)
                    count += 1

                    new_dir_fixed_t = output_dir + "_t%s_r%s_fixed_t" % (i, j)
                    if not os.path.exists(new_dir_fixed_t):
                        os.system("mkdir %s" % new_dir_fixed_t)


                    new_dir_fixed_r = output_dir + "_t%s_r%s_fixed_r" % (i, j)
                    if not os.path.exists(new_dir_fixed_r):
                        os.system("mkdir %s" % new_dir_fixed_r)

                    out_fixed_t = new_dir_fixed_t + "/" + protein
                    out_fixed_r = new_dir_fixed_r + "/" + protein


                    script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                  "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40 -ft"
                                  % (
                                  SELECTON_PATH_T_AND_THAN_R, aln, best_tree, out_fixed_t, out_fixed_t, out_fixed_t
                                  , out_fixed_t, out_fixed_t, out_fixed_t, i, j, add_omega, 1 - prob, alpha, beta,
                                  kappa), alias="fixed_t")
                    script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                  "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40 -fr"
                                  % (
                                      SELECTON_PATH_T_AND_THAN_R, aln, best_tree, out_fixed_r, out_fixed_r, out_fixed_r
                                  , out_fixed_r, out_fixed_r, out_fixed_r, i, j, add_omega, 1 - prob, alpha, beta,
                                  kappa), alias="fixed_r")


    if action == "t" or action == "e":
        for index, row in best_results.iterrows():
            #print(row["model"])
            best_model_path = dir + "/" + row["model"]
            protein = row["Protein"]
            best_model_path = "/sternadi/home/volume1/taliakustin/HIV_selecton/YangModel"
            best_tree = glob.glob(best_model_path + "/%s*tree.txt" % protein)[0]

            #print(dir + "/*" + protein + "*.aln")
            aln = glob.glob(dir + "/codon_aln/*" + protein + "*.aln")[0]
            #print (aln)
            best_result_file = glob.glob(best_model_path + "/*%s*output.txt" %protein)[0]
            print(best_result_file)

            (likelihood, kappa, theta, rho, alpha, beta, add_omega, prob) = extract_selecton_final_params_single_file(best_result_file)

            print(best_result_file, aln, best_tree)

            if action == "t":
                new_dir_first_t = output_dir + "_first_t"
                new_dir_first_r = output_dir + "_first_r"


                if not os.path.exists(new_dir_first_t):
                    os.system("mkdir %s" % new_dir_first_t)
                if not os.path.exists(new_dir_first_r):
                    os.system("mkdir %s" % new_dir_first_r)

                out_first_t = new_dir_first_t + "/" + protein
                out_first_r = new_dir_first_r + "/" + protein



                script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                              "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40"
                              % (SELECTON_PATH_T_AND_THAN_R, aln, best_tree, out_first_t, out_first_t, out_first_t
                                 , out_first_t, out_first_t, out_first_t, theta, rho, add_omega, 1-prob, alpha, beta, kappa), alias="best_t")
                script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                              "-c %s_color.txt -t %s_tree.txt -y %s -z %s -w %f -p %f -a %f -x %f -k %f -j 40"
                              % (SELECTON_PATH_R_AND_THAN_T, aln, best_tree, out_first_r, out_first_r, out_first_r
                                 , out_first_r, out_first_r, out_first_r, theta, rho, add_omega, 1 - prob, alpha, beta, kappa), alias="best_r")
            elif action == "e":
                if not os.path.exists(output_dir):
                    os.system("mkdir %s" % output_dir)
                out = output_dir + "/" + protein

                script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                              "-c %s_color.txt -t %s_tree.txt -w %f -p %f -a %f -x %f -k %f -j 40"
                              % (SELECTON_WEB_PATH, aln, best_tree, out, out, out
                                 , out, out, out, add_omega, 1 - prob, alpha,
                                 beta, kappa), alias="web_selecton_best_non_web_selecton_params")


if __name__ == "__main__":
    main()

