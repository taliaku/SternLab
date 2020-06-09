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


SELECTON_WEB_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton/programs/selecton/selecton"
SELECTON_WEB_WRT_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton-wrt/programs/selecton-wrt/selecton"
SELECTON_WEB_WRTM_PATH = "/sternadi/home/volume1/taliakustin/selecton-web/selecton-wrtm/programs/selecton-wrtm/selecton"


def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-r", "--results", dest="results", help="results file")
    parser.add_option("-i", "--input_dir", dest="input_dir", help="dir with input files of best models")
    parser.add_option("-d", "--dir", dest="output_dir", help="output dir path (/sternadi/nobackup/volume1/phyVirus/HIV/selecton)")
    parser.add_option("-n", "--name", dest="output_dir_name", help="output dir name in the /sternadi/nobackup/volume1/phyVirus/HIV/selecton directory")
    parser.add_option("-a", "--action", dest="action_for_script", help="action for script:\n"
                            "w - write best results to file\n"
                            "a - run wrt results on YangModel_small results - fixed theta and rho = 1\n"
                            "i - run wrt result on best results iterate on ")
    parser.add_option("-f", "--fasta", dest="fasta_dir", help="directory with fasta alignemnts with .aln extention")
    parser.add_option("-t", "--tree", dest="tree_dir", default = None, help="directory with trees - will not be used in iterations that take into account last best parameters")


    possible_actions = ["w", "a", "i", "m", "n", "p"]

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
    aln_dir = options.fasta_dir
    tree_dir = options.tree_dir
    aln_dir = check_dirname(aln_dir)
    if tree_dir != None:
        tree_dir = check_dirname(tree_dir)
    if action not in possible_actions:
        raise Exception("actions must be one of: %s" % " ".join(possible_actions))

    print("results file: %s" % results_file)

    results = pd.read_csv(results_file, index_col=None)
    best_lnL = results.groupby("protein")["log-likelihood"].transform(max) == results["log-likelihood"]
    best_results = results[best_lnL].groupby("protein", as_index=False).first()

    if action == "w":
        print(best_results)
        best_results.to_csv("/sternadi/home/volume1/taliakustin/best_results.csv")

    if action == "a":
        small_Yang = results.loc[results["model"] =="YangModel_small"]
        for index, row in best_results.iterrows():
            best_model_path = input_dir + "/" + row["model"]
            protein = row["protein"]

            aln = glob.glob(aln_dir + "/*%s*.aln" % protein)[0]
            tree = glob.glob(tree_dir + "/*%s*tree.txt" % protein)[0]


            if not os.path.exists(output_dir):
                os.system("mkdir %s" % output_dir)
            out = output_dir + "/" + protein

            script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                          "-c %s_color.txt -t %s_tree.txt -v 1 -y 1 -z 3 -j 40 -ft -fr"
                          % (SELECTON_WEB_WRT_PATH, aln, tree, out, out, out, out, out, out),
                          alias="wrt_tr_fixed")





    if action == "i":
        print(best_results)

        #rt_dic = [0.1] + list(frange(0.5, 5.5, 0.5))
        rt_dic = [0.1] + list(frange(0.4, 5.6, 0.5))
        for index, row in best_results.iterrows():
            best_model_path = input_dir + "/" + row["model"]
            protein = row["protein"]
            best_tree = glob.glob(best_model_path + "/%s*tree.txt" % protein)[0]
            aln = glob.glob(aln_dir + "/*" +  protein + "*.aln")[0]
            best_result_file = glob.glob(best_model_path + "/*%s*output.txt" % protein)[0]
            parameters = extract_selecton_final_params_single_file(best_result_file)
            for key in parameters.keys():
                parameters[key] = float(parameters[key])
            print(parameters)
            for i in rt_dic:
                for j in rt_dic:
                    new_dir_GA = output_dir + "_t%s_r%s_GA" % (i, j)
                    new_dir_AG = output_dir + "_t%s_r%s_AG" % (i, j)
                    new_dir_CT = output_dir + "_t%s_r%s_CT" % (i, j)
                    new_dir_TC = output_dir + "_t%s_r%s_TC" % (i, j)


                    #if not os.path.exists(new_dir_GA):
                    #    os.system("mkdir %s" % new_dir_GA)

                    if not os.path.exists(new_dir_AG):
                        os.system("mkdir %s" % new_dir_AG)
                    """

                    if not os.path.exists(new_dir_CT):
                        os.system("mkdir %s" % new_dir_CT)
                    if not os.path.exists(new_dir_TC):
                        os.system("mkdir %s" % new_dir_TC)
                    """

                    out_GA = new_dir_GA + "/" + protein
                    out_AG = new_dir_AG + "/" + protein
                    out_CT = new_dir_CT + "/" + protein
                    out_TC = new_dir_TC + "/" + protein

                    """
                    script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                  "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -z 3"
                                  % (SELECTON_WEB_WRT_PATH, aln, best_tree, out_GA, out_GA, out_GA
                                     , out_GA, out_GA, out_GA, i, j, parameters["additional_omega_category"],
                                     1 - parameters["prob(additional_omega_category)"],
                                     parameters["alpha"], parameters["beta"], parameters["kappa"]), alias="GA_wrt_forth_iteration")
                    """
                    script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                  "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -z 0"
                                  % (SELECTON_WEB_WRT_PATH, aln, best_tree, out_AG, out_AG, out_AG
                                     , out_AG, out_AG, out_AG, i, j, parameters["additional_omega_category"],
                                     1 - parameters["prob(additional_omega_category)"],
                                     parameters["alpha"], parameters["beta"], parameters["kappa"]),
                                  alias="AG_wrt_fourth_iteration")

                    """
                    
                    script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                  "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -v 7"
                                  % (SELECTON_WEB_WRT_PATH, aln, best_tree, out_CT, out_CT, out_CT
                                     , out_CT, out_CT, out_CT, i, j, parameters["additional_omega_category"],
                                     1 - parameters["prob(additional_omega_category)"],
                                     parameters["alpha"], parameters["beta"], parameters["kappa"]),
                                  alias="wrt_first_iteration_CT")

                    script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                  "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -v 11"
                                  % (SELECTON_WEB_WRT_PATH, aln, best_tree, out_TC, out_TC, out_TC
                                     , out_TC, out_TC, out_TC, i, j, parameters["additional_omega_category"],
                                     1 - parameters["prob(additional_omega_category)"],
                                     parameters["alpha"], parameters["beta"], parameters["kappa"]),
                                  alias="wrt_first_iteration_TC")
                    """


    if action == "f":
        print(best_results)

        #rt_dic = [0.1] + list(frange(0.5, 5.5, 0.5))
        r_dic = [0.1] + list(frange(0.4, 5.6, 0.5))
        for index, row in best_results.iterrows():
            best_model_path = input_dir + "/" + row["model"]
            protein = row["protein"]
            best_tree = glob.glob(best_model_path + "/%s*tree.txt" % protein)[0]
            aln = glob.glob(aln_dir + "/*" +  protein + "*.aln")[0]
            best_result_file = glob.glob(best_model_path + "/*%s*output.txt" % protein)[0]
            parameters = extract_selecton_final_params_single_file(best_result_file)
            for key in parameters.keys():
                parameters[key] = float(parameters[key])
            print(parameters)
            for i in rt_dic:
                new_dir_AG = output_dir + "_t%s_r%s_AG" % (i, j)



                #if not os.path.exists(new_dir_GA):
                #    os.system("mkdir %s" % new_dir_GA)

                if not os.path.exists(new_dir_AG):
                    os.system("mkdir %s" % new_dir_AG)


                out_GA = new_dir_GA + "/" + protein
                out_AG = new_dir_AG + "/" + protein
                out_CT = new_dir_CT + "/" + protein
                out_TC = new_dir_TC + "/" + protein


                script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                              "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -z 0"
                              % (SELECTON_WEB_WRT_PATH, aln, best_tree, out_AG, out_AG, out_AG
                                 , out_AG, out_AG, out_AG, i, j, parameters["additional_omega_category"],
                                 1 - parameters["prob(additional_omega_category)"],
                                 parameters["alpha"], parameters["beta"], parameters["kappa"]),
                              alias="AG_wrt_fourth_iteration")



    if action == "m":
        print(best_results)

        rt_dic = [0.1] + list(frange(0.5, 5.5, 0.5))
        for index, row in best_results.iterrows():
            best_model_path = input_dir + "/" + row["model"]
            protein = row["protein"]
            best_tree = glob.glob(best_model_path + "/%s*tree.txt" % protein)[0]
            aln = glob.glob(aln_dir + "/*" +  protein + "*.aln")[0]
            best_result_file = glob.glob(best_model_path + "/*%s*output.txt" % protein)[0]
            parameters = extract_selecton_final_params_single_file(best_result_file)
            for key in parameters.keys():
                parameters[key] = float(parameters[key])
            print(parameters)
            for i in rt_dic:
                for j in rt_dic:
                    for y in rt_dic:
                        new_dir_GA = output_dir + "_t%s_r%s_m%s_GA" % (i, j, y)



                        #if not os.path.exists(new_dir_GA):
                        #    os.system("mkdir %s" % new_dir_GA)

                        if not os.path.exists(new_dir_GA):
                            os.system("mkdir %s" % new_dir_GA)


                        out_GA = new_dir_GA + "/" + protein



                        script_runner("%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                      "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -z 3 -g %s"
                                      % (SELECTON_WEB_WRTM_PATH, aln, best_tree, out_GA, out_GA, out_GA
                                         , out_GA, out_GA, out_GA, i, j, parameters["additional_omega_category"],
                                         1 - parameters["prob(additional_omega_category)"],
                                         parameters["alpha"], parameters["beta"], parameters["kappa"], y),
                                      alias="GA_wrtm_first_iteration")

    if action == "n":
        print(best_results)

        rt_dic = [0.1] + list(frange(0.5, 5.5, 0.5))
        for index, row in best_results.iterrows():
            best_model_path = input_dir + "/" + row["model"]
            protein = row["protein"]
            best_tree = glob.glob(best_model_path + "/%s*tree.txt" % protein)[0]
            aln = glob.glob(aln_dir + "/*" + protein + "*.aln")[0]
            best_result_file = glob.glob(best_model_path + "/*%s*output.txt" % protein)[0]
            parameters = extract_selecton_final_params_single_file(best_result_file)
            for key in parameters.keys():
                parameters[key] = float(parameters[key])
            print(parameters)
            for i in rt_dic:

                new_dir_GA = output_dir + "_m%s_GA" % (i)

                # if not os.path.exists(new_dir_GA):
                #    os.system("mkdir %s" % new_dir_GA)

                if not os.path.exists(new_dir_GA):
                    os.system("mkdir %s" % new_dir_GA)

                out_GA = new_dir_GA + "/" + protein

                script_runner(
                    "%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                    "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -z 3 -g %s"
                    % (SELECTON_WEB_WRTM_PATH, aln, best_tree, out_GA, out_GA, out_GA
                       , out_GA, out_GA, out_GA, parameters["theta"], parameters["rho"], parameters["additional_omega_category"],
                       1 - parameters["prob(additional_omega_category)"],
                       parameters["alpha"], parameters["beta"], parameters["kappa"], i),
                    alias="GA_wrtm_all_best_previous_params")


    if action == "p":
        print(best_results)
        count = 0

        rt_dic = [0.1] + list(frange(0.5, 5.5, 0.5))
        for index, row in best_results.iterrows():
            protein = row["protein"]

            aln = glob.glob(aln_dir + "/*%s*.aln" % protein)[0]
            tree = glob.glob(tree_dir + "/*%s*tree.txt" % protein)[0]

            alpha = [0.6, 1, 1.7]
            beta = [0.7, 1, 2]
            kappa = [2, 3, 4]
            theta = [0.5, 1, 1.5]
            rho = [0.5, 1, 2]
            omega = [1.2, 2, 3.3]
            prob = [0.01, 0.05, 0.2]
            mu = [0.5, 1, 2]
            for a in alpha:
                for x in beta:
                    for k in kappa:
                        for t in theta:
                            for r in rho:
                                for w in omega:
                                    for p in prob:
                                        for g in mu:
                                            
                                            new_dir_GA = output_dir + "_a%s_b%s_k%s_t%s_r%s_w%s_p%s_m%s" % (str(a), str(x), str(k), str(t), str(r),
                                                                                                            str(w), str(p), str(g))


                                            if not os.path.exists(new_dir_GA):
                                                os.system("mkdir %s" % new_dir_GA)

                                            out_GA = new_dir_GA + "/" + protein


                                            script_runner(
                                                "%s -i %s -u %s -l %s_log.txt -r %s_result.txt -o %s_output.txt -s %s_rasmul.txt "
                                                "-c %s_color.txt -t %s_tree.txt -v %s -y %s -w %f -p %f -a %f -x %f -k %f -j 40 -z 3 -g %s"
                                                % (SELECTON_WEB_WRTM_PATH, aln, tree, out_GA, out_GA, out_GA
                                                   , out_GA, out_GA, out_GA, t, r,
                                                  w,
                                                   1 - p,
                                                   a, x, k, g),
                                                alias="GA_wrtm_itrate_on_initial_parameters_third_run")


                                            count += 1


        print(count)





if __name__ == "__main__":
    main()

