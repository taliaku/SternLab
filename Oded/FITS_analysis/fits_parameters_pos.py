#! /usr/local/python/anaconda_python-3.6.1

"""
@Author: odedkushnir

"""


import glob
import pandas as pd
import subprocess
import os
from Context_analysis_RV import checkKey


def fitness_parameters_pos(input_dir, output_dir):
    """
    Creates the fitness parameters file for each mutation
    :param input_dir: the output directory for fitness -mutation
    :param output_dir: the input directory for FITS
    :return: Dictionary for all mutation and their mutation rates
    """
    mutation_lst = ["AG", "AG_adar", "AG_nonadar", "UC", "GA", "CU"]
    rate_dict = {}
    for mutation in mutation_lst:
        file = input_dir + "/all.txt"
        data = pd.read_csv(file, sep="\t")
        data["inferred_mu"] = data["inferred_mu"].map(lambda x: str(x).lstrip('*'))
        data["inferred_mu"] = pd.to_numeric(data["inferred_mu"], errors='coerce')
        rate = data["inferred_mu"].median()
        data["rev_inferred_mu"] = data["rev_inferred_mu"].map(lambda x: str(x).lstrip('*'))
        data["rev_inferred_mu"] = pd.to_numeric(data["rev_inferred_mu"], errors='coerce')
        adar_rev = data["rev_inferred_mu"].median()
        rate_dict[mutation] = rate
        if (mutation == "AG_adar") | (mutation == "AG_nonadar"):
            rev_mutation = mutation[::-1]
            rate_dict[rev_mutation] = adar_rev
        print(rate_dict)

    for mutation in rate_dict:
        if (mutation == "AG_adar") | (mutation == "AG_nonadar"):
            with open(output_dir + "parameters_fitness_%s.txt" % mutation, "w+") as out_handler:
                out_handler.writelines("verbose 1\n"
                                       "N 756000\n"
                                       "fitness_prior smoothed_composite\n"
                                       "mutation_rate0_1 %s\n"
                                       "mutation_rate1_0 %s\n"
                                       "min_fitness_allele0 1.0\n"
                                       "max_fitness_allele0 1.0\n"
                                       "min_fitness_allele1 0.0\n"
                                       "max_fitness_allele1 2.0\n"
                                       "num_samples_from_prior 100000\n"
                                       "acceptance_rate 0.01" % (rate_dict["AG_nonadar"], rate_dict["radanon_GA"]))
        elif (mutation != "rada_GA") & (mutation != "radanon_GA"):
            rev_mutation = mutation[::-1]
            with open(output_dir + "parameters_fitness_%s.txt" % mutation, "w+") as out_handler:
                out_handler.writelines("verbose 1\n"
                                        "N 756000\n"
                                        "fitness_prior smoothed_composite\n"
                                        "mutation_rate0_1 %s\n"
                                        "mutation_rate1_0 %s\n"
                                        "min_fitness_allele0 1.0\n"
                                        "max_fitness_allele0 1.0\n"
                                        "min_fitness_allele1 0.0\n"
                                        "max_fitness_allele1 2.0\n"
                                        "num_samples_from_prior 100000\n"
                                        "acceptance_rate 0.01" % (rate_dict[mutation], rate_dict[rev_mutation]))


    return rate_dict

def main():
    virus = "PV"
    passages = "p3-p8"

    input_dir = "/Users/odedkushnir/Projects/fitness/test/"
    output_dir = "/Users/odedkushnir/Projects/fitness/test"

    fitness_parameters_pos(input_dir, output_dir)

if __name__ == "__main__":
    main()
