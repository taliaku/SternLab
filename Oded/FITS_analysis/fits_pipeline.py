#!/powerapps/share/python-anaconda-3.2019.7/bin/python

"""
@Author: odedkushnir

"""

import pandas as pd
import os
import numpy as np
import sys, argparse
import subprocess
from Context_analysis import Context_analysis_RV
from time import sleep
from FITS_analysis import fits_fitness_united
from FITS_analysis import fits_mutation_united
from FITS_analysis import fits_parameters_pos
from pbs_jobs import create_array_pbs_cmd


def check_pbs(job_id):
    """
    :param job_id: The PBS job id
    :return: "Done!", when the job is done
    """
    status = "Running..."
    try:
        process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
        while process != "":
            process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
            sleep(0.05)
        print("")
    except (subprocess.CalledProcessError):
        process = ""
    if process == "":
        status = "Done"
    return status

def submit(cmdfile):
    cmd = "/opt/pbs/bin/qsub " + cmdfile
    result = os.popen(cmd).read()
    return result.split(".")[0]

def fits_adjustments(data):
    data.sort_values(by=['pos'])
    grouped = data.groupby(['pos', 'gen'])
    df_all = pd.DataFrame(columns=["gen", "base", "freq", "pos"])
    for group in grouped:
        df = pd.concat([group[1]]*2, ignore_index=True)
        df.at[0, 'freq'] = 1 - df.at[1, 'freq']
        df.at[0, 'base'] = 0
        df.at[1, 'base'] = 1
        df_all = pd.concat([df_all, df], ignore_index=True)
    df_all["pos"] = df_all["pos"].astype(int)
    return df_all


def fits_data_construction(input_dir, output_dir, from_passage, to_passage, quality, start_position, without_passage=None):
    all = pd.read_csv(input_dir + "%s_data_mutation.csv" % (quality))
    all = all.loc[all.Pos >= start_position]
    all = all[all["label"] != "RVB14-Next-RNA Control"]
    all = all[all["label"] != "RVB14-p1"]

    # all = all[all["label"] != "%s-RNA Control" % virus]
    # all = all[all["label"] != "%s-%s" % (virus, from_passage)]
    # all = all[all["label"] != "%s-%s" % (virus, from_passage)]
    all["passage"] = all["label"].apply(lambda x: x.split("-")[-1].split("p")[-1])
    all["passage"] = np.where(all["passage"] == "RNA Control", 0, all["passage"])
    all["passage"] = all["passage"].astype(int)
    all = all[all["passage"] >= from_passage]
    all = all[all["passage"] <= to_passage]
    all = all[all["passage"] != without_passage]
    all["Type"] = all["Type"].fillna("NonCodingRegion")
    # all["pval"] = all["pval"].fillna(1)
    # Filtering
    # all["Frequency"] = np.where(all["pval"] > 0.01, 0, all["Frequency"])
    # all["Frequency"] = np.where(all["Prob"] < 0.95, 0, all["Frequency"])

    syn = all[all['Type'] == "Synonymous"]

    df_dict = {"syn": syn, "all": all}

    transitions = ["A>G", "G>A", "U>C", "C>U"]
    for key in df_dict:
        for mutation in transitions:
            df_mutation = Context_analysis_RV.checkKey(df_dict, key)
            df_mutation = df_mutation[df_mutation["Mutation"] == mutation]
            quarte_allelic_mapping = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
            if mutation == "A>G":
                df_mutation["ADAR_like"] = df_mutation.Prev.str.contains('UA') | df_mutation.Prev.str.contains('AA')
                adar = df_mutation[df_mutation["ADAR_like"] == True]
                nonadar = df_mutation[df_mutation["ADAR_like"] == False]
                adar_dict = {"adar": adar, "nonadar": nonadar}
                for adar_key in adar_dict:
                    like = Context_analysis_RV.checkKey(adar_dict, adar_key)
                    like['Base'] = like['Base'].apply(lambda x: quarte_allelic_mapping[x])
                    like['Ref'] = like['Ref'].apply(lambda x: quarte_allelic_mapping[x])
                    like = like.rename(columns={'passage': 'gen', 'Base': 'base', 'Frequency': 'freq', 'Ref': 'ref',
                                                'Read_count': 'read_count', 'Rank': 'rank', 'Pos': 'pos'})
                    like = like[['gen', 'base', 'freq', 'pos']]
                    like = like.sort_values(['pos', 'gen', 'freq'])
                    mutation = mutation.replace(">", "")
                    like.to_csv(output_dir + "%s_mutations_%s_%s.csv" % (key, mutation, adar_key), index=False)

                    final_df = fits_adjustments(like)
                    output_file = output_dir + 'final_%sMutations_sorted_%s_%s.txt' % (key, mutation, adar_key)
                    final_df.to_csv(output_file, index=False, sep="\t")

            df_mutation['Base'] = df_mutation['Base'].apply(lambda x: quarte_allelic_mapping[x])
            df_mutation['Ref'] = df_mutation['Ref'].apply(lambda x: quarte_allelic_mapping[x])
            df_mutation = df_mutation.rename(columns={'passage': 'gen', 'Base': 'base', 'Frequency': 'freq', 'Ref': 'ref',
                                    'Read_count': 'read_count', 'Rank': 'rank', 'Pos': 'pos'})
            df_mutation = df_mutation[['gen', 'base', 'freq', 'pos']]
            df_mutation = df_mutation.sort_values(['pos', 'gen', 'freq'])
            mutation = mutation.replace(">", "")
            output_file = output_dir + '%s_mutations_%s.csv' % (key, mutation)
            df_mutation.to_csv(output_file, index=False)
            final_df = fits_adjustments(df_mutation)
            output_file = output_dir + 'final_%sMutations_sorted_%s.txt' % (key, mutation)
            final_df.to_csv(output_file, index=False, sep="\t")



def main(args):
    virus = args.virus
    passages = args.passages
    without_passage = args.without
    from_passage = int(passages.split("p")[1].split("-")[0])
    to_passage = int(passages.split("-p")[1])
    no_generations = to_passage
    input_dir = args.input_dir # directory contains q23_data_mutation.csv - "/Volumes/STERNADILABHOME$/volume3/okushnir/Cirseq/PV/OPV/39mixMOI/"
    pipline_quality = args.quality
    mutation_type = args.mutation_type
    fits_dir = input_dir + "fits/"
    try:
        os.mkdir(fits_dir)
    except OSError:
        print("Creation of the directory %s failed" % fits_dir)
    else:
        print("Successfully created the directory %s " % fits_dir)

    fits_input_dir = fits_dir + "input/"
    try:
        os.mkdir(fits_input_dir)
    except OSError:
        print("Creation of the directory %s failed" % fits_input_dir)
    else:
        print("Successfully created the directory %s " % fits_input_dir)

    fits_input_pass_dir = fits_input_dir + "%s/" % passages
    try:
        os.mkdir(fits_input_pass_dir)
    except OSError:
        print("Creation of the directory %s failed" % fits_input_pass_dir)
    else:
        print("Successfully created the directory %s " % fits_input_pass_dir)

    mutation_lst = ["CU", "AG", "AG_adar", "AG_nonadar", "GA", "UC"]
    fits_input_dir = input_dir + "fits/input/%s/" % passages

    # 1. Create fits dataset from data_mutation.csv file
    print("Creating fits dataset from data_mutation.csv file...")
    start_pos_dict = {"OPV": 3832, "RVB14": 3635, "CVB3": 3745, "PV": 3833} # start from 2B for all viruses
    start_position = Context_analysis_RV.checkKey(start_pos_dict, virus)
    fits_data_construction(input_dir, fits_input_pass_dir, from_passage, to_passage, pipline_quality, start_position, without_passage)

    # break the dataset to positions datasets
    print("breaking the dataset to positions sub-datasets")

    for mutation in mutation_lst:
        df_mutation = pd.read_csv(fits_input_pass_dir + "final_%sMutations_sorted_%s.txt" % (mutation_type, mutation), sep="\t")
        fits_input_pass_mutation_dir = fits_input_pass_dir + "%s" % mutation
        try:
            os.mkdir(fits_input_pass_mutation_dir)
        except OSError:
            print("Creation of the directory %s failed" % fits_input_pass_mutation_dir)
        else:
            print("Successfully created the directory %s " % fits_input_pass_mutation_dir)
        for position in df_mutation["pos"].iloc[0:-1]:
            df_pos = df_mutation.groupby(["pos"]).get_group(position)
            df_pos.to_csv(fits_input_pass_dir + "%s/final_%sMutations_sorted_%s.txt" % (mutation, mutation_type,  position), index=False,
                              sep="\t")
    # 2. Run FITS_jobarray_mutation.cmd
    if virus == "OPV":
        no_generations = (to_passage*2)+1
    with open(fits_input_dir+"parameters_mutation.txt", "w") as parameter_mutation:
        parameter_mutation.write("# basic parameters\n"
                                 "N 1000000\n"
                                 "num_alleles 2\n"
                                 "num_generations %s\n"
                                 "min_log_mutation_rate0_0 -9\n"
                                 "min_log_mutation_rate0_1 -9\n"
                                 "min_log_mutation_rate1_0 -9\n"
                                 "min_log_mutation_rate1_1 -9\n"
                                 "max_log_mutation_rate0_0 -3\n"
                                 "max_log_mutation_rate0_1 -3\n"
                                 "max_log_mutation_rate1_0 -3\n"
                                 "max_log_mutation_rate1_1 -3\n"
                                 "# bottleneck\n"
                                 "bottleneck_interval 2\n"
                                 "bottleneck_size 1000000\n"
                                 "fitness_allele0 1.0\n"
                                 "fitness_allele1 1.0\n"
                                 "num_samples_from_prior 100000\n"
                                 "acceptance_rate 0.01" % (no_generations))
    try:
        os.mkdir(input_dir + "fits/output/")
    except OSError:
        print("Creation of the directory %s failed" % input_dir + "fits/output/")
    else:
        print("Successfully created the directory %s " % input_dir + "fits/output/")
    try:
        os.mkdir(input_dir + "fits/output/mutation/")
    except OSError:
        print("Creation of the directory %s failed" % input_dir + "fits/output/mutation/")
    else:
        print("Successfully created the directory %s " % input_dir + "fits/output/mutation/")

    fits_mutation_output_dir = input_dir + "fits/output/mutation/%s/" % passages
    try:
        os.mkdir(fits_mutation_output_dir)
    except OSError:
        print("Creation of the directory %s failed" % fits_mutation_output_dir)
    else:
        print("Successfully created the directory %s " % fits_mutation_output_dir)

    cmds = "module load gcc/gcc-8.2.0\n" \
           "VIRUS=$VIRUS\n" \
           "PASSAGES=$PASSAGES\n" \
           "fits_dir=$fits_dir\n" \
           "PARAM='mutation'\n" \
           "namesj=('AG' 'UC' 'GA' 'CU' 'AG_adar' 'AG_nonadar')\n"\
           "j=$(($[PBS_ARRAY_INDEX-1]%${#namesj[@]}))\n"\
           "i=$(($[PBS_ARRAY_INDEX-1-j]/${#namesj[@]}))\n" \
           "element=$(ls ${fits_dir}/input/${PASSAGES}/${namesj[j]}| grep 'syn'| grep -oh '[0-9]'*)\n"\
           "namesi=($element)\n" \
           "for element in '${namesi[@]}';do echo '$element';done\n" \
           "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}\n" \
           "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}\n" \
           "/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits1.3.3 -$PARAM ${fits_dir}/input/${PASSAGES}/" \
            "parameters_${PARAM}.txt ${fits_dir}/input/${PASSAGES}/${namesj[j]}/final_synMutations_sorted_" \
            "${namesi[i]}.txt ${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/posterior_${PARAM}_syn_${namesi[i]}.txt " \
            "${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/summary_${PARAM}_syn_${namesi[i]}.txt"
    part_lst = []
    for mutation in mutation_lst:
        list = os.listdir(fits_input_dir + mutation) # dir is your directory path
        number_files = len(list)
        part_lst.append(number_files)
    n = max(part_lst)
    print(n)
    a1 = 1
    d = len(mutation_lst)
    jnum = a1 + d * (n - 1)
    print(jnum)

    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/new_FITS_jobarray_mutation_%s.cmd" % virus
    create_array_pbs_cmd(cmd_file, jnum, alias="Fits_mutation", gmem=3, cmds=cmds)
    print("qsub -v VIRUS='%s',PASSAGES='%s', fits_dir='%s' %s" % (virus, passages, fits_dir, cmd_file))
    job_id = submit("-v VIRUS='%s',PASSAGES='%s, fits_dir=%s' %s" % (virus, passages, fits_dir, cmd_file))
    job_id = job_id.split("[")[0]
    print("Running new_FITS_jobarray_mutation_%s.cmd, job_id:%s" % (virus, job_id))
    status = check_pbs(job_id)
    if status == "Done":
        print("Done!")

    ##For all position at once
    ## job_id = submit("/sternadi/home/volume3/okushnir/Cluster_Scripts/FITS_jobarray_mutation_%s.cmd" % virus)
    ## job_id = job_id.split("[")[0]
    ## print("Running FITS_jobarray_mutation.cmd, job_id:%s"%job_id)
    ## status = check_pbs(job_id)

    # 3. Run fits_mutation_united
        for mutation in mutation_lst:
            output_mutation_dir = input_dir + "fits/output/mutation/%s/%s/" % (passages, mutation)
            output_file = output_mutation_dir + "all.txt"
            print("Creating the mutation conjugated report of: %s" % output_mutation_dir)
            fits_mutation_united(output_mutation_dir, output_file)

    # 4. Run fits_parameters.py
    print("Creating fitness parameters")
    fits_parameters_pos.fitness_parameters_pos(input_dir=fits_mutation_output_dir, output_dir=fits_input_dir)

    # 5. Run new_FITS_jobarray_fitness.cmd
    try:
        os.mkdir(input_dir + "fits/output/fitness/")
    except OSError:
        print("Creation of the directory %s failed" % (input_dir + "fits/output/fitness/"))
    else:
        print("Successfully created the directory %s " % (input_dir + "fits/output/fitness/"))

    if mutation_type == "all":
        cmds = "module load gcc/gcc-8.2.0\n" \
               "VIRUS=$VIRUS\n" \
               "PASSAGES=$PASSAGES\n" \
               "fits_dir=$fits_dir\n" \
               "PARAM='fitness'\n" \
               "namesj=('AG' 'UC' 'GA' 'CU' 'AG_adar' 'AG_nonadar')\n"\
               "j=$(($[PBS_ARRAY_INDEX-1]%${#namesj[@]}))\n"\
               "i=$(($[PBS_ARRAY_INDEX-1-j]/${#namesj[@]}))\n" \
               "element=$(ls ${fits_dir}/input/${PASSAGES}/${namesj[j]}| grep 'all'| grep -oh '[0-9]'*)\n"\
               "namesi=($element)\n" \
               "for element in '${namesi[@]}';do echo '$element';done\n" \
               "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}\n" \
               "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}\n" \
               "/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits1.3.3 -$PARAM ${fits_dir}/input/${PASSAGES}/" \
                "parameters_${PARAM}_${namesj[j]}.txt ${fits_dir}/input/${PASSAGES}/${namesj[j]}/final_allMutations_sorted_" \
                "${namesi[i]}.txt ${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/posterior_${PARAM}_all_${namesi[i]}.txt " \
                "${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/summary_${PARAM}_all_${namesi[i]}.txt"
    if mutation_type == "syn":
        cmds = "module load gcc/gcc-8.2.0\n" \
               "VIRUS=$VIRUS\n" \
               "PASSAGES=$PASSAGES\n" \
               "fits_dir=$fits_dir\n" \
               "PARAM='fitness'\n" \
               "namesj=('AG' 'UC' 'GA' 'CU' 'AG_adar' 'AG_nonadar')\n"\
               "j=$(($[PBS_ARRAY_INDEX-1]%${#namesj[@]}))\n"\
               "i=$(($[PBS_ARRAY_INDEX-1-j]/${#namesj[@]}))\n" \
               "element=$(ls ${fits_dir}/input/${PASSAGES}/${namesj[j]}| grep 'syn'| grep -oh '[0-9]'*)\n"\
               "namesi=($element)\n" \
               "for element in '${namesi[@]}';do echo '$element';done\n" \
               "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}\n" \
               "mkdir ${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}\n" \
               "/sternadi/home/volume1/talzinger/FITS_Analyses/FITS_bin/fits1.3.3 -$PARAM ${fits_dir}/input/${PASSAGES}/" \
                "parameters_${PARAM}_${namesj[j]}.txt ${fits_dir}/input/${PASSAGES}/${namesj[j]}/final_synMutations_sorted_" \
                "${namesi[i]}.txt ${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/posterior_${PARAM}_syn_${namesi[i]}.txt " \
                "${fits_dir}/output/${PARAM}/${PASSAGES}/${namesj[j]}/summary_${PARAM}_syn_${namesi[i]}.txt"
    part_lst = []
    for mutation in mutation_lst:
        list = os.listdir(fits_input_dir + mutation) # dir is your directory path
        number_files = len(list)
        part_lst.append(number_files)
    n = max(part_lst)
    print(n)
    a1 = 1
    d = len(mutation_lst)
    jnum = a1 + d * (n - 1)
    print(jnum)

    cmd_file = "/sternadi/home/volume3/okushnir/Cluster_Scripts/new_FITS_jobarray_fitness_%s.cmd" % virus
    create_array_pbs_cmd(cmd_file, jnum, alias="Fits_fitness", gmem=3, cmds=cmds)
    print("qsub -v VIRUS='%s',PASSAGES='%s', fits_dir='%s' %s" % (virus, passages, fits_dir, cmd_file))
    job_id = submit("-v VIRUS='%s',PASSAGES='%s, fits_dir=%s' %s" % (virus, passages, fits_dir, cmd_file))
    job_id = job_id.split("[")[0]
    print("Running new_FITS_jobarray_fitness_%s.cmd, job_id:%s" % (virus, job_id))
    status = check_pbs(job_id)
    if status == "Done":
    # 6. Run fits_fitness_united
        for mutation in mutation_lst:
            output_fitness_dir = input_dir + "fits/output/fitness/%s/%s/" % (passages, mutation)
            output_file = output_fitness_dir + "all.txt"
            print("Creating the fitness conjugated report of: %s" % output_fitness_dir)
            fits_fitness_united(output_fitness_dir, output_file, mutation_type)
    # 7. Run fits_plotter.py

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("virus", type=str, help="name of the virus, RVB14")
    parser.add_argument("passages", type=str, help="from which passages, p0-p12")
    parser.add_argument("without", type=int, help="Exclude passage no.")
    parser.add_argument("input_dir", type=str, help="the path to the directory that contains data_mutation.csv")
    parser.add_argument("quality", type=str, help="what is the prefix for the data_mutation.csv file; quality of the pipline ; for example: q38")
    parser.add_argument("mutation_type", type=str, help="all/syn")
    args = parser.parse_args(sys.argv[1:])
    main(args)
