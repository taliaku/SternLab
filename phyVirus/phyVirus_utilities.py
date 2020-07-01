#! /usr/local/python_anaconda/bin/python3.4

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')

import pandas as pd
from file_utilities import check_filename, check_dirname
from phyVirus.get_baltimore import get_baltimore_classifiaction
from phylogenetic_utilities import total_branch_length, average_branch_length, \
    maximum_branch_length, get_average_tip_to_root_distance, root_at_midpoint, get_median_branch_length,\
    get_median_leaf_length, get_branch_length_info, get_number_of_nodes
from Bio import SeqIO
import glob
import pbs_runners
from collections import Counter
import os
from seqFileTools import convert_fasta_to_phylip
from PAML_utilities import *
import seqFileAnalyzer
import seqFileTools
import fastml_utilities
import tqdm


def generate_sequence_count_file(files, output):
    check_filename(output, Truefile=False)
    seq_count = pd.DataFrame()
    for f in files:
        base = f.split("/")[-1].split(".aln.best.fas")[0].split(".fasta")[0]
        family = base.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        count = seqFileTools.count_seqs_in_fasta(f)
        seq_count = seq_count.append({"baltimore": baltimore, "family": family, "basename": base, "count": count},
                                     ignore_index=True)
    seq_count.to_csv(output)


def get_tree_statistics(files, output):
    check_filename(output, Truefile=False)
    df = pd.DataFrame()
    for f in tqdm.tqdm(files):
        if "check" in f:
            continue
        base = f.split("/")[-1].split(".aln.best.fas")[0].split(".fasta")[0].split(".aln.best.phy")[0]
        family = base.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        total_branch = total_branch_length(f)
        average_branch = average_branch_length(f)
        maximum_branch = maximum_branch_length(f)
        average_tip_to_root_distance = get_average_tip_to_root_distance(f)
        median_branch_length = get_median_branch_length(f)
        median_leaf_length = get_median_leaf_length(f)
        node_number = get_number_of_nodes(f)
        df = df.append({"baltimore": baltimore, "family": family, "base": base,
                       "total_branch":total_branch, "average_branch":average_branch,
                       "maximum_branch":maximum_branch, "tip_to_root_distance":average_tip_to_root_distance,
                        "median_branch_length":median_branch_length, "median_leaf_length":median_leaf_length, "nodes":node_number},
                       ignore_index=True)

    df.to_csv(output)
    return(df)

def get_fasta_statistics():
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta/*fasta")
    df = pd.DataFrame()
    for f in tqdm.tqdm(files):
        if "check" in f:
            continue
        basename = f.split("/")[-1].split(".fasta")[0]
        family = basename.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        seq_num = seqFileTools.get_number_of_seqs(f)
        if len(glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/" + "codon_aln/%s*phy" % basename)) == 1:
            codon_aln = True
        else:
            codon_aln = False
        df = df.append({"base":basename, "family":family, "baltimore":baltimore, "nodes":seq_num, "codon_aln":codon_aln}, ignore_index=True)
    df.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/sequence_info.csv")


def get_alignment_statistics():
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln/*fas")
    df = pd.DataFrame()
    for f in tqdm.tqdm(files):
        if "check" in f:
            continue
        basename = f.split("/")[-1].split(".fasta")[0]
        family = basename.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        seq_num = seqFileTools.get_number_of_seqs(f)
        if len(glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/" + "codon_aln/%s*phy" % basename)) == 1:
            codon_aln = True
        else:
            codon_aln = False
        df = df.append({"base":basename, "family":family, "baltimore":baltimore, "nodes":seq_num, "codon_aln":codon_aln}, ignore_index=True)
    df.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_info.csv")



def get_branch_length_info_all_trees(files, output):
    check_filename(output, Truefile=False)
    output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/branch_length_info/"
    df = pd.DataFrame()
    for f in tqdm.tqdm(files):
        if "check" in f:
            continue
        base = f.split("/")[-1].split(".aln.best.fas")[0].split(".fasta")[0].split(".aln.best.phy")[0]
        family = base.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        output_file = output_dir + base + ".branch_lengths_info.csv"
        temp_df = get_branch_length_info(f)
        temp_df["base"] = base
        temp_df["family"] = family
        temp_df["baltimore"] = baltimore
        temp_df.to_csv(output_file)
        df = df.append(temp_df, ignore_index=True)
    df.to_csv(output)





def return_numerated_to_names(fasta_file, output_dir="/sternadi/home/volume3/taliakustin/phyVirus_final/with_names"):
    base = fasta_file.split("/")[-1].split(".fasta")[0].split(".with_names")[0]
    ids = {}
    count = 1
    seqs = list(SeqIO.parse(fasta_file, "fasta"))
    for s in seqs:
        ids[count] = s.description
        count += 1
    print(base)
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_final/%s*fas" % base) + glob.glob(
        "/sternadi/home/volume3/taliakustin/phyVirus_final/%s*fasta" % base) \
            + glob.glob(os.path.dirname(fasta_file) + "/%s*fasta" % base) + glob.glob(os.path.dirname(fasta_file) + "/%s*fas" %base)
    print(files)

    for f in files:
        if f == fasta_file:
            continue
        print(f)
        outfile = output_dir + "/" + f.split("/")[-1]
        seqs = list(SeqIO.parse(f, "fasta"))
        for s in seqs:
            des = int(s.description)
            s.description = ids[des]
            s.id = ids[des]
        SeqIO.write(seqs, outfile, format="fasta")
    files = glob.glob("%s/%s*fas" % (output_dir, base))
    for f in files:
        pbs_runners.phyml_runner(f, phylip=False)


def give_names_to_files(files, output_dir="/sternadi/home/volume3/taliakustin/phyVirus_final/new_names/"):
    #files - only fasta files
    for f in files:
        with open(f, "r") as input_file:
            data = input_file.readlines()
            first = []
            second = []
            for l in data:
                if ">" in l:
                    fir = l.split("|")[0].split(">")[-1]
                    sec = l.split("|")[1]
                    if fir != "NA":
                        first.append(fir)
                    if sec != "NA":
                        second.append(sec)
            if first != []:
                name = Counter(first).most_common()[0][0]
            elif second != []:
                name = Counter(second).most_common()[0][0]
            else:
                continue
            name = name.replace("/", "_")
            name = name.upper()
            family = f.split("/")[-1].split("_")[0]
        counter = 1
        new_base = "%s_%s_%i.fasta" % (family, name, counter)
        while os.path.isfile("/sternadi/home/volume3/taliakustin/phyVirus_final/names/%s" % new_base) or os.path.isfile("/sternadi/home/volume3/taliakustin/phyVirus_final/new_names/%s" % new_base)\
                or os.path.isfile(output_dir + "/%s" % new_base):
            counter += 1
            new_base = "%s_%s_%i.fasta" % (family, name, counter)
        print(f, new_base)
        to_move = glob.glob(f.split("fasta")[0] + "*")
        if len(to_move) == 2 or len(to_move) == 5:
            for i in to_move:
                if "fasta" in i:
                    os.system("cp %s %s/%s" % (i, output_dir, new_base))
                else:
                    new_base_without_fasta = new_base.split("fasta")[0]
                    output = "aln.best" + i.split("aln.best")[-1]
                    os.system("cp %s %s/%s%s" % (
                    i, output_dir, new_base_without_fasta, output))
        else:
            print(f)



def run_analysis_for_new_files_in_phyVirus_analysis(file):
    #Prerequisites - fasta file, alignemnet file and tree file in directories:
    #/sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta
    #/sternadi/home/volume3/taliakustin/phyVirus_analysis/aln
    #/sternadi/home/volume3/taliakustin/phyVirus_analysis/tree
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    file = check_filename(file)
    basename = file.split("/")[-1].split("fasta")[0]
    basename_4jobs = basename[:-1]
    aln = glob.glob(phyVirus_analysis_dir + "aln/%s*fas" % basename)[0]
    tree = glob.glob(phyVirus_analysis_dir + "tree/%s*tree.txt" % basename)[0]
    ran_something = False
    #add dinuc_freqs info
    if not os.path.isfile(phyVirus_analysis_dir + "dinuc_freqs/%sdinuc_averaged_odds_ratio" % basename):
        print("gets dinucleotide data")
        seqFileAnalyzer.get_dinucleotide_odds_ratio(file, output_dir=phyVirus_analysis_dir + "dinuc_freqs/")
        ran_something = True

    #add fasta freq info
    if not os.path.isfile(phyVirus_analysis_dir + "nuc_freqs/%sfasta_base_counts_info.csv" % basename):
        print("gets nuc data")
        seqFileAnalyzer.analyze_nuc_frequencies_from_fasta(file, output_dir=phyVirus_analysis_dir + "nuc_freqs/")
        ran_something = True

    #add phy alignment file
    if glob.glob(phyVirus_analysis_dir + "aln/%s*phy" % basename) == []:
        print("make phy alignment")
        aln_phy = convert_fasta_to_phylip(aln)
        ran_something = True
    else:
        aln_phy =  glob.glob(phyVirus_analysis_dir + "aln/%s*phy" % basename)[0]

    #add rooted file
    rooted_tree = phyVirus_analysis_dir + "rooted_trees/%saln.best.phy_phyml_rooted_tree.txt" % basename
    if not os.path.isfile(rooted_tree):
        print("roots tree")
        rooted_tree = root_at_midpoint(tree, outfile=rooted_tree)
        ran_something = True

    #run codon_alignment
    if len(glob.glob(phyVirus_analysis_dir + "fasta/%s*codon_aln*" % basename)) == 1: #if the codon aln file is in the fasta dir
        old_codon_aln = glob.glob(phyVirus_analysis_dir + "fasta/%s*codon_aln*" % basename)[0]
        codon_aln = phyVirus_analysis_dir + "codon_aln/" +  old_codon_aln.split("/")[-1]
        os.system("mv %s %s" % (old_codon_aln, codon_aln))
        print("gets nucleotide data")
        seqFileAnalyzer.analyze_nuc_frequencies_and_wobble_freqs_from_aln(codon_aln, output_dir=phyVirus_analysis_dir + "nuc_freqs/")
        ran_something = True
    elif len(glob.glob(phyVirus_analysis_dir + "codon_aln/%s*codon_aln*" % basename)) == 1:
        codon_aln = glob.glob(phyVirus_analysis_dir + "codon_aln/%s*codon_aln*" % basename)[0]
        if len(glob.glob(phyVirus_analysis_dir + "nuc_freqs/%s*codon_aln*" % basename)) == 0:
            seqFileAnalyzer.analyze_nuc_frequencies_and_wobble_freqs_from_aln(codon_aln, output_dir=phyVirus_analysis_dir + "nuc_freqs/")
            ran_something = True
    else: #have to run codon_aln
        codon_aln_job = pbs_runners.prank_codon_runner(file, alias=basename_4jobs)
        old_codon_aln = phyVirus_analysis_dir + "fasta/%scodon_aln.best.fas" % basename
        codon_aln = phyVirus_analysis_dir + "codon_aln/" +  old_codon_aln.split("/")[-1]
        print("runs codon alignment and nucleotide data")
        pbs_runners.script_runner("mv %s %s" % (old_codon_aln, codon_aln), run_after_job=codon_aln_job, alias=basename_4jobs, cmdname="mv")
        pbs_runners.script_runner("/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_analyze_nuc_frequencies_and_wobble_freqs_from_aln.py "
                                     "-a %s -o %s" % (codon_aln, phyVirus_analysis_dir + "nuc_freqs"), run_after_job=codon_aln_job, alias=basename_4jobs, cmdname="analyze_nuc_freqs")
        ran_something = True

    gtr_output = phyVirus_analysis_dir + "/gtr/%smlb7" % basename
    """
    #baseml and dirSel
    if not os.path.isfile(phyVirus_analysis_dir + "/gtr/%smlb7" % basename) and not os.path.isfile(phyVirus_analysis_dir + "/unrest/%smlb8" % basename):
        print("runs PAML")
        baseml_job7, baseml_job8 = run_baseml_on_aln_files(aln_phy, t=tree, t_rooted=rooted_tree, job_name=basename_4jobs, output_dirs=True)

        #run dirSel
        print("runs dirSel")
        dirsel_alternative_job = pbs_runners.script_runner("/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
                                  "-a %s -t %s  -m alternative -g %s"
                                  % (aln, tree, gtr_output), run_after_job=baseml_job7, job_name=basename_4jobs, cmdname="dirSel")
        dirsel_null_job = pbs_runners.script_runner("/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
                                  "-a %s -t %s  -m fixed_beta -g %s"
                                  % (aln, tree, gtr_output), run_after_job=baseml_job7, job_name=basename_4jobs, cmdname="dirSel_beta0")

        #run dirSel analysis
        print("runs dirSel analysis")
        pbs_runners.script_runner("/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel_analysis.py -b %s" % basename,
                                  run_after_job=dirsel_alternative_job, job_name=basename_4jobs, cmdname="dirSel_analysis")
        ran_something = True
    elif not os.path.isfile(phyVirus_analysis_dir + "/dirSel/%sdirSel.log" % basename):
        # run dirSel
        print("runs dirSel")
        dirsel_alternative_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m alternative -g %s"
            % (aln, tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel")
        dirsel_null_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m fixed_beta -g %s"
            % (aln, tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel_beta0")

        # run dirSel analysis
        print("runs dirSel analysis")
        pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel_analysis.py -b %s" % basename,
            run_after_job=dirsel_alternative_job, job_name=basename_4jobs, cmdname="dirSel_analysis")
        ran_something = True

    elif not os.path.isfile(phyVirus_analysis_dir + "/dirSel_analysis/%sdirSel.aln" % basename):
        print("runs dirSel analysis")
        pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel_analysis.py -b %s" % basename,
            job_name=basename_4jobs, cmdname="dirSel_analysis")
        ran_something = True
    """
    """
    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_extInt/%sdirSel_extInt.log" % basename):
        #run dirSel internal-external branches.
        print("runs dirSel-extInt")
        dirsel_extInt_job = pbs_runners.script_runner(
             "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
             "-a %s -t %s  -m extInt -g %s -d %s"
             % (aln, tree, gtr_output, "/sternadi/home/volume1/taliakustin/software/directional_selection_extInt/programs/directionalSelection/directionalSelection"),
             job_name=basename_4jobs, cmdname="dirSel-extInt")
        print(dirsel_extInt_job)
        ran_something = True

    
    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_extInt_null/%sdirSel_extIntNull.log" % basename):
        # run dirSel internal-external branches null model (S=0, fixed=True).
        print("runs dirSel-extInt-null")
        dirsel_extInt_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m dirSel_extIntNull -g %s -d %s"
            % (aln, tree, gtr_output,
               "/sternadi/home/volume1/taliakustin/software/directional_selection_extInt/programs/directionalSelection/directionalSelection"),
            job_name=basename_4jobs, cmdname="dirSel-extInt-null")
        print(dirsel_extInt_job)
        ran_something = True
    """
    """
    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_A_internal/%sdirSel_A_internal.log" % basename):
        # run dirSel internal-external branches only with selection towards A.
        print("runs dirSel_A_internal")
        dirsel_A_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m dirSel_A_internal -g %s -d %s"
            % (aln, tree, gtr_output,
               "/sternadi/home/volume1/taliakustin/software/directional_selection_A_internal_selection/programs/directionalSelection/directionalSelection"),
            job_name=basename_4jobs, cmdname="dirSel-A")
        print(dirsel_A_job)
        ran_something = True

    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_A_internal_null/%sdirSel_A_internal_null.log" % basename):
        # run dirSel internal-external branches only with selection towards A.
        print("runs dirSel_A_internal_null")
        dirsel_A_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m dirSel_A_internal_null -g %s -d %s"
            % (aln, tree, gtr_output,
               "/sternadi/home/volume1/taliakustin/software/directional_selection_A_internal_selection/programs/directionalSelection/directionalSelection"),
            job_name=basename_4jobs, cmdname="dirSel-A")
        print(dirsel_A_job)
        ran_something = True

    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_A_external/%sdirSel_A_external.log" % basename):
        # run dirSel internal-external branches only with selection towards A.
        print("runs dirSel_A_external")
        dirsel_A_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m dirSel_A_external -g %s -d %s"
            % (aln, tree, gtr_output,
               "/sternadi/home/volume1/taliakustin/software/directional_selection_A_external_selection/programs/directionalSelection/directionalSelection"),
            job_name=basename_4jobs, cmdname="dirSel-A")
        print(dirsel_A_job)
        ran_something = True

    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_A_external/%sdirSel_A_external_null.log" % basename):
        # run dirSel internal-external branches only with selection towards A.
        print("runs dirSel_A_external_null")
        dirsel_A_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m dirSel_A_external_null -g %s -d %s"
            % (aln, tree, gtr_output,
               "/sternadi/home/volume1/taliakustin/software/directional_selection_A_external_selection/programs/directionalSelection/directionalSelection"),
            job_name=basename_4jobs, cmdname="dirSel-A")
        print(dirsel_A_job)
        ran_something = True

    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_A_external_null/%sdirSel_A_internal_fixed_beta.log" % basename):
        # run dirSel internal-external branches only with selection towards A.
        print("runs dirSel_A_internal_fixed_beta")
        dirsel_A_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m dirSel_A_internal_fixed_beta -g %s -d %s"
            % (aln, tree, gtr_output,
               "/sternadi/home/volume1/taliakustin/software/directional_selection_A_internal_selection/programs/directionalSelection/directionalSelection"),
            job_name=basename_4jobs, cmdname="dirSel-A")
        print(dirsel_A_job)
        ran_something = True

    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_A_external_null/%sdirSel_A_external_fixed_beta.log" % basename):
        # run dirSel internal-external branches only with selection towards A.
        print("runs dirSel_A_external_fixed_beta")
        dirsel_A_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m dirSel_A_external_fixed_beta -g %s -d %s"
            % (aln, tree, gtr_output,
               "/sternadi/home/volume1/taliakustin/software/directional_selection_A_external_selection/programs/directionalSelection/directionalSelection"),
            job_name=basename_4jobs, cmdname="dirSel-A")
        print(dirsel_A_job)
        ran_something = True

    if not ran_something:
        print("didn't run any analysis on file %s" % basename)
    """


def run_dirSel_phyVirus_analysis(file):
    #Prerequisites - fasta file, alignemnet file and tree file in directories:
    #/sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta
    #/sternadi/home/volume3/taliakustin/phyVirus_analysis/aln
    #/sternadi/home/volume3/taliakustin/phyVirus_analysis/tree
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    file = check_filename(file)
    basename = file.split("/")[-1].split("fasta")[0]
    basename_4jobs = basename[:-1]
    aln = glob.glob(phyVirus_analysis_dir + "aln/%s*fas" % basename)[0]
    tree = glob.glob(phyVirus_analysis_dir + "tree/%s*tree.txt" % basename)[0]
    rooted_tree = glob.glob(phyVirus_analysis_dir + "rooted_trees/%s*rooted_tree.txt" % basename)[0]
    gtr_output = phyVirus_analysis_dir + "/gtr/%smlb7" % basename


    if not os.path.isfile(phyVirus_analysis_dir + "/dirSel_midpoint/%sdirSel.log" % basename):
        # run dirSel
        print("runs dirSel midpoint A internal")
        dirsel_alternative_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m midpoint -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_internal/programs/directionalSelection/directionalSelection"
            % (aln, rooted_tree, gtr_output), alias=basename_4jobs, cmdname="dirSel_selection_toA_internal")
        #dirsel_null_job = pbs_runners.script_runner(
        #    "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
        #    "-a %s -t %s  -m midpoint_null -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_internal/programs/directionalSelection/directionalSelection"
        #    % (aln, rooted_tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel_selection_toA_internal_null")
        """
        # run dirSel analysis
        print("runs dirSel analysis")
        pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel_analysis.py -b %s" % basename,
            run_after_job=dirsel_alternative_job, job_name=basename_4jobs, cmdname="dirSel_analysis")
        ran_something = True
        """


def run_dirSel_phyVirus_analysis_HIV(file):
    # Prerequisites - fasta file, alignemnet file and tree file in directories:
    # /sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta
    # /sternadi/home/volume3/taliakustin/phyVirus_analysis/aln
    # /sternadi/home/volume3/taliakustin/phyVirus_analysis/tree
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    file = check_filename(file)
    basename = file.split("/")[-1].split("fasta")[0]
    basename_4jobs = basename[:-1]
    aln = glob.glob(phyVirus_analysis_dir + "aln/%s*fas" % basename)[0]
    rooted_tree = glob.glob(phyVirus_analysis_dir + "rooted_trees/%s*rooted_tree.txt" % basename)[0]
    gtr_output = phyVirus_analysis_dir + "/gtr/%smlb7" % basename

    print(basename)
    #run midpoint rooting regular dirSel
    """
    dirsel_alternative_job = pbs_runners.script_runner(
        "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
        "-a %s -t %s  -m midpoint -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root/programs/directionalSelection/directionalSelection"
        % (aln, rooted_tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel_selection_toA_internal")
    dirsel_null_job = pbs_runners.script_runner(
        "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
        "-a %s -t %s  -m midpoint_null -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root/programs/directionalSelection/directionalSelection"
        % (aln, rooted_tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel_selection_toA_internal_null")
    """
    #run midpoint rooting - selection to A internal
    """
    dirsel_alternative_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m toA_internal -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_internal/programs/directionalSelection/directionalSelection"
            % (aln, rooted_tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel_selection_toA_internal")
    """
    dirsel_null_job = pbs_runners.script_runner(
        "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
        "-a %s -t %s  -m toA_internal_null -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_internal/programs/directionalSelection/directionalSelection"
        % (aln, rooted_tree, gtr_output), alias=basename_4jobs, cmdname="dirSel_selection_toA_internal_null")
    #run midpoint rooting - selection to A external
    """
    dirsel_alternative_job = pbs_runners.script_runner(
            "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
            "-a %s -t %s  -m toA_external -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_external//programs/directionalSelection/directionalSelection"
            % (aln, rooted_tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel_selection_toA_external")
    dirsel_null_job = pbs_runners.script_runner(
        "/sternadi/home/volume1/taliakustin/SternLab/phyVirus/run_dirSel.py "
        "-a %s -t %s  -m toA_external_null -g %s -d /sternadi/home/volume1/taliakustin/software/directional_selection_no_root_selection_toA_external//programs/directionalSelection/directionalSelection"
        % (aln, rooted_tree, gtr_output), job_name=basename_4jobs, cmdname="dirSel_selection_toA_external_null")
    """

def cp_files_from_new_names_to_phyVirus_analysis(file):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    file = check_filename(file)
    if seqFileTools.count_seqs_in_fasta(file) < 10:
        return
    basename = file.split("/")[-1].split("fasta")[0]
    dir = os.path.dirname(file)
    files = glob.glob(dir + "/%s*" % basename)

    for f in files:
        if "mlb" in f:
            if not is_complete_mlb_file(f):
                print("PAML is still running of this file: %s" % basename)
                print(f)
                #return

    new_file = phyVirus_analysis_dir + "fasta/" + os.path.basename(file)
    if os.path.isfile(new_file):
        print("this file was already moved: %s" % file)
        return

    print("Moving files: %s" % basename)

    for f in files:
        if ".numerate_index.txt" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "numerate_index/"))
        elif "rooted" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "rooted_trees/"))
        elif ".aln.best.fas" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "aln/"))
        elif ".fasta" in f:
            os.system("cp %s %s" % (f,  phyVirus_analysis_dir + "fasta/"))
        elif ".aln.best.phy_phyml_tree.txt" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "tree/"))
        elif ".aln.best.phy_phyml_stats.txt" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "tree/"))
        elif "aln.best.phy" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "aln/"))
        elif "codon_aln.best.fas" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "codon_aln/"))
        elif ".aln.best.phy_phyml_tree.txt.rooted.tree" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "rooted_trees/"))
        elif "mlb7" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "gtr/"))
        elif "ctl7" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "gtr/"))
        elif "mlb8" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "unrest/"))
        elif "ctl8" in f:
            os.system("cp %s %s" % (f, phyVirus_analysis_dir + "unrest/"))
        else:
            print("unknown file: %s" % f)





def get_alignment_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = fasta.split("/")[-1].split("fasta")[0]
    aln = glob.glob(phyVirus_analysis_dir + "aln/%s*fas" % basename)[0]
    return aln


def get_alignment_from_basename(basename):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    aln = glob.glob(phyVirus_analysis_dir + "aln/%s*fas" % basename)[0]
    return aln

def get_alignment_phy_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = fasta.split("/")[-1].split("fasta")[0]
    aln = glob.glob(phyVirus_analysis_dir + "aln/%s*phy" % basename)[0]
    return aln

def get_alignment_phy_from_basename(basename):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    aln = glob.glob(phyVirus_analysis_dir + "aln/%s*phy" % basename)[0]
    return aln

def get_codon_alignment_phy_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = fasta.split("/")[-1].split("fasta")[0]
    aln = glob.glob(phyVirus_analysis_dir + "codon_aln/%s*phy" % basename)[0]
    return aln

def get_codon_alignment_phy_from_basename(basename):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    aln = glob.glob(phyVirus_analysis_dir + "codon_aln/%s*phy" % basename)[0]
    return aln

def get_codon_alignment_fas_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = fasta.split("/")[-1].split("fasta")[0]
    aln = glob.glob(phyVirus_analysis_dir + "codon_aln/%s*fas" % basename)
    if aln == []:
        return aln
    return aln[0]

def get_rooted_tree_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = fasta.split("/")[-1].split("fasta")[0]
    rooted_tree = glob.glob(phyVirus_analysis_dir + "rooted_trees/%s*rooted_tree.txt" % basename)[0]
    return rooted_tree

def get_rooted_tree_from_basename(basename):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    rooted_tree = glob.glob(phyVirus_analysis_dir + "rooted_trees/%s*rooted_tree.txt" % basename)[0]
    return rooted_tree


def get_tree_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = fasta.split("/")[-1].split("fasta")[0]
    tree = glob.glob(phyVirus_analysis_dir + "tree/%s*tree.txt" % basename)[0]
    return tree


def get_labled_tree_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = fasta.split("/")[-1].split("fasta")[0]
    tree = glob.glob(phyVirus_analysis_dir + "rtl/%s*tree.txt" % basename)[0]
    return tree


def get_basename(fasta):
    return fasta.split("/")[-1].split("fasta")[0]


def run_fastml_rooted_tree(files):
    outdir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/fastml"
    for f in files:
        aln = get_alignment_from_fasta(f)
        rooted_tree = get_rooted_tree_from_fasta(f)
        pbs_runners.fastml_runner(aln, rooted_tree, outdir=outdir, fastml_path="/sternadi/home/volume1/taliakustin/software/phylogenyCode/programs/fastml/fastml", optBranchLen=False)


def run_fastml_codon_rooted_tree(files):
    outdir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/fastml_codon_midpoint_root"
    for f in files:
        if "Herpes" in f or "Pox" in f or "check" in f or "Reo" in f:
            continue
        aln = get_codon_alignment_fas_from_fasta(f)
        if aln == []:
            continue
        rooted_tree = get_rooted_tree_from_fasta(f)
        pbs_runners.fastml_runner(aln, rooted_tree, outdir=outdir, fastml_path="/sternadi/home/volume1/taliakustin/software/phylogenyCode/programs/fastml/fastml", optBranchLen=False, alias="fastml_codon")


def get_gtr_from_fasta(fasta):
    phyVirus_analysis_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/"
    fasta = check_filename(fasta)
    basename = get_basename(fasta)
    gtr_output = glob.glob(phyVirus_analysis_dir + "/gtr/%smlb7" % basename)[0]
    return gtr_output




def run_fastml(files):
    outdir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/fastml_unrooted_tree"
    for f in files:
        aln = get_alignment_from_fasta(f)
        tree = get_tree_from_fasta(f)
        pbs_runners.fastml_runner(aln, tree, outdir=outdir, fastml_path="/sternadi/home/volume1/taliakustin/software/phylogenyCode/programs/fastml/fastml", optBranchLen=False, alias="fastml_unrooted")




def remove_human_seqs_from_SIV_aln(file):
    file = check_filename(file)
    output = file.split(".fasta")[0] + "_no_human.fasta"
    seqs = list(SeqIO.parse(file, "fasta"))
    new_seqs = []
    for s in seqs:
        if s.id[0] != "H":
            new_seqs.append(s)
    SeqIO.write(new_seqs, output, "fasta")



def write_codeml_files():
    files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/fasta/*fasta")
    outdir1 = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/c1/"
    outdir2 = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/c2/"
    for f in files:
        if "Herpes" in f or "Pox" in f or "Reo" in f or "check" in f:
            continue
        base = get_basename(f)
        try:
            aln = get_codon_alignment_phy_from_fasta(f)
        except:
            print("no codon aln for {}".format(base))
            continue
        tree = get_labled_tree_from_fasta(f)
        ctl1 = "{}{}ctl".format(outdir1, base)
        ctl2 = "{}{}ctl".format(outdir2, base)
        mlc1 = "{}{}mlc".format(outdir1, base)
        mlc2 = "{}{}mlc".format(outdir2, base)

        write_ctl_codeml_deleterious_load_file(ctl1, aln, tree, mlc1, 0)
        write_ctl_codeml_deleterious_load_file(ctl2, aln, tree, mlc2, 2)