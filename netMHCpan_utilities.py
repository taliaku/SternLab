#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
from file_utilities import check_filename
from Bio.Seq import Seq
from scipy.stats import ttest_ind
import numpy as np
import os
from functools import reduce
from optparse import OptionParser
from phyVirus.phyVirus_utilities import *
from Bio import SeqIO
import glob
from general_utilities import merge_two_dicts

def netMHCpan_to_csv(res_file, to_return=False):
    print(res_file)
    res_file = check_filename(res_file)
    output = res_file + ".csv"
    base = res_file.split("/")[-1].split(".")[0]
    if os.path.isfile(output):
        return
    df = pd.DataFrame()
    res = open(res_file, "r").read()
    if res == "/sternadi/home/volume1/taliakustin/software/netMHCpan-4.0/Linux_x86_64\n":
        print("Empty res file")
        return pd.DataFrame()
    if "cannot be found in allelenames list" in res:
        return pd.DataFrame()
    res = res.split("-----------------------------------------------------------------------------------")[2]
    res = res.split("\n")[1:-1]
    for r in res:
        r = r.strip()
        r = r.split(" ")
        new_r = []
        for i in r:
            if i!="":
                new_r.append(i)
        pos = new_r[0]
        allele = new_r[1]
        basename = new_r[10]
        rank = new_r[12]
        bind = new_r[11]
        peptide = new_r[2]
        core = new_r[3]
        icore = new_r[9]
        df = df.append({"pos":pos, "base":base, "basename":basename, "allele":allele, "rank":rank, "bind":bind, "peptide":peptide, "core":core, "icore":icore}, ignore_index = True)
    if to_return== True:
        return(df)
    df.to_csv(output)


def make_peptide_csv_for_each_seq(file):
    base = get_basename(file).split(".from_codon_aln_")[0].split("-translated")[0]
    output = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/{}_peptide_info.csv".format(base)
    res = pd.DataFrame()
    seqs = list(SeqIO.parse(file, "fasta"))
    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)

    for seq in seqs:
        print(seq.id)
        for i in range(0, len(seq), 3):
            pep_nuc = seq.seq[i:i+(9*3)]
            pep_aa = pep_nuc.translate()
            if len(pep_aa) != 9:
                continue
            res = res.append({"pep_nuc":str(pep_nuc), "pep_aa":(str(pep_aa)), "seq_name":seq.id, "base":base, "family":family, "baltimore":baltimore}, ignore_index=True)

    res.to_csv(output)


def go_over_netMHCpan_results(fasta_file):
    #from_codon_aln_fasta fasta file
    seqs = list(SeqIO.parse(fasta_file, "fasta"))
    basename = get_basename(fasta_file).split(".from_codon_aln_")[0]
    family = basename.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)
    output_freqs = fasta_file.split(".from_codon_aln_fasta")[0] + ".netMHCpan_freqs.csv"
    output_counts = fasta_file.split(".from_codon_aln_fasta")[0] + ".netMHCpan_counts.csv"
    if os.path.isfile(output_freqs):
        print("file exists")
        return
    #res_files = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/%s.*.netMHCpan_results.csv" % basename)
    #if len(seqs) != len(res_files):
    #    raise ("PROBLEM with file")
    nucs = ["A", "C", "G", "T"]
    nuc_counts = {"strong":{"A":0, "C":0, "G":0, "T":0}, "weak":{"A":0, "C":0, "G":0, "T":0}, "none":{"A":0, "C":0, "G":0, "T":0}}
    for seq_num in range(len(seqs)):
        res_file = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/%s.%i.netMHCpan_results.csv" % (basename, seq_num))[0]
        seq = seqs[seq_num]
        for n in nucs:
            nuc_counts["none"][n] += seq.seq.count(n)
        res = pd.read_csv(res_file)
        if res.empty:
            continue
        for i in range(len(seq)):
            for j in [9]: #do only peptides in the length of 9 (can be done for 10 and 11 as well)
                    pep_nuc = seq[i:i+(j*3)]
                    pep_aa = pep_nuc.seq.translate()
                    if not res.loc[res["peptide"] == str(pep_aa)].empty:
                        if (len(res.loc[res["peptide"] == str(pep_aa)])) > 1:
                            print(res.loc[res["peptide"] == str(pep_aa)])

                            bind = float(list(res.loc[res["peptide"] == str(pep_aa)].bind)[0])
                            print(bind)
                        else:
                            bind = float(res.loc[res["peptide"] == str(pep_aa)].bind)
                        if bind > 0.5:
                            property = "strong"
                        else:
                            property = "weak"
                        for n in nucs:
                            nuc_counts["none"][n] -= pep_nuc.seq.count(n)
                            nuc_counts[property][n] += pep_nuc.seq.count(n)
    freqs = {}
    counts_output = pd.DataFrame()
    freqs_output = pd.DataFrame()
    for t in nuc_counts:
        freqs[t] = {}
        for n in nucs:
            freqs[t][n] = float(nuc_counts[t][n]) / sum(nuc_counts[t].values())
            counts_output = counts_output.append({"basename":basename, "family":family, "baltimore":baltimore,
                                    "type":t, "nuc":n, "count":nuc_counts[t][n]}, ignore_index=True)
            freqs_output = freqs_output.append({"basename": basename, "family": family, "baltimore": baltimore,
                                    "type": t, "nuc": n, "count": freqs[t][n]}, ignore_index=True)

    freqs_output.to_csv(output_freqs)
    counts_output.to_csv(output_counts)



def count_amico_acids_in_netMHCpan_results(fasta_files):
    df = pd.DataFrame()
    for f in fasta_files:
        basename = get_basename(f)
        family = basename.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        res_files = glob.glob(
        "/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/%s*.netMHCpan_results.csv" % (
        basename))
        output_dic = {"basename":basename, "family":family, "baltimore":baltimore}
        codons = {}
        for r in res_files:
            res = pd.read_csv(r)
            for pep in res["peptide"]:
                if len(pep) == 9:
                    for p in pep:
                        if p in codons:
                            codons[p] += 1
                        else:
                            codons[p] = 1
            full_output = merge_two_dicts(output_dic, codons)
        df = df.append(full_output, ignore_index=True)
    print(df)


def merge_peptide_file_with_netMHCpan4_results(net_file):
    output = net_file.split("_results.csv")[0] + "_freqs.csv"
    if os.path.isfile(output):
        print("output file exists")
        return
    net = pd.read_csv(net_file)
    if net.empty:
        print("empty")
        return
    all_pep = pd.read_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/peptide_info_all.csv")
    uniq_net = net.groupby(["bind", "peptide", "rank", "allele"]).size().reset_index().rename(columns={0:'count'})
    uniq_net["property"] = "weak"
    uniq_net.loc[uniq_net.bind > 0.5, "property"] = "strong"

    merged = pd.merge(uniq_net, all_pep, left_on=["peptide"], right_on=["pep_aa"], how="inner")
    merged["A"] = merged["pep_nuc"].apply(lambda x: x.count("A"))
    merged["C"] = merged["pep_nuc"].apply(lambda x: x.count("C"))
    merged["G"] = merged["pep_nuc"].apply(lambda x: x.count("G"))
    merged["U"] = merged["pep_nuc"].apply(lambda x: x.count("T"))
    summed = merged.groupby(["family", "property", "allele"])["A", "C", "G", "U"].apply(lambda x: x.astype(int).sum()).reset_index()
    summed["sum"] = summed["A"] + summed["C"] + summed["G"] + summed["U"]
    summed["count_A"] = summed["A"]
    summed["count_C"] = summed["C"]
    summed["count_G"] = summed["G"]
    summed["count_U"] = summed["U"]
    summed["A"] = summed["count_A"] / summed["sum"]
    summed["C"] = summed["count_C"] / summed["sum"]
    summed["G"] = summed["count_G"] / summed["sum"]
    summed["U"] = summed["count_U"] / summed["sum"]

    summed.to_csv(output)
    print(f"saved summary to {output}")


def merge_peptide_file_with_netMHCpan4_results_get_null_results(net_file):
    output = net_file.split("_results.csv")[0] + "_freqs_not_MHC.csv"
    if os.path.isfile(output):
        print("output file exists")
        return
    net = pd.read_csv(net_file)
    if net.empty:
        print("empty")
        return
    all_pep = pd.read_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/peptide_info_all.csv")
    uniq_net = net.groupby(["bind", "peptide", "rank", "allele"]).size().reset_index().rename(columns={0:'count'})
    uniq_net["property"] = "weak"
    uniq_net.loc[uniq_net.bind > 0.5, "property"] = "strong"

    merged = pd.merge(uniq_net, all_pep, left_on=["peptide"], right_on=["pep_aa"], how="right")
    null = merged[merged.isnull().any(axis=1)]


    null["A"] = null["pep_nuc"].apply(lambda x: x.count("A"))
    null["C"] = null["pep_nuc"].apply(lambda x: x.count("C"))
    null["G"] = null["pep_nuc"].apply(lambda x: x.count("G"))
    null["U"] = null["pep_nuc"].apply(lambda x: x.count("T"))
    summed = null.groupby(["family", "property", "allele"])["A", "C", "G", "U"].apply(lambda x: x.astype(int).sum()).reset_index()
    summed["sum"] = summed["A"] + summed["C"] + summed["G"] + summed["U"]
    summed["count_A"] = summed["A"]
    summed["count_C"] = summed["C"]
    summed["count_G"] = summed["G"]
    summed["count_U"] = summed["U"]
    summed["A"] = summed["count_A"] / summed["sum"]
    summed["C"] = summed["count_C"] / summed["sum"]
    summed["G"] = summed["count_G"] / summed["sum"]
    summed["U"] = summed["count_U"] / summed["sum"]

    summed.to_csv(output)
    print(f"saved summary to {output}")


def go_over_MHC_results_by_netMHC_file(net_file):
    nucs_list = ["A", "C", "G", "T"]
    net = pd.read_csv(net_file)
    host_info = pd.read_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/host_info_per_seq.csv")
    protein_info = pd.read_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/phyVirus_database_protein_info.csv")
    if net.empty:
        print("empty")
        return
    output = net_file.split("_results.csv")[0] + "_freqs_by_allele.csv"
    if os.path.isfile(output):
        print("exists")
        return
    allele = net.allele[0]
    results = pd.DataFrame()
    fastas = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/*.from_codon_aln_fasta")
    for f in fastas:
        if "Reo" in f or "Herpes" in f or "Pox" in f or "Nairo" in f or "Nairo" in f:
            continue
        base = f.split("/")[-1].split(".from_codon_aln_fasta")[0]
        family = base.split("_")[0]
        baltimore = get_baltimore_classifiaction(family)
        seqs  = list(SeqIO.parse(f, "fasta"))
        protein_info.loc[protein_info.base == base, "protein_type"]
        print(base, family, len(seqs))
        for s in seqs:
            print(s.id)
            length = len(s.seq)
            seq_host = host_info.loc[(host_info.seq == s.id) & (host_info.family == family), "host_class"]
            i = 0
            nucleotide_counts = {"None":{"A":0, "C":0, "G":0, "T":0}, "Weak":{"A":0, "C":0, "G":0, "T":0}, "Strong":{"A":0, "C":0, "G":0, "T":0}}
            while i != length - 27:
                nucs = s.seq[i:i+27]
                pep = str(nucs.translate())
                pep_res = net.loc[(net.peptide==pep) & (net.basename == s.id)]
                if pep_res.empty:
                    i = i+3
                    to_count = s.seq[i:i+3]
                    for n in nucs_list:
                        nucleotide_counts["None"][n] += to_count.count(n)
                else:
                    to_add = 27
                    if float(pep_res.bind) > 0.5:
                        property = "Strong"
                    else:
                        property = "Weak"
                    for n in nucs_list:
                        nucleotide_counts[property][n] += nucs.count(n)
                    for j in range(0,27, 3):
                        nucs = s.seq[i+j:i + j + 27]
                        pep = str(nucs.translate())
                        pep_res = net.loc[(net.peptide == pep) & (net.base == base) & (net.basename == s.id)]
                        if not pep_res.empty:
                            to_add = j + 27
                    i = i + to_add


            for p in nucleotide_counts.keys():
                nucs_to_add = nucleotide_counts[p]
                nucs_to_add["base"] = base
                nucs_to_add["baltimore"] = baltimore
                nucs_to_add["seq_host"] = seq_host
                nucs_to_add["family"] = family
                nucs_to_add["seq_name"] = s.id
                results = results.append(nucs_to_add, ignore_index=True)
            print(results)

    results["allele"] = allele
    results.to_csv(output)

def go_over_MHC_results_by_fasta_file(f):


    nucs_list = ["A", "C", "G", "T"]
    host_info = pd.read_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/host_info_per_seq.csv")
    protein_info = pd.read_csv(
        "/sternadi/home/volume3/taliakustin/phyVirus_analysis/phyVirus_database_protein_info.csv")

    base = f.split("/")[-1].split(".from_codon_aln_fasta")[0]

    output = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/superfamily/{base}_freqs.csv"
    """
    if os.path.isfile(output):
        output_text = open(output, "r").read()
        if output_text != '""\n':
            print("output exists")
            return
    else:
        print("still running!!!")
        return
    """

    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)
    seqs = list(SeqIO.parse(f, "fasta"))
    protein_type = protein_info.loc[protein_info.base == base, "protein_type"].max()
    net_mhc_results_file = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/netMHC_no_dups/{base}_netMHC_results.csv"
    net_mhc_results_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/superfamily/superfamily_netMHC_prediction.csv"
    net_mhc_results = pd.read_csv(net_mhc_results_file)
    net_mhc_results = net_mhc_results.loc[net_mhc_results.base != "check"]
    #net_mhc_results["base"] = net_mhc_results["basename"].apply(lambda x: "_".join(x.split("_")[:-1]))
    results = pd.DataFrame()


    for s in seqs:
        length = len(s.seq)
        seq_host = host_info.loc[(host_info.seq == s.id) & (host_info.family == family), "host_class"].max()
        seq_id = s.id
        if "8" not in seq_id:
            continue
        seq_id_for_search = seq_id[:15]
        if "." in seq_id_for_search:
            seq_id_for_search = seq_id_for_search.replace(".", "_")
        net_mhc_for_seq = net_mhc_results.loc[net_mhc_results.basename == seq_id_for_search]

        alleles_for_seq = list(np.unique(net_mhc_for_seq.allele))
        for a in alleles_for_seq:
            print(a)
            nucleotide_counts = {"None": {"A": 0, "C": 0, "G": 0, "T": 0}, "Weak": {"A": 0, "C": 0, "G": 0, "T": 0},
                                 "Strong": {"A": 0, "C": 0, "G": 0, "T": 0}}
            i = 0
            while i < length - 26:
                #print(i)
                nucs = s.seq[i:i + 27]
                if (len(nucs) != 27):
                    print(nucs, i)
                pep = str(nucs.translate())
                pep_res = net_mhc_for_seq.loc[(net_mhc_for_seq.peptide==pep) & (net_mhc_for_seq.allele == a) & (net_mhc_for_seq.basename == s.id)]

                if pep_res.empty:
                    i = i+3
                    to_count = s.seq[i:i+3]
                    for n in nucs_list:
                        nucleotide_counts["None"][n] += to_count.count(n)
                else:
                    print(pep_res)
                    pep_res = pep_res.iloc[0]
                    to_add = 27
                    if float(pep_res.bind) > 0.5:
                        property = "Strong"
                    else:
                        property = "Weak"

                    j = 3
                    nucs = nucs

                    while(pep_res.empty == False):
                        if i+j+27 > length-27:
                            continue
                        next_nucs = s.seq[i + j:i + j + 27]
                        if (len(nucs) != 27):
                            print(nucs, i, j)
                        pep = str(next_nucs.translate())


                    for j in range(3,27, 3):
                        if i+j+27 > length-27:
                            continue
                        next_nucs = s.seq[i+j:i + j + 27]
                        if (len(nucs) != 27):
                            print(nucs, i, j)
                        pep = str(next_nucs.translate())
                        pep_res = net_mhc_for_seq.loc[(net_mhc_for_seq.peptide == pep) & (net_mhc_for_seq.allele == a) & (net_mhc_for_seq.basename == s.id)]
                        if not pep_res.empty:
                            for n in nucs_list:
                                nucleotide_counts[property][n] += nucs[:j].count(n)
                            pep_res_new = pep_res.iloc[0]
                            print("overlap {}".format(j))
                            print(pep_res_new)
                            if float(pep_res_new.bind) > 0.5:
                                property = "Strong"
                            else:
                                property = "Weak"
                            for n in nucs_list:
                                nucleotide_counts[property][n] += next_nucs.count(n)
                            to_add = j + 27
                        else: # empty
                            for n in nucs_list:
                                nucleotide_counts[property][n] += nucs.count(n)
                    i = i + to_add


            for p in nucleotide_counts.keys():
                nucs_to_add = nucleotide_counts[p]
                nucs_to_add["base"] = base
                nucs_to_add["baltimore"] = baltimore
                nucs_to_add["seq_host"] = seq_host
                nucs_to_add["family"] = family
                nucs_to_add["seq_name"] = seq_id
                nucs_to_add["allele"] = a
                nucs_to_add["property"] = p
                nucs_to_add["protein_type"] = protein_type

                results = results.append(nucs_to_add, ignore_index=True)
    results.to_csv(output)

def go_over_MHC_results_by_fasta_with_set(f):
    nucs_list = ["A", "C", "G", "T"]
    aas_list = ['F', 'L', 'I', 'M', 'V', 'S', 'P', 'T', 'A', 'Y', 'H', 'Q', 'N', 'K', 'D', 'E', 'C', 'W', 'R', 'G']
    host_info = pd.read_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/host_info_per_seq.csv")
    protein_info = pd.read_csv(
        "/sternadi/home/volume3/taliakustin/phyVirus_analysis/phyVirus_database_protein_info.csv")
    base = f.split("/")[-1].split(".from_codon_aln_fasta")[0].split(".all_alleles.csv")[0].split("-translated.fasta")[0]
    output_nucs = f"/sternadi/nobackup/volume1/talia_temp/netMHC_by_seq_250/{base}_nuc_freqs.csv"
    output_aas = f"/sternadi/nobackup/volume1/talia_temp/netMHC_by_seq_250/{base}_aas_freqs.csv"

    seq_file = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/{base}.from_codon_aln_fasta"

    if os.path.isfile(output_nucs):
        print("file exists")
        return
    family = base.split("_")[0]
    baltimore = get_baltimore_classifiaction(family)
    seqs = list(SeqIO.parse(seq_file, "fasta"))
    protein_type = protein_info.loc[protein_info.base == base, "protein_type"].max()
    #net_mhc_results_file = f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/netMHC_no_dups/{base}_netMHC_results.csv"
    net_mhc_results_file = f"/sternadi/nobackup/volume1/talia_temp/netMHC_by_seq_250/{base}.all_alleles.csv"
    #net_mhc_results_file = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/superfamily/superfamily_netMHC_prediction.csv"
    net_mhc_results = pd.read_csv(net_mhc_results_file)
    # net_mhc_results["base"] = net_mhc_results["basename"].apply(lambda x: "_".join(x.split("_")[:-1]))
    results_nucs = pd.DataFrame()
    results_aas = pd.DataFrame()

    for s in seqs:
        sequence_length = len(s.seq)
        seq_host = host_info.loc[(host_info.seq == s.id) & (host_info.family == family), "host_class"].max()
        seq_id = s.id
        seq_len_name = len(seq_id)
        seq_id_for_search = seq_id[:15]
        if "." in seq_id_for_search:
            seq_id_for_search = seq_id_for_search.replace(".", "_")
        #net_mhc_for_seq = net_mhc_results.loc[net_mhc_results.basename == seq_id_for_search]
        net_mhc_for_seq = net_mhc_results.loc[net_mhc_results.ID == seq_id_for_search]

        alleles_for_seq = list(np.unique(net_mhc_for_seq.allele))
        seq_translated = s.seq.translate()
        seq_translated_len = len(seq_translated)
        seq_nucs = {}
        seq_aas= {}
        for n in nucs_list:
            seq_nucs[n] = s.seq.count(n)
        for n in aas_list:
            seq_aas[n] = seq_translated.count(n)
        for a in alleles_for_seq:
            strong_nucs_set = set()
            weak_nucs_set = set()
            strong_aas_set = set()
            weak_aas_set = set()
            strong_nucs = {"A": 0, "C": 0, "G": 0, "T": 0}
            weak_nucs = {"A": 0, "C": 0, "G": 0, "T": 0}
            none_nucs = {"A": 0, "C": 0, "G": 0, "T": 0}
            strong_aas = {"F": 0, "L": 0, "I": 0, "M": 0, "V": 0, "S": 0, "P": 0, "T": 0, "A": 0, "Y": 0, "H": 0, "Q": 0, "N": 0,
               "K": 0, "D": 0, "E": 0, "C": 0, "W": 0, "R": 0, "S": 0, "G": 0}
            weak_aas = {"F": 0, "L": 0, "I": 0, "M": 0, "V": 0, "S": 0, "P": 0, "T": 0, "A": 0, "Y": 0, "H": 0, "Q": 0, "N": 0,
               "K": 0, "D": 0, "E": 0, "C": 0, "W": 0, "R": 0, "S": 0, "G": 0}
            none_aas = {"F": 0, "L": 0, "I": 0, "M": 0, "V": 0, "S": 0, "P": 0, "T": 0, "A": 0, "Y": 0, "H": 0, "Q": 0, "N": 0,
               "K": 0, "D": 0, "E": 0, "C": 0, "W": 0, "R": 0, "S": 0, "G": 0}
            temp_res_nucs = pd.DataFrame()
            temp_res_aas = pd.DataFrame()
            #net_mhc_for_allele = net_mhc_results.loc[(net_mhc_results.basename == seq_id_for_search) & (net_mhc_results.allele == a)]
            net_mhc_for_allele = net_mhc_results.loc[(net_mhc_results.ID == seq_id_for_search) & (net_mhc_results.allele == a)]
            for index, row in net_mhc_for_allele.iterrows():
                #pos = row.pos
                pos = row.Pos
                #bind = row.bind
                property = row.property
                nucs_poss_to_add = list(range((pos)*3, (pos+9)*3))
                aas_poss_to_add = list(range((pos), (pos+9)))

                if property == "strong":
                    for p in aas_poss_to_add:
                        if p >= seq_translated_len:
                            continue
                        strong_aas_set.add(p)
                elif property == "weak":
                    for p in aas_poss_to_add:
                        if p >= seq_translated_len:
                            continue
                        weak_aas_set.add(p)

                if seq_len_name > 15:
                    pep = row.Peptide
                    pep_in_seq = s.seq[nucs_poss_to_add[0]:nucs_poss_to_add[-1]+1].translate()
                    if pep == pep_in_seq:
                        if property == "strong":
                            for p in nucs_poss_to_add:
                                if p >= sequence_length:
                                    continue
                                strong_nucs_set.add(p)
                        elif property == "weak":
                            for p in nucs_poss_to_add:
                                if p >= sequence_length:
                                    continue
                                weak_nucs_set.add(p)
                else:
                    if property == "strong":
                        for p in nucs_poss_to_add:
                            if p >= sequence_length:
                                continue
                            strong_nucs_set.add(p)
                    elif property == "weak":
                        for p in nucs_poss_to_add:
                            if p >= sequence_length:
                                continue
                            weak_nucs_set.add(p)



            weak_nucs_set.difference_update(strong_nucs_set)
            weak_aas_set.difference_update(strong_aas_set)


            for char in strong_nucs_set:
                if s.seq[char] not in ["A", "C", "G", "T"]:
                    continue
                strong_nucs[s.seq[char]] += 1
            for char in weak_nucs_set:
                if s.seq[char] not in ["A", "C", "G", "T"]:
                    continue
                weak_nucs[s.seq[char]] += 1
            for n in nucs_list:
                none_nucs[n] = seq_nucs[n] - strong_nucs[n] - weak_nucs[n]
            for char in strong_aas_set:
                if seq_translated[char] not in aas_list:
                    continue
                strong_aas[seq_translated[char]] += 1
            for char in weak_aas_set:
                if seq_translated[char] not in aas_list:
                    continue
                weak_aas[seq_translated[char]] += 1
            for n in aas_list:
                none_aas[n] = seq_aas[n] - strong_aas[n] - weak_aas[n]



            strong_nucs["property"] = "strong"
            weak_nucs["property"] = "weak"
            none_nucs["property"] = "none"
            strong_aas["property"] = "strong"
            weak_aas["property"] = "weak"
            none_aas["property"] = "none"

            temp_res_nucs = temp_res_nucs.append(strong_nucs, ignore_index=True)
            temp_res_nucs = temp_res_nucs.append(weak_nucs, ignore_index=True)
            temp_res_nucs = temp_res_nucs.append(none_nucs, ignore_index=True)
            temp_res_nucs["base"] = base
            temp_res_nucs["baltimore"] = baltimore
            temp_res_nucs["seq_host"] = seq_host
            temp_res_nucs["family"] = family
            temp_res_nucs["seq_name"] = seq_id
            temp_res_nucs["allele"] = a
            temp_res_nucs["protein_type"] = protein_type
            results_nucs = results_nucs.append(temp_res_nucs, ignore_index=True)



            temp_res_aas = temp_res_aas.append(strong_aas, ignore_index=True)
            temp_res_aas = temp_res_aas.append(weak_aas, ignore_index=True)
            temp_res_aas = temp_res_aas.append(none_aas, ignore_index=True)
            temp_res_aas["base"] = base
            temp_res_aas["baltimore"] = baltimore
            temp_res_aas["seq_host"] = seq_host
            temp_res_aas["family"] = family
            temp_res_aas["seq_name"] = seq_id
            temp_res_aas["allele"] = a
            temp_res_aas["protein_type"] = protein_type
            results_aas = results_aas.append(temp_res_aas, ignore_index=True)

    results_nucs.to_csv(output_nucs)
    results_aas.to_csv(output_aas)




def main():
    go_over_MHC_results_by_fasta_with_set("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/Flavi_NS4B_1.from_codon_aln_fasta")

if __name__ == "__main__":
    main()

