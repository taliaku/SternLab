#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import re
import math
import os
from file_utilities import check_filename
from pbs_runners import script_runner
from Bio import SeqIO
import pandas as pd



def read_bowtie2_summary_paired_and_unpaired(file, results={}):
    file = check_filename(file)
    summary = open(file, "r").read()
    if not "overall alignment rate" in summary:
        print("no bowtie2 results")
        return {}
    results["file"] = summary.split("reads; of these:")[0].split("\n")[-2].split(">")[-1].strip()
    base = summary.split("reads; of these:")[0].split("\n")[-2].split(">")[-1].strip().split("/")[-1].split(".bacteria_output.bam")[0]
    results["base"] = base.split("_L00")[0]
    results["lane"] = base.split("L00")[1].split("_")[0]
    results["patient"] = results["base"].split("_")[0]
    results["genetic_material"] =  results["base"].split("_")[1]
    results["S"] = results["base"].split("_")[2]
    results["total_read_count"] = int(summary.split("reads; of these:")[0].split("\n")[-1].strip())
    info = summary.split("reads; of these:")[1].strip().split("\n")
    info = [x.strip() for x in info]
    results["pair_count"] = int(info[0].split(" ")[0])
    results["pair_percent"] = float(info[0].split("(")[1].split("%")[0])
    results["concordantly_0_times_count"] = int(info[1].split(" ")[0])
    results["concordantly_0_times_percent"] = float(info[1].split("(")[1].split("%")[0])

    results["concordantly_1_time_count"] = int(info[2].split(" ")[0])
    results["concordantly_1_time_percent"] = float(info[2].split("(")[1].split("%")[0])

    results["concordantly_more_than_1_time_count"] = int(info[3].split(" ")[0])
    results["concordantly_more_than_1_time_percent"] = float(info[3].split("(")[1].split("%")[0])

    results["discordantly_1_time_count"] = int(info[6].split(" ")[0])
    results["discordantly_1_time_percent"] = float(info[6].split("(")[1].split("%")[0])

    results["unpaired_mates_count"] = int(info[9].split(" ")[0])

    results["mates_0_count"] = int(info[10].split(" ")[0])
    results["mates_0_percent"] = float(info[10].split("(")[1].split("%")[0])

    results["mates_1_count"] = int(info[11].split(" ")[0])
    results["mates_1_percent"] = float(info[11].split("(")[1].split("%")[0])

    results["mates_more_than_1_count"] = int(info[12].split(" ")[0])
    results["mates_more_than_1_percent"] = float(info[12].split("(")[1].split("%")[0])

    results["overall_alignment_rate"] = float(info[13].split(" ")[0].split("%")[0])
    results["total_reads_mapped"] =  results["mates_more_than_1_count"] +  results["mates_1_count"] + 2*results["discordantly_1_time_count"] +2*results["concordantly_more_than_1_time_count"] + 2*results["concordantly_1_time_count"]


    return (results)


def read_bowtie2_summary(file, results={}):
    file = check_filename(file)
    summary = open(file, "r").read()
    if not "overall alignment rate" in summary:
        print("no bowtie2 results")
        return {}
    results["file"] = summary.split("reads; of these:")[0].split("\n")[-2].split(">")[-1].strip()
    base = summary.split("reads; of these:")[0].split("\n")[-2].split(">")[-1].strip().split("/")[-1].split(
        ".bacteria_output.bam")[0]
    results["base"] = base.split("_L00")[0]
    results["lane"] = base.split("L00")[1].split("_")[0]
    results["patient"] = results["base"].split("_")[0]
    results["genetic_material"] = results["base"].split("_")[1]
    results["total_read_count"] = int(summary.split("reads; of these:")[0].split("\n")[-1].strip())
    info = summary.split("reads; of these:")[1].strip().split("\n")
    info = [x.strip() for x in info]
    results["pair_count"] = int(info[0].split(" ")[0])
    results["pair_percent"] = float(info[0].split("(")[1].split("%")[0])
    results["concordantly_0_times_count"] = int(info[1].split(" ")[0])
    results["concordantly_0_times_percent"] = float(info[1].split("(")[1].split("%")[0])
    results["concordantly_1_time_count"] = int(info[2].split(" ")[0])
    results["concordantly_1_time_percent"] = float(info[2].split("(")[1].split("%")[0])
    results["concordantly_more_than_1_time_count"] = int(info[3].split(" ")[0])
    results["concordantly_more_than_1_time_percent"] = float(info[3].split("(")[1].split("%")[0])
    results["overall_alignment_rate"] = float(info[4].split(" ")[0].split("%")[0])
    results["reads_aligned"] = results["concordantly_1_time_count"] + results["concordantly_more_than_1_time_count"]
    return(results)


# TODO - consider chaning to conda samtools instead of using a module from server
def sam_2_bam(sam_file, bam_file=None):
    """
    convert sam file into bam file
    """
    if not bam_file:
        bam_file = sam_file.replace(".sam", ".bam")
    script_runner(f"module load samtools/samtools-1.6\nsamtools view -S -b {sam_file} > {bam_file}")
    return bam_file

def bam_2_sorted_bam(bam_file, sorted_file=None):
    """
    sort bam file
    """

    if not sorted_file:
        sorted_file = bam_file.replace(".bam", ".sorted.bam")
    script_runner(f"module load samtools/samtools-1.6\nsamtools sort {bam_file} -o {sorted_file}")
    return bam_file  

def bowtie_db_2_taxonomy(db_fasta_path, out=None):
    """
    map tax id to identifier in output file
    """
    res = []
    for record in SeqIO.parse(db_fasta_path, "fasta"):
        name = record.description.split('| ')[-1].split(',')[0]
        identifier = record.id
        res.append(tuple(name, identifier))
    df = pd.DataFrame(res, columns=["name","id"])
    if out:
        df.to_csv(out, index=False)

    return df

def summarize_bam(bam_file, out=None):
    """
    run samtools summary on an indexed version of a bam file
    :param bam_file: sorted bam file
    """
    if not out:
        out = "./summerized_bam"
    script_runner(f"module load samtools/samtools-1.6\nsamtools idxstats {bam_file} | awk '$4 != 0' | head -n-1 | sort -k3 -V > {out}")
    return bam_file  

def merge_bam_and_summary(summary, tax_mapping, out=None):
    """
    merged summary file from bam with following tax mapping from db used to map reads
    """
    bam_summary = pd.read_table(summary, names=["id", "seq_len","num_mapped","num_unmapped"])
    mapping = pd.read_csv(mapping)

    merged = bam_summary.merge(mapping, on="id")
    if out:
        merged.to_csv(out, index=False)
    return merged


        

