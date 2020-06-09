#! /usr/local/python_anaconda/bin/python3.4

import pandas as pd
import re
import math
import os
from file_utilities import check_filename



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

