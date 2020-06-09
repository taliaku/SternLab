#! python/python-anaconda3.2019.7

import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import subprocess

PROBABILITY_LIMIT = 0.85

def FindFilesInDir(dir_path, file_type):
	file_path = dir_path + "/*" + file_type
	list_of_files = sorted(glob.glob(file_path))
	num_of_files = len(list_of_files)
	if num_of_files > 0:
		for file in list_of_files:
			size = os.path.getsize(file)
			if size == 0:
				raise Exception("Unexpected error, some of the " + file_type + " files in " + dir_path + " are empty\n")
       
	return list_of_files

def summarize_stats (dir_path, sample_basename_pattern, freqs_file_path, Coverage):
	os.chdir(dir_path)
	pipeline_summary = dir_path + "/Summary.txt"
	with open(pipeline_summary, "a") as o:
		file_type = ".fasta"
		input_fasta_files = FindFilesInDir(dir_path, file_type)
		if len(input_fasta_files) > 0:
			try:
				num_reads = subprocess.getoutput("grep '^>' *fasta -h | wc -l")
				o.write("Total number of reads: {}\n".format(int(num_reads)))
			except:
				print("\nWarning. Unable to summarize num of reads per sample " + sample_basename_pattern + ", or cannot open or write to file " + pipeline_summary + "\n")

		file_type = ".stats"
		input_stats_files = FindFilesInDir(dir_path, file_type)
		if len(input_stats_files) > 0:		
			try:
				num_contributing_reads = subprocess.getoutput("grep 'Number of reads contributing to frequency counts' *stats -h | awk -F '\t' '{sum+=$2}; END{print sum}'")
				num_contributing_bases = subprocess.getoutput("grep 'Number of bases contributing to frequency counts' *stats -h | awk -F '\t' '{sum+=$2}; END{print sum}'")
				num_non_contributing_reads = subprocess.getoutput("grep 'Number of non contributing reads' *stats -h | awk -F '\t' '{sum+=$2}; END{print sum}'")
				num_mapped_reads = int(num_contributing_reads) + int(num_non_contributing_reads)
				percentage_mapped = round((int(num_mapped_reads)/int(num_reads))*100,2)
				percentage_contribution = round((int(num_contributing_reads)/int(num_mapped_reads))*100,2)
				o.write("Number of reads mapped to reference: {}\n".format(int(num_mapped_reads)))
				o.write("% of mapped reads: {}\n".format(percentage_mapped))
				o.write("Total number of reads contributing to frequency count: {}\n".format(int(num_contributing_reads)))
				o.write("% of reads contributing to frequency count: {}\n".format(percentage_contribution))
				o.write("Total number of bases contributing to frequency count: {}\n\n".format(int(num_contributing_bases)))
				
				get_repeats = []
				get_repeats.append(subprocess.getoutput("grep -h 'repeats, ' *.stats | awk -F \" \" '{print $1}'"))
				repeats = get_repeats[0].split("\n")
				repeats = list(set(repeats))
				repeats.sort(key = int)
				for repeat in repeats:
					num_contributing_reads_by_repeat = subprocess.getoutput("grep '" + repeat + " repeats, ' *.stats | awk -F \" \" '{sum+=$3}; END {print sum}'")
					num_contributing_bases_by_repeat = subprocess.getoutput("grep '" + repeat + " repeats, ' *.stats | awk -F \" \" '{sum+=$7}; END {print sum}'")
					o.write("{} repeats contributing {} reads and {} bases to frequency count\n".format(repeat, int(num_contributing_reads_by_repeat), int(num_contributing_bases_by_repeat)))
				print(freqs_file_path)
				freqs_data = pd.read_csv(freqs_file_path, sep=',')
				mutation_data = freqs_data[(freqs_data["rank"] != 0) & (freqs_data["coverage"] > Coverage) & (freqs_data["probability"] >= PROBABILITY_LIMIT)]
				print(f"Sum of Mutations: {len(mutation_data)}")
				o.write(f"Sum of Mutations: {len(mutation_data)}")
			except:
				print("\nWarning. Unable to summarize pipeline statistics, or cannot open or write to file " + pipeline_summary + "\n")
		else:
			print("\nWarning. Missing stats files. Cannot summarize pipeline results\n")
	
def ref_length(ref_FilePath):
	try:
		with open(ref_FilePath, 'rt') as ref_file:
			ReadLines = ref_file.readlines()
	except:
		raise Exception("Cannot open ref file " + ref_FilePath + "\n")

	Lines = len(ReadLines)
	if Lines < 2:
		raise Exception("Unexpected error, empty or missing lines in ref file " + ref_FilePath + "\n")

	if ReadLines[0].startswith(">"):
		ref_genome = " "
		for i in range(1, Lines):
			if i == 1:
				ref_genome = ReadLines[i].strip().upper()
			else:
				ref_genome += ReadLines[i].strip().upper()
	else:
		raise Exception("first line in ref fasta file " + ref_FilePath + " does not start with >\n")

	ref_genome_length = len(ref_genome)
	return ref_genome_length

def count_coverage_positions (dir_path, freqs_file_path, Coverage, ref_FilePath):
	pipeline_summary = dir_path + "/Summary.txt"
	with open(pipeline_summary, "a") as o:
		try:
			ref_genome_length = ref_length(ref_FilePath)
			o.write("\nNumber of positions in reference genome: {}\n".format(int(ref_genome_length)))
			coverage_key = "x" + str(Coverage)
			data = pd.read_csv(freqs_file_path, sep = ',')
			data = data.drop_duplicates("ref_position")
			data[coverage_key] = np.where(data["coverage"] > Coverage, 1, 0)
			coverage_stats = pd.DataFrame.sum(data, axis = 0).get(key = coverage_key)
			o.write("Number of positions with min coverage x{}: {}\n".format(str(Coverage), int(coverage_stats)))
			percentage_coverage = round((int(coverage_stats)/int(ref_genome_length))*100,2)
			o.write("% of positions with min coverage x{}: {}\n".format(str(Coverage), percentage_coverage))
		except:
			print("\nWarning. Unable to count positions with " + Coverage + " coverage, or cannot open or write to file " + pipeline_summary + "\n")

def create_mutation_rate_csv(dir_path, freqs_file_path, Coverage, sample_basename_pattern, PROBABILITY_LIMIT):
	mutation_rates = dir_path + "/" + sample_basename_pattern + 'mutations_frequencies.csv'
	with open(mutation_rates, "w") as mutations_file:
		try:
			data = pd.read_csv(freqs_file_path, sep=',')
			data["mutation"] = data["ref_base"] + ">" + data["base"]
			data = pd.DataFrame.groupby(data[(data["rank"] != 0) & (data["coverage"] > Coverage) & (data["probability"] >= PROBABILITY_LIMIT)], by=["mutation"]).sum()  # & (data["base"] != "-") to also remove deletions
			data["mutation_frequency"] = round(data["base_counter"] / data["coverage"], 6)
			data.drop(columns=["ref_position", "base_counter", "coverage", "frequency", "probability", "rank"], inplace=True)
			mutations_file.write(data.to_csv())
		except:
			pipeline_summary = dir_path + "/Summary.txt"
			with open(pipeline_summary, "a") as o:
				o.write("\nUnable to create mutation rate csv file. No mutations were found in positions with >=" + str(Coverage) + " coverage answering the requested number of repeats and q-score\n")

def pipeline_statistics_plots(dir_path, sample_basename_pattern, freqs_file_path, Coverage, PROBABILITY_LIMIT, lower_ylim=10**-5, upper_ylim=10**-1):
	try:
		os.chdir(dir_path)

		fig, axes = plt.subplots(figsize=(25,10), ncols=3, nrows=2)
		plt.suptitle(sample_basename_pattern + "Pipeline statistics", fontsize=18)
		plt.subplots_adjust(wspace=0.2,hspace=0.25)
		sns.set_style("whitegrid")

		get_lengths = []
		get_lengths.append(subprocess.getoutput("awk -F '\t' '{print $8}' *blast"))
		lengths = get_lengths[0].split("\n")
		lengths = [int(i) for i in lengths]
		sns.distplot(lengths, axlabel="Length (bp)", kde=True, color="royalblue", ax=axes[0][0])
		axes[0][0].set_ylabel("Density")
		axes[0][0].title.set_text('Length distribution')

		get_read_ids = []
		get_read_ids.append(subprocess.getoutput("awk -F '\t' '{print $1}' *blast"))
		read_ids = get_read_ids[0].split("\n")
		MATCH_STATISTICS = {}
		for read_id in read_ids:
			if read_id not in MATCH_STATISTICS:
				MATCH_STATISTICS[read_id] = 0
			MATCH_STATISTICS[read_id] += 1
		matches = list(MATCH_STATISTICS.values())
		match_max = max(matches)
		sns.distplot(matches, axlabel="Num of blast matches", kde=False, color="deeppink", ax=axes[0][1])
		axes[0][1].set_xticks(list(range(1, match_max + 1)))
		axes[0][1].set_ylabel("Count")
		axes[0][1].title.set_text('Blast match distribution')

		data = pd.read_csv(freqs_file_path, sep=',')
		data["mutation"] = data["ref_base"] + ">" + data["base"]
		mutation_data = data[(data["rank"] != 0) & (data["coverage"] > Coverage) & (data["probability"] >= PROBABILITY_LIMIT)]
		sns.boxplot("mutation", "frequency", data=mutation_data, ax=axes[1][0])
		axes[1][0].set_yscale("log")
		axes[1][0].set_ylim(lower_ylim, upper_ylim)
		axes[1][0].set_xlabel("Mutation")
		axes[1][0].set_ylabel("Frequency")
		axes[1][0].title.set_text("Mutations frequencies by group")

		sns.scatterplot(x="ref_position", y="frequency", data=mutation_data, ax=axes[1][1])
		axes[1][1].set_yscale("log")
		axes[1][1].set_ylim(lower_ylim, upper_ylim)
		axes[1][1].set_xlabel("Ref position")
		axes[1][1].set_ylabel("Frequency")
		axes[1][1].title.set_text(f"Mutations frequencies by position (sum of mutations: {len(mutation_data)})")

		data = data.drop_duplicates("ref_position")
		InsertionIndex = data[data['ref_base'] == "-"].index
		data.drop(InsertionIndex, inplace=True)
		sns.lineplot(x='ref_position', y='coverage', data=data, color='lightseagreen', ax=axes[0][2])
		axes[0][2].set_yscale("log")
		axes[0][2].set_xlabel("Ref position")
		axes[0][2].set_ylabel("Coverage")
		axes[0][2].title.set_text("Coverage")

		plt.savefig(os.path.join(dir_path, sample_basename_pattern + 'pipeline_plots.png'), format='png', bbox_inches="tight", pad_inches=0.5)
		plt.close()

	except:
		pipeline_summary = dir_path + "/Summary.txt"
		with open(pipeline_summary, "a") as o:
			o.write("\nUnable to create pipeline statistics plots\n")

def main(args):
	dir_path = args.output_dir
	if not os.path.isdir(dir_path):
		raise Exception("input directory " + dir_path + " does not exist or is not a valid directory\n")

	Coverage = args.coverage
	try:
		Coverage = int(Coverage) 
	except:
		raise Exception("Unexpected error, number of reads per file " + Coverage + " is not a valid integer value\n")

	find_file_name = FindFilesInDir(dir_path, ".part1.fasta")
	if len(find_file_name) > 0:
		if "L00" in find_file_name[0]:
			sample_basename_pattern = os.path.basename(find_file_name[0].split("L00")[0])
		else:
			raise Exception("Unexpected error, was not able to find a common path for sample name. Unable to perform Join step\n")
	else:
		raise Exception("Unexpected error, was not able to find *part1.fasta files in directory " + dir_path + ". Unable to perform Join step\n")

	file_type = ".stats"
	stats_files = FindFilesInDir(dir_path, file_type)
	file_type = ".fasta"
	fasta_files = FindFilesInDir(dir_path, file_type)
	file_type = "merge.freqs.csv"
	freqs_file_path = FindFilesInDir(dir_path, file_type)
	print(f"freqs_file_path: {freqs_file_path}")
	if len(stats_files) > 0 or len(fasta_files) > 0:
		summarize_stats(dir_path, sample_basename_pattern, freqs_file_path[0], Coverage)
	else:
		print("Warning. stats files or fasta files are missing in directory " + dir_path + ". Cannot perform summary analysis\n")

	ref_FilePath = args.ref
	if not (os.path.isfile(ref_FilePath) and os.path.splitext(ref_FilePath)[1] == '.fasta'):
		raise Exception("Unexpected error, " + ref_FilePath + " does not exist, is not a file or is not a fasta file\n")

	
	if len(freqs_file_path) == 1:
		create_mutation_rate_csv(dir_path, freqs_file_path[0], Coverage, sample_basename_pattern, PROBABILITY_LIMIT)
		count_coverage_positions(dir_path, freqs_file_path[0], Coverage, ref_FilePath)
	else:
		print("Warning. number of merge.freqs.csv file in directory " + dir_path + " is different than 1. Cannot perform summary analysis of coverage and mutation rates\n")

	file_type = ".blast"
	blast_files = FindFilesInDir(dir_path, file_type)
	if len(blast_files) > 0 and len(freqs_file_path) == 1:
		pipeline_statistics_plots(dir_path, sample_basename_pattern, freqs_file_path[0], Coverage, PROBABILITY_LIMIT)
	else:
		print("Warning. blast files or freqs file are missing in directory " + dir_path + ". Cannot create summary pipeline plots\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory containing stats and merge.freqs.csv file", required=True)
	parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference seq file (fasta)", required=True)
	parser.add_argument("-c", "--coverage", type=int, help="coverage cut-off for statistics, default = 10000", required=False, default=10000)
	args = parser.parse_args()
	main(args)
