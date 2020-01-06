#! python/python-anaconda3.2019.7

import glob
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
import os
import subprocess
plt.switch_backend('agg')

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

def summarize_stats (dir_path):
	pipeline_summary = dir_path + "/Summary.txt"
	with open(pipeline_summary, "a") as o:	
		os.chdir(dir_path)
		file_type = ".stats"
		input_stats_files = FindFilesInDir(dir_path, file_type)
		if len(input_stats_files) > 0:		
			try:
				#only_once_reads = subprocess.getoutput("grep '1\t' *.stats -h | awk '{sum+=$2}; END {print sum}'")
				num_reads = subprocess.getoutput("cat *fasta | grep '^>' | wc -l")
				num_contributing_reads = subprocess.getoutput("grep 'Number of reads contributing to frequency counts' *stats -h | awk -F '\t' '{sum+=$2}; END{print sum}'" )
				num_contributing_bases = subprocess.getoutput("grep 'Number of bases contributing to frequency counts' *stats -h | awk -F '\t' '{sum+=$2}; END{print sum}'" )
				num_non_contributing_reads = subprocess.getoutput("grep 'Number of non contributing reads' *stats -h | awk -F '\t' '{sum+=$2}; END{print sum}'" ) 
				num_mapped_reads = int(num_contributing_reads) + int(num_non_contributing_reads)
				percentage_mapped = round((int(num_mapped_reads)/int(num_reads))*100,2)
				percentage_contribution = round((int(num_contributing_reads)/int(num_mapped_reads))*100,2)
				o.write("Total number of reads: {}\n".format(int(num_reads)))
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
			except:
				print("\nWarning. Unable to summarize pipeline statistics, or cannot open or write to file " + pipeline_summary + "\n")
		else:
			print("\nWarning. Missing stats files. Cannot summarize pipeline results\n")
	
def count_coverage_positions (dir_path, freqs_file_path, Coverage):
	pipeline_summary = dir_path + "/Summary.txt"
	with open(pipeline_summary, "a") as o:	
		coverage_key = "x" + str(Coverage)
		try:
			data = pd.read_csv(freqs_file_path, sep = ',')
			data = data.drop_duplicates("ref_position")
			data[coverage_key] = np.where(data["coverage"] > Coverage, 1, 0)
			coverage_stats = pd.DataFrame.sum(data, axis = 0).get(key = coverage_key)	
			o.write("\nNumber of positions with min coverage x{}: {}\n".format(Coverage, coverage_stats))	
		except:
			print("\nWarning. Unable to summarize pipeline statistics, or cannot open or write to file " + pipeline_summary + "\n")

def create_coverage_plot (dir_path, freqs_file_path):
	try:
		data = pd.read_csv(freqs_file_path, sep= ',')
		data = data.drop_duplicates("ref_position")
		sns.set_style("whitegrid", {'axes.axisbelow': False})
		plt.plot(data['ref_position'].values, data['coverage'].values, color='teal', zorder = 5)
		plt.title("Coverage", fontsize=16)
		plt.xlabel("Position in the genome (bp)")
		plt.yscale("log")
		plt.ylabel("Read count")
		plt.savefig(os.path.join(dir_path, 'coverage.png'), format='png')
		plt.close()
	except:
		pipeline_summary = dir_path + "/Summary.txt"	
		with open(pipeline_summary, "a") as o:
			o.write("Unable to create coverage plot\n")

def create_mutation_rate_plot (dir_path, freqs_file_path, Coverage, lower_ylim = 10**-5, upper_ylim = 10**-1):
	try:
		data = pd.read_csv(freqs_file_path, sep = ',')
		data["mutation"] = data["ref_base"] + ">" + data["base"]
		sns.set_style("whitegrid")
		ax=sns.boxplot("mutation", "frequency", data=data[(data["rank"] != 0) & (data["coverage"] > Coverage)]) 
		#ax=sns.boxplot("mutation", "frequency", data=data[(data["rank"] != 0) & (data["ref_base"] != data["base"]) & (data["coverage"] > Coverage)])
		ax.set_yscale("log")
		ax.set_ylim(lower_ylim, upper_ylim)
		plt.title("Mutation Rates", fontsize=16)
		plt.savefig(os.path.join(dir_path, 'mutation_rate.png'), format='png')
		plt.close()
		data = pd.DataFrame.groupby(data[(data["rank"] != 0) & (data["coverage"] > Coverage)], by = ["mutation"]).sum()	#& (data["base"] != "-") to also remove deletions
		#data = pd.DataFrame.groupby(data[(data["rank"] != 0) & (data["ref_base"] != data["base"]) & (data["coverage"] > Coverage)], by = ["mutation"]).sum()	#& (data["base"] != "-") to also remove deletions
		data["mutation_frequency"] = round(data["base_counter"]/data["coverage"],6)
		pd.DataFrame.drop(data, columns=["ref_position","base_counter","coverage","frequency","probability","rank"], inplace=True)
		mutation_rates = dir_path + "/mutation_rates.csv"
		with open(mutation_rates, "w") as o:
			o.write(data.to_csv())

	except:	
		pipeline_summary = dir_path + "/Summary.txt"	
		with open(pipeline_summary, "a") as o:
			o.write("Unable to create mutation rate plot. No mutations were found in positions with >=" + str(Coverage) + " coverage answering the requested number of repeats and q-score\n")

def main(args):
	dir_path = args.output_dir
	if not os.path.isdir(dir_path):
		raise Exception("input directory " + dir_path + " does not exist or is not a valid directory\n")

	Coverage = args.coverage
	try:
		Coverage = int(Coverage) 
	except:
		raise Exception("Unexpected error, number of reads per file " + Coverage + " is not a valid integer value\n") 
	
	file_type = ".stats"
	stats_files = FindFilesInDir(dir_path, file_type)
	if len(stats_files) > 0:
		summarize_stats(dir_path)
	else:
		print("Warning. stats files are missing in directory " + dir_path + ". Cannot perform summary analysis\n")
	
	file_type = "merge.freqs.csv"
	freqs_file_path = FindFilesInDir(dir_path, file_type)
	if len(freqs_file_path) == 1:
		count_coverage_positions(dir_path, freqs_file_path[0], Coverage)
		create_coverage_plot(dir_path, freqs_file_path[0])
		create_mutation_rate_plot(dir_path, freqs_file_path[0], Coverage)
	else:
		print("Warning. number of merge.freqs.csv file in directory " + dir_path + " is different than 1. Cannot perform summary analysis\n")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory containing stats and merge.freqs.csv file", required=True)
	parser.add_argument("-c", "--coverage", type=int, help="coverage cut-off for statistics, default = 10000", required=True)
	args = parser.parse_args()
	main(args)
