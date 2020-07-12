#! python/python-anaconda3.2019.7

import argparse
import pandas as pd
import glob
import os
import numpy as np

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
	
def create_ref_seq (ref_FilePath):	
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

	REF_GENOME = {}
	ref_genome_length = len(ref_genome)
	for i in range(ref_genome_length):
		if ref_genome[i] in ['A','C','G','T']:
			REF_GENOME[float(i+1)] = [ref_genome[i]]
		else:
			raise Exception("Found a non valid DNA letter in position " + ref_genome[i+1] + " of the reference genome\n")
		
	return REF_GENOME


def wrangle_freqs_df(data):
	data = pd.DataFrame.groupby(data, level=[0, 1, 2]).sum()
	data["frequency"] = round(data["base_counter"] / data["coverage"], 6)
	data["probability"] = round(1 - np.power(10, np.log10(1 - data["frequency"] + 1e-07) * (data["coverage"] + 1)), 2)
	# Get ranks
	pd.DataFrame.reset_index(data, level=[0, 1, 2], inplace=True)
	pd.DataFrame.sort_values(data, by=['ref_position', 'frequency', 'base'], ascending=[True, False, False],
							 inplace=True)
	data["coverage_to_set_rank"] = np.where(data["coverage"] > 0, 1, 0)
	data["rank"] = (pd.Series.cumsum(pd.Series(data["coverage_to_set_rank"])) - 1) % 5
	del data["coverage_to_set_rank"]
	# Create freqs file
	pd.DataFrame.set_index(data, keys="ref_position", inplace=True)
	return data


def JoinFreqs(dir_path, sample_basename_pattern, REF_GENOME, Minimal_insertion_coverage = 1):
	file_type = ".freqs"
	input_freqs_files = FindFilesInDir(dir_path, file_type)
	if len(input_freqs_files) == 0:
		raise Exception("Unexpected error, no freqs files were found in directory " + dir_path + "\n")
	dfs = []
	for path in input_freqs_files:
		with open(path, "r") as source_file:
			df = pd.read_csv(path, sep= '\t', index_col = [0, 1, 2], usecols=[0, 1, 2, 3, 4])
		dfs.append(df)
	data = pd.concat(dfs)
	data = wrangle_freqs_df(data)
	csv_file_path = dir_path + "/" + sample_basename_pattern + "merge.freqs.csv"
	try:
		with open(csv_file_path, 'w') as csv_file:
			csv_file.write(data.to_csv())
	except:
		raise Exception("Unexpected error, cannot write into file " + csv_file_path + "\n")

	#Adjust dataframe for merge with original ref
	pd.DataFrame.reset_index(data, inplace = True)
	data = data.drop_duplicates("ref_position")
	pd.DataFrame.set_index(data, keys = "ref_position", inplace = True)
	#Get indexes of insertions with <= minimal coverage to get a "cleaner" consensus ref
	InsertionCoverageIndex = data[(data['coverage'] <= Minimal_insertion_coverage) & (data['ref_base'] == "-")].index
	data.drop(InsertionCoverageIndex, inplace=True)

	#Create dataframe of original ref
	df_ref = pd.DataFrame.from_dict(REF_GENOME, orient='index')
	pd.DataFrame.rename_axis(df_ref, "ref_position", inplace = True)
	pd.DataFrame.rename(df_ref, columns={0: "ref_base"}, inplace = True)

	#Create consensus with indels
	consensus_df = pd.merge(df_ref, data, how = 'outer', on = 'ref_position', sort = False)
	pd.DataFrame.rename_axis(consensus_df, "ref_position", inplace = True)
	pd.DataFrame.sort_index(consensus_df, inplace = True)
	consensus_df["base_consensus"] = np.where(consensus_df["coverage"] > 0, consensus_df["base"], consensus_df["ref_base_x"])
	indels_consensus = ''.join(consensus_df["base_consensus"].values).replace('-','')
	consensus_file_path = dir_path + "/" + sample_basename_pattern + "consensus_with_indels.fasta"
	try:
		with open(consensus_file_path, 'w') as consensus_file:
			consensus_file.write(">" + os.path.basename(consensus_file_path).split(".fasta")[0] + "\n")
			consensus_file.write(indels_consensus)
	except:
		raise Exception("Unexpected error, cannot write into file " + consensus_file_path + "\n")

	#Create consensus without insertions to stay in frame with original ref
	consensus_df["base_consensus"] = np.where(consensus_df["base_consensus"] == '-', consensus_df["ref_base_x"], consensus_df["base_consensus"])  # ignores deletions, leaves original ref base
	InsertionIndex = consensus_df[consensus_df['ref_base_y'] == "-"].index	#removes insertions
	consensus_df.drop(InsertionIndex, inplace=True)
	no_indels_consensus = ''.join(consensus_df["base_consensus"].values)
	consensus_file_path = dir_path + "/" + sample_basename_pattern + "consensus_without_indels.fasta"
	try:
		with open(consensus_file_path, 'w') as consensus_file:
			consensus_file.write(">" + os.path.basename(consensus_file_path).split(".fasta")[0] + "\n")
			consensus_file.write(no_indels_consensus)
	except:
		raise Exception("Unexpected error, cannot write into file " + consensus_file_path + "\n")

def links_between_mutations(dir_path, sample_basename_pattern, min_value = 1):
	file_type = ".good_mutations"
	input_good_mutations_files = FindFilesInDir(dir_path, file_type)
	if len(input_good_mutations_files) == 0:
		raise Exception("Unexpected error, no good mutations files were found in directory " + dir_path + "\n")
	
	os.chdir(dir_path)
	good_mutations_all_file = dir_path + "/" + sample_basename_pattern + "good_mutations_all.txt"
	os.system("cat *.good_mutations > " + good_mutations_all_file)
	os.system("sort -Vk2,2 -Vk1,1 " + good_mutations_all_file + " -o " + good_mutations_all_file + "\n")	#sort file for linked mutations analysis is a must. Do not change sort values!!!
	try:
		with open(good_mutations_all_file, 'rt') as mutations_file:
			ReadLines = mutations_file.readlines()
	except:		
		raise Exception("Cannot open good mutations file\n")
	
	MUTATIONS_LINKS = {}
	if len(ReadLines) > 1: 
		i = 0
		while i < len(ReadLines)-1:
			while ReadLines[i].split("\t")[1].strip() != ReadLines[i+1].split("\t")[1].strip():
				i += 1
				if i == len(ReadLines)-1:
					break
			if i < len(ReadLines)-1:
				positions = ReadLines[i].split("\t")[0].strip()
				mutations = ReadLines[i].split("\t")[2].strip()
				while ReadLines[i].split("\t")[1].strip() == ReadLines[i+1].split("\t")[1].strip():
					positions = positions + "+" + ReadLines[i+1].split("\t")[0].strip()
					mutations = mutations + "+" + ReadLines[i+1].split("\t")[2].strip()
					i += 1
					if i == len(ReadLines)-1:
						break

				if positions not in MUTATIONS_LINKS:
					MUTATIONS_LINKS[positions] = {}
				if mutations not in MUTATIONS_LINKS[positions]:
					MUTATIONS_LINKS[positions][mutations] = [0,[]]
				MUTATIONS_LINKS[positions][mutations][0] += 1
				MUTATIONS_LINKS[positions][mutations][1].append(ReadLines[i].split("\t")[1].strip())

	linked_mutations_file = dir_path + "/" + sample_basename_pattern + "linked_mutations.txt"
	with open(linked_mutations_file, 'wt') as linked_mutations:
		for positions in MUTATIONS_LINKS:
			for mutations in MUTATIONS_LINKS[positions]:  
				if MUTATIONS_LINKS[positions][mutations][0] > min_value:
					try:
						linked_mutations.write(positions + "\t" + mutations + "\t" + str(MUTATIONS_LINKS[positions][mutations][0]) + "\t" + str(MUTATIONS_LINKS[positions][mutations][1]) + "\n") 
					except:
						raise Exception("Cannot write to file " + linked_mutations_file + "\n")	
					
	os.system("sort -k3nr " + linked_mutations_file + " -o " + linked_mutations_file + "\n")	

def main(args):
	dir_path = args.out_dir
	if not os.path.isdir(dir_path):
		raise Exception("Directory " + dir_path + " does not exist or is not a valid directory\n")

	find_files = FindFilesInDir(dir_path, ".part1.fasta")
	if len(find_files) > 0:
		if "L00" in find_files[0]:
			sample_basename_pattern = os.path.basename(find_files[0].split("L00")[0])
		else:
			raise Exception("Unexpected error, was not able to find a common path for sample name. Unable to perform Join step\n")
	else:
		raise Exception("Unexpected error, was not able to find *part1.fasta files in directory " + dir_path + ". Unable to perform Join step\n")

	ref_FilePath = args.ref
	if not (os.path.isfile(ref_FilePath) and os.path.splitext(ref_FilePath)[1] == '.fasta'):
		raise Exception("Unexpected error, " + ref_FilePath + " does not exist, is not a file or is not a fasta file\n")

	Min_Coverage = 1000
	Coverage = args.coverage
	if Coverage != None:
		if Coverage < Min_Coverage:
			print("\nWarning. Running pipeline with coverage smaller than " + str(Min_Coverage) + "\n")

	REF_GENOME = create_ref_seq(ref_FilePath)
	JoinFreqs(dir_path, sample_basename_pattern, REF_GENOME, Coverage)
	links_between_mutations(dir_path, sample_basename_pattern)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", "--out_dir", type=str, help = "a path to an output directory in which the results will be saved", required=True)
	parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference seq file (fasta)", required=True)
	parser.add_argument("-c", "--coverage", type=int, help="coverage cut-off for statistics analysis, default=10000", required=False, default=10000)
	args = parser.parse_args()
	main(args)
	

