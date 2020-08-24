# for each good reads file create list of good reads.
# for each item in list look for match in dict (or list) of reads.
# find index of read and count in dict.
# get count of items add to summary file.

import argparse
import os
import glob
import subprocess

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
    
def main(args):
	dir_path = args.output_dir.strip()
	if not os.path.isdir(dir_path):
		raise Exception("Directory " + dir_path + " does not exist or is not a valid directory path\n")
        
	os.chdir(dir_path)
    
	#file_type = ".good_reads"
	#good_reads_files = FindFilesInDir(dir_path, file_type)
	#num_of_input_files = len(input_fasta_files)
	#print(len(good_reads_files))
    
	get_read_ids = []
	get_read_ids.append(subprocess.getoutput("awk -F '\n' '{print $1}' *.good_reads"))
	read_ids = get_read_ids[0].split("\n")
		
	#fasta_file = file.split("part")[0] + "part" + file.split("part")[1].split(".")[0] + ".fasta"
	get_read_ids_fasta = []
	get_read_ids_fasta.append(subprocess.getoutput("awk -F '\n' '{print $1}' *.fasta"))
	read_ids_fasta = get_read_ids_fasta[0].split("\n")
	read_ids_fasta_keys = read_ids_fasta[::2]
	read_ids_fasta_values = read_ids_fasta[1::2]
	read_ids_fasta_dictionary = dict(zip(read_ids_fasta_keys, read_ids_fasta_values))
    
	primerID_STATISTICS = {}
	extension_seq = "TCCTGGTCGAGCTGGAC" #add seq
	primer_ID_length = 12
	N_ID = primer_ID_length * "N"
        
	for read_id in read_ids:
		read_id = read_id.replace("@",">")        
		extension_index = read_ids_fasta_dictionary[read_id].find(extension_seq)
		if extension_index > -1:
			primerID_start = extension_index+len(extension_seq)
			primerID_end = extension_index+len(extension_seq)+primer_ID_length
			primerID = read_ids_fasta[read_ids_fasta.index(read_id)+1][primerID_start:primerID_end+1]
			if primerID not in primerID_STATISTICS:
				primerID_STATISTICS[primerID] = 0
			primerID_STATISTICS[primerID] += 1
		
	if N_ID in primerID_STATISTICS:
		print("N_ID was found")
		del primerID_STATISTICS [N_ID]
    
	#os.chdir(dir_path)
	pipeline_summary = dir_path + "/Summary.txt"
	with open(pipeline_summary, "a") as o:
		o.write("\nTotal number of primerID contributing to frequency count: {}\n\n".format(len(primerID_STATISTICS.keys())))
	print("Number of primerIDs for sample " + os.path.basename(dir_path) + " " + str(len(primerID_STATISTICS.keys())))
        
        
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory", required=True)	
	args = parser.parse_args()
	main(args)