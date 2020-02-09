#! python/python-anaconda3.2019.7

import argparse
import glob
import os
import time
	
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

def SplitToSmallerFiles(dir_path, FilePath, Num_reads_per_file): 
	FileSuffix = os.path.splitext(FilePath)[1]
	if FileSuffix == '.gz': 
		os.system("zcat " + FilePath + " > " + dir_path + "/" + os.path.basename(FilePath).split(".gz")[0])
		FilePath = dir_path + "/" + os.path.basename(FilePath).split(".gz")[0]
	elif FileSuffix == '.fastq':
		pass
	else:
		raise Exception("Unexpected error. Wrong file. Input file should be either fastq or gz file and not " + FileSuffix + "\n")        
	
	try:
		Lines = int(os.popen("awk 'END {print NR}' " + FilePath).read())        
	except:
		raise Exception ("Unexpected error, number of lines in file " + FilePath + " is not an integer value\n")
			
	if Lines == 0 or Lines % 4 != 0:
		raise Exception ("Unexpected error, file " + FilePath + " is empty, or number of lines in file does not divide by 4\n")

	Num_reads_in_file = Lines/4
	Num_files = int(Num_reads_in_file/Num_reads_per_file)+1
	if (Num_reads_in_file % Num_reads_per_file) == 0:
		Num_files -= 1	

	ReadCount = 1
	File_counter = 1
	Bases = ['A', 'T', 'C', 'G', 'N']    
	ASCII_signs = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', \
            '!', '"', '#', '$', '%', '&',"'", '(', ')', '*', '+', ',', '-', '.', '/', ':', ';', '<', '=', '>', '@', '?'] 
	
	split_file_basename = dir_path + "/" + os.path.basename(FilePath).split(".fastq")[0]
	with open(FilePath, 'rt') as ReadRecords:
		while File_counter <= Num_files: 
			Curr_out_file = split_file_basename + ".part" + str(File_counter) + ".fasta"
			Curr_qual_out_file = split_file_basename + ".part" + str(File_counter) + ".qual"
			try:
				with open(Curr_out_file, 'at') as FastaFile:
					try:
						with open(Curr_qual_out_file, 'at') as QualityFile:
							while (ReadCount <= Num_reads_per_file*File_counter) and (ReadCount <= Num_reads_in_file):
								for i in range(4):
									if i == 0:
										read_id = ReadRecords.readline().strip()
										if "R1" in os.path.basename(FilePath):
											read_id = read_id.split(" ")[0].strip()
										if not read_id.startswith("@"):
											raise Exception("Unexpected error, missing @. Cannot identify read id in first line of read id\n")
									if i == 1:
										seq_line = ReadRecords.readline().strip()
										seq_line_length = len(seq_line)
										for position in range(seq_line_length):
											if not seq_line[position] in Bases:
												raise Exception("Unexpected error at position " + str(position+1) + " in read id " + read_id + ". Cannot identify a valid DNA letter\n")
									if i == 2:
										plus = ReadRecords.readline().strip()
										if not plus.startswith("+"):
											raise Exception("Unexpected error, missing +. Cannot identify plus sign in third line of read id\n")
									if i == 3:
										qual_line = ReadRecords.readline().strip()
										qual_line_length = len(qual_line)
										for position in range(qual_line_length):
											if not qual_line[position] in ASCII_signs:
												raise Exception("Unexpected error at position " + str(position+1) + " in read id " + read_id + ". Cannot identify a valid ASCII letter\n")

								try:
									QualityFile.write(read_id + "\n" + qual_line + "\n")
								except:
									raise Exception("Unexpected error, cannot write into file " + Curr_qual_out_file + "\n")

								read_id = read_id.replace("@",">")
								try:
									FastaFile.write(read_id + "\n" + seq_line + "\n")
								except:
									raise Exception("Unexpected error, cannot write into file " + Curr_out_file + "\n")

								ReadCount += 1
					except:
						raise Exception("Unexpected error, cannot open file " + Curr_qual_out_file + "\n")
			except:
				raise Exception("Unexpected error, cannot open file " + Curr_out_file + "\n")
			
			File_counter += 1 
			
	time.sleep(10)

	file_type = os.path.basename(FilePath).split(".fastq")[0] + "*.fasta"
	fasta_files = FindFilesInDir(dir_path, file_type)
    
	file_type = os.path.basename(FilePath).split(".fastq")[0] + "*.qual"
	quality_files = FindFilesInDir(dir_path, file_type)
	
	if len(fasta_files) != Num_files or len(quality_files) != Num_files: 
		raise Exception("Unexpected error, number of fasta and / or quality output files does not match expected number of output files " + str(Num_files) + "\n")

def main(args):
	dir_path = args.out_dir
	if not os.path.isdir(dir_path):
		raise Exception("/ninput directory " + dir_path + " does not exist or is not a valid directory\n")
	
	FilePath = args.file
	if not (os.path.isfile(FilePath) and ((os.path.splitext(FilePath)[1] == '.gz') or (os.path.splitext(FilePath)[1] == '.fastq'))):
		raise Exception("/ninput file " + FilePath + " does not exist or is not a valid gz. or fastq. file\n")

	Min_Num_reads_per_file = 10000
	Max_Num_reads_per_file = 40000
	Num_reads_per_file = args.num_reads
	if Num_reads_per_file != None:
		if Num_reads_per_file < Min_Num_reads_per_file:
			print("\nWarning. Running pipeline with less than " + str(Min_Num_reads_per_file) + " reads per split file\n")
		if Num_reads_per_file > Max_Num_reads_per_file:
			print("\nWarning. Running pipeline with more than " + str(Max_Num_reads_per_file) + " reads per split file\n")

	SplitToSmallerFiles(dir_path, FilePath, Num_reads_per_file)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-o", "--out_dir", type=str, help = "a path to an output directory where split files are saved", required=True)
	parser.add_argument("-f", "--file", type=str, help = "a path to a gz or fastq file to split to smaller fasta and qual files", required=True)
	parser.add_argument("-n", "--num_reads", type=int, help = "number of reads per split file, default=25000", required=False, default=25000)
	args = parser.parse_args()
	main(args)
