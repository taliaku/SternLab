#! /powerapps/share/Modules/Centos7/modulefiles/python/python-anaconda3.2019.7

import argparse
import time
import datetime
#import pbs_jobs
#import pipeline_pbs_jobs
import os
import glob

''' 
Pipeline generating frequency files from raw sequencing data files, either fastq.gz or fastq.
0.	Merge files. Recieving a full path to a sample with sample name written as pattarn (dir1/dir2/sample1_\*) and an input_dir, where merged file will be created.
1.  Split files. Recieving an input direcory containing the (merged or single read) files, either fastq.gz or fastq, and a dir_path where split files will be created. 
	If needed, convert fastq.gz to fastq files, and split into N equally sized smaller fasta files containing 25K reads as default, or Num_of_reads_per_file as input. 
2.  Blast. Running blast on each of the N fasta files, either sequences vs. reference (default) mode or reference against sequences (as indexed db) mode.
    Task of choice, either blastn (default), megablast or dc-megablast (dc-megablast available in sequence to reference mode only).
3.  BaseCalling. Running base calling on each of the N blast files and creating a frequency output file based on the requested parameters (repeats and q-score).
4.  Join. Combine N frequency files to a single frequencies file. Create linked mutations file. Create concensus.
5.  Summerize. Create basic plots and summerize data from N stats files.
6.	Warp up. zip files.
'''

def create_pbs_cmd(cmdfile, alias, jnum, gmem, cmds, load_python=True, queue="adis"):
	with open(cmdfile, 'w') as o:
		o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
		o.write("#PBS -q %s\n" % queue)
		o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
		o.write("#PBS -N "+ alias+"\n")
		o.write("#PBS -o %s\n" % "/".join(cmdfile.split("/")[:-1]))
		o.write("#PBS -e %s\n" % "/".join(cmdfile.split("/")[:-1]))
		if gmem:
			mem=gmem*1000
			o.write("#PBS -l mem="+str(mem)+"mb\n")
		if jnum:
			if jnum != 1:
				o.write("#PBS -J 1-"+str(jnum)+"\n\n")
		o.write("id\n")
		o.write("date\n")
		o.write("hostname\n")
		if load_python:
			o.write("module load python/python-anaconda3.2019.7\n")       
		o.write("\n")
		o.write(cmds)
		o.write("\n")
		o.write("date\n")
	o.close()

def submit(cmdfile):
	cmd = "/opt/pbs/bin/qsub " + cmdfile
	result = os.popen(cmd).read()
	if 'power' in result:
		return result.split(".")[0]
	else:
		print(cmdfile + " was not submitted")	

def Sleep (alias, job_id, sleep_max = 1200000, sleep_quantum = 10):
	i = 0 
	process = os.popen("qstat -t " + job_id + " | wc -l").read()
	try:
		process = int(process)
	except:
		process = 0
	
	while process > 0 and i <= sleep_max: 
		time.sleep(sleep_quantum)
		i += sleep_quantum
		print ("Running...")
		process = os.popen("qstat -t " + job_id + " | wc -l").read()
		try:
			process = int(process)
		except:
			process = 0
		
	if process > 0: 
		raise Exception(alias + " stage was not completed. Max sleep time reached\n")

def FindFilesInDir(dir_path, file_type):
	file_path = dir_path + "/*" + file_type
	list_of_files = sorted(glob.glob(file_path))
	num_of_files = len(list_of_files)
	if num_of_files > 0:
		for file in list_of_files:
			size = os.path.getsize(file)
			if size == 0:
				time.sleep(15)
				size = os.path.getsize(file)
			if size == 0:
				raise Exception("Unexpected error, some of the " + file_type + " files in " + dir_path + " are empty\n")
       
	return list_of_files
	
def create_array(files_list):
	array = '(' 
	for i in range(len(files_list)):
		array += files_list[i]
		if i != len(files_list)-1:
			array += " "
	array += ')'
	
	return array

def merge_fastq_files(sample_dir_path, sample_basename_pattern, number_of_N, dir_path):
	alias = "MergeFiles"
	script_path = "/sternadi/home/volume1/shared/SternLab/scripts/merge_fastq_files.py"
	
	#sample_dir_path = os.path.dirname(sample_name_pattarn)
	#sample_basename_pattern = os.path.basename(sample_name_pattarn) + "*"
	files_to_merge = FindFilesInDir(sample_dir_path, sample_basename_pattern)	
	num_of_files_to_merge = len(files_to_merge)
	if num_of_files_to_merge == 0 or num_of_files_to_merge % 2 != 0:
		raise Exception("Unexpected error, number of files to merge with pattern " + sample_basename_pattern + " in directory " + sample_dir_path + " is zero, or does not divide by 2\n")
	
	num_of_expected_merged_files = int(num_of_files_to_merge/2)
	files_to_merge_list = []
	for i in range(num_of_expected_merged_files):
		#if 'R1_001' in os.path.basename(files_to_merge[2*i]):
		output_file_basename = os.path.basename(files_to_merge[2*i]).replace('R1_001','merged')
		#else:
			#output_file_basename = os.path.basename(files_to_merge[2*i]).replace('R1','merged')		
		files_to_merge_list.extend((files_to_merge[2*i], files_to_merge[2*i+1], dir_path + "/" + output_file_basename))
	
	if num_of_expected_merged_files == 1:
		f = '0'
		e = '1'
		o = '2'
		gmem = 2
	else:
		f = '$((PBS_ARRAY_INDEX*3))-3'
		e = '$((PBS_ARRAY_INDEX*3))-2'
		o = '$((PBS_ARRAY_INDEX*3))-1'
		gmem = 7
	
	array = create_array(files_to_merge_list)
	
	cmd1 = 'declare -a FILENAMES\n'
	cmd2 = 'FILENAMES=' + array + "\n\n"
	cmd3 = "python " + script_path + " -f ${FILENAMES[" + f + "]} -e ${FILENAMES[" + e + "]} -o ${FILENAMES[" + o + "]} -r " + str(number_of_N) + "\n"
	cmds = cmd1 + cmd2 + cmd3

	cmdfile = dir_path + "/merge_files.cmd"	
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_expected_merged_files, gmem=gmem, cmds=cmds, load_python=True)
	job_id = submit(cmdfile)
	print(job_id)
	Sleep(alias, job_id)
	
	time.sleep(10)
	file_type = "merged*"
	merged_files = FindFilesInDir(dir_path, file_type)
	if len(merged_files) != num_of_expected_merged_files: 
		raise Exception("Unexpected error, number of merged files does not match expected number of merged files " + str(num_of_expected_merged_files) + "\n")
	
def toFastaAndSplit(dir_path, input_files, Num_reads_per_file):
	alias = "toFastaAndSplit"
	script_path = "/sternadi/home/volume2/yaara/temp/temp/ToFastaAndSplit_ver7.py"
	
	num_of_input_files = len(input_files)
	if num_of_input_files == 0:
		raise Exception("Unexpected error, there are no input files in input directory " + dir_path + "\n") #move to main? what about 2 samples types, gz and fastq?
	elif num_of_input_files == 1:
		index = '0'
		gmem = 2
	else:
		index = '$PBS_ARRAY_INDEX-1'
		gmem = 7
  
	array = create_array(input_files)
	cmd1 = 'declare -a FILENAMES\n'
	cmd2 = 'FILENAMES=' + array + "\n\n"
	cmd3 = "python " + script_path + " -o " + dir_path + " -f ${FILENAMES[" + index + "]} -n " + str(Num_reads_per_file) + "\n" #{arr[3]}
	cmds = cmd1 + cmd2 + cmd3
	
	cmdfile = dir_path + "/FastaAndSplit.cmd"
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_input_files, gmem=gmem, cmds=cmds, load_python=True)
	job_id = submit(cmdfile)
	print(job_id)
	Sleep(alias, job_id)
	
	time.sleep(10)
	file_type = ".fasta"
	fasta_files = FindFilesInDir(dir_path, file_type)    
	file_type = ".qual"
	quality_files = FindFilesInDir(dir_path, file_type)
	if len(fasta_files) != len(quality_files) or len(fasta_files) == 0 or len(quality_files) == 0: 
		raise Exception("Unexpected error, number of fasta and / or quality output files does not match expected number of output files\n")

def Blast (dir_path, ref_genome, task, mode, e_value, ID_blast):
	alias = "Blast"
	blast_dir = "/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin"

	file_type = ".fasta"
	input_fasta_files = FindFilesInDir(dir_path, file_type)
	num_of_input_files = len(input_fasta_files)
	if num_of_input_files == 0:
		raise Exception("Unexpected error, there are no fasta files in directory" + dir_path + "\n") 
	elif num_of_input_files == 1:
		index = '0'
		gmem = 2
	else:
		index = '$PBS_ARRAY_INDEX-1'
		gmem = 7

	array = create_array(input_fasta_files)
	cmd1 = 'declare -a FILENAMES\n'
	cmd2 = 'FILENAMES=' + array + "\n\n"
	
	if mode in ["SeqtoRef", "SR", "sr"]:
		outfile = '${FILENAMES[' + index + ']}.SR.' + task + '.blast'
		cmd3 = blast_dir + "/blastn -query ${FILENAMES[" + index + "]} -task " + task + " -out " + outfile + " -subject " + ref_genome + " -outfmt '6 qseqid sseqid qstart qend qstrand sstart send sstrand length btop' " \
				"-num_alignments 1000000 -dust no -soft_masking F -perc_identity " + str(ID_blast) + " -evalue " + str(e_value) + "\n"
		cmds = cmd1 + cmd2 + cmd3
	
	elif mode in ["ReftoSeq", "RS", "rs"]:	
		outfile = '${FILENAMES[' + index + ']}.RS.' + task + '.blast'
		cmd3 = blast_dir + '/makeblastdb -in ${FILENAMES[' + index + ']} -dbtype nucl\n'
		cmd4 = blast_dir + "/blastn -query " + ref_genome + " -task " + task + " -out " + outfile + " -db ${FILENAMES[" + index + "]} -outfmt '6 qseqid sseqid qstart qend qstrand sstart send sstrand length btop' " \
				"-num_alignmnts 1000000 -dust no -soft_masking F -perc_identity " + str(ID_blast) + " -evalue " + str(e_value) + "\n"
		cmds = cmd1 + cmd2 + cmd3 + cmd4
	
	else:
		raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")
	
	cmdfile = dir_path + "/Blast.cmd"
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_input_files, gmem=gmem, cmds=cmds, load_python=True)
	job_id = submit(cmdfile)
	print(job_id)
	Sleep(alias, job_id)
	
	time.sleep(10)
	file_type = ".blast"
	blast_files = FindFilesInDir(dir_path, file_type)
	if len(blast_files) != len(input_fasta_files): 
		raise Exception("Unexpected error, number of blast output files " + str(len(blast_files)) + " does not match number of input fasta files " + str(len(input_fasta_files)) + "\n")				
		
def BaseCall(dir_path, ref_genome, min_num_repeats, q_score, mode, Protocol):
	alias = "BaseCalling"
	script_path = "/sternadi/home/volume2/yaara/temp/temp/BaseCall_33.py"
		
	file_type = ".blast"
	input_blast_files = FindFilesInDir(dir_path, file_type)
	num_of_input_files = len(input_blast_files)
	if num_of_input_files == 0:
		raise Exception("Unexpected error, there are no blast files in directory" + dir_path + "\n") 
	elif num_of_input_files == 1:
		index = '0'
		gmem = 2
	else:
		index = '$PBS_ARRAY_INDEX-1'
		gmem = 7
	
	array = create_array(input_blast_files)
	cmd1 = 'declare -a FILENAMES\n'
	cmd2 = 'FILENAMES=' + array + "\n\n"
	cmd3 = "python " + script_path + " -f ${FILENAMES[" + index + "]} -r " + ref_genome + " -q " + str(q_score) + " -rep " + str(min_num_repeats) + " -m " + mode + " -pr " + Protocol
	cmds = cmd1 + cmd2 + cmd3
	    
	cmdfile = dir_path + "/BaseCalling.cmd"
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_input_files, gmem=gmem, cmds=cmds, load_python=True)
	job_id = submit(cmdfile)
	print(job_id)
	Sleep(alias, job_id)
	
	time.sleep(30)	
	file_type = ".freqs"
	freqs_files = FindFilesInDir(dir_path, file_type)	
	if len(input_blast_files) != len(freqs_files):  
		raise Exception("Unexpected error, number of freqs output files " + str(len(freqs_files)) + " does not match number of input blast files " + str(len(input_blast_files)) + "\n")				
	
def Join (dir_path, path, ref_genome, Coverage):
	alias = "Join"
	script_path = "/sternadi/home/volume2/yaara/temp/temp/Join6.py"
	
	file_type = ".freqs"
	input_freqs_files = FindFilesInDir(dir_path, file_type)
	if len(input_freqs_files) == 0:
		raise Exception("Unexpected error, there are no freqs files in directory " + dir_path + "\n") 
	
	cmdfile = dir_path + "/Join.cmd"
	cmds = "python " + script_path + " -o " + dir_path + " -p " + path + " -r " + ref_genome + " -c " + str(Coverage)
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=1, gmem=2, cmds=cmds, load_python=True)
	job_id = submit(cmdfile)
	print(job_id) 
	Sleep(alias, job_id)
	
	time.sleep(10)
	file_type = "merge.freqs.csv"
	merge_file = FindFilesInDir(dir_path, file_type)
	if len(merge_file) != 1:
		raise Exception("Unexpected error, merge.freqs.csv file does not exist in directory " + dir_path + " or is different than 1\n")
		
def Summary (dir_path, Coverage):
	alias = "Summary"
	script_path = "/sternadi/home/volume2/yaara/temp/temp/summary7.py"
	
	file_type = "merge.freqs.csv"
	input_csv_files = FindFilesInDir(dir_path, file_type)
	if len(input_csv_files) != 1:
		raise Exception("Unexpected error, number of merge.freqs.csv file in directory " + dir_path + " is different than 1. Cannot perform summary analysis\n")
		
	file_type = ".stats"
	input_stats_files = FindFilesInDir(dir_path, file_type)
	if len(input_stats_files) < 1:
		raise Exception("Unexpected error, there are no stats files in directory " + dir_path + ". Cannot perform summary analysis\n") 
	
	cmdfile = dir_path + "/Summary.cmd"
	cmds = "python " + script_path + " -o " + dir_path + " -c " + str(Coverage)
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=1, gmem=2, cmds=cmds, load_python=True)
	job_id = submit(cmdfile)
	print(job_id)    
	Sleep(alias, job_id)
	
	time.sleep(10)
	file_type = "Summary.txt"
	summary_file = FindFilesInDir(dir_path, file_type)
	if len(summary_file) != 1:
		raise Exception("Unexpected error, number of Summary.txt file in directory" + dir_path + " is different than 1\n")
		
def main(args):
	
	pipeline_path = '/sternadi/home/volume2/yaara/temp/temp/runner_ver14.py'
	
	start_stage = args.start
	if start_stage != None:
		try:
			start_stage = int(start_stage) 
		except:
			raise Exception("Unexpected error, start_stage " + start_stage + " is not a valid value\n")
		if not start_stage in [0,1,2,3,4,5,6]:
			raise Exception("Unexpected error, start_stage " + start_stage + " is not a valid value\n") 	
    
	end_stage = args.end
	if end_stage != None:
		try:
			end_stage = int(end_stage) 
		except:
			raise Exception("Unexpected error, end_stage " + end_stage + " is not a valid integer value between 0-6")   
		if not end_stage in [0,1,2,3,4,5,6]:
			raise Exception("Unexpected error, end_stage is not a valid integer value between 0-6\n")
    
	if (start_stage > end_stage):
		raise Exception ("Unexpected error, start stage " + str(start_stage) + " is larger than end stage " + str(end_stage) + "\n")
			
	path = args.path
	if path != None:
		sample_dir_path = os.path.dirname(path)
		if not os.path.isdir(sample_dir_path):
			raise Exception("Directory for merge files " + sample_dir_path + " does not exist or is not a valid directory path\n")
		sample_basename_pattern = os.path.basename(path) + "*"
	
	dir_path = args.output_dir
	number_of_N = args.num_of_N
	if start_stage == 0:		
		if number_of_N != None:
			try:
				number_of_N = int(number_of_N) 
			except:
				raise Exception("Unexpected error, number of Ns per merge file " + str(number_of_N) + " is not a valid integer value\n")  
		if not os.path.isdir(dir_path):
			try:
				os.system("mkdir -p " + dir_path)
			except:
				raise Exception("failed to create directory " + dir_path + "\n")
			#if not os.path.isdir(dir_path):
				#raise Exception("Directory " + dir_path + " does not exist or is not a valid directory path\n")	
	
	if not os.path.isdir(dir_path):	#check, consider single read which would start here!!! 
		raise Exception("Directory " + dir_path + " does not exist or is not a valid directory path\n")

	ref_genome = args.ref
	if not (os.path.isfile(ref_genome) and os.path.splitext(ref_genome)[1] == '.fasta'):
		raise Exception("Unexpected error, " + ref_genome + " reference genome file does not exist, is not a file or is not a fasta file\n")
    	
	q_score = args.q_score
	if q_score != None:
		try:
			q_score = int(q_score) 
		except:
			raise Exception("Unexpected error, q_score should be an integer value between 0-40\n")
		if q_score < 0 or q_score > 40: 
			raise Exception("Unexpected error, q-score value " + q_score + " is not valid, should be an integer value between 0-40\n")
		if q_score < 16:
			print("Warning, running pipeline with q-score value of " + q_score + "\n")
	
	blast_id = args.blast_id
	if blast_id != None:
		try:
			blast_id = int(blast_id) 
		except:
			raise Exception("Unexpected error, identity % for blast should be an integer value\n")
		if blast_id < 0 or blast_id > 100:
			raise Exception("Unexpected error, identity % for blast is not a valid value: " + blast_id + " \n")
	
	e_value = args.evalue
	if e_value != None:
		try:
			e_value = float(e_value) 
		except:
			raise Exception("Unexpected error, e-value " + e_value + " is not a valid value\n")
		if e_value < 1e-7:
			print("Warning, running pipeline with e_value < " + e_value + "\n")
	
	task = args.blast_task
	if task != None:
		if task not in ["megablast", "blastn", "dc-megablast"]:
			raise Exception("Unexpected error, blast task has to be 'blastn', 'megablast' or 'dc-megablast'\n") 	
			
	mode = args.blast_mode
	if mode != None:
		if mode not in ["ReftoSeq", "RS", "rs", "sr", "SR", "SeqtoRef"]:
			raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n") 

	min_num_repeats = args.repeats
	if min_num_repeats != None:
		try:
			min_num_repeats = int(args.repeats)
		except:
			raise Exception("Unexpected error, min number of repeats " + min_num_repeats + " is not a valid integer value\n")
		if min_num_repeats < 1: 
			raise Exception ("Unexpected error, min number of repeats is less than 1\n")
		if min_num_repeats < 2:  
			print("Warning. Running pipeline with min number of repeats less than 2\n")
		if min_num_repeats > 2:  
			print("Warning. Running pipeline with min number of repeats bigger than 2\n")
			
	Num_reads_per_file = args.num_reads
	if Num_reads_per_file != None:	
		try:
			Num_reads_per_file = int(Num_reads_per_file) 
		except:
			raise Exception("Unexpected error, number of reads per file " + Num_reads_per_file + " is not a valid integer value\n")   
        
	Coverage = args.coverage
	if Coverage != None:
		try:
			Coverage = int(Coverage) 
		except:
			raise Exception("Unexpected error, number of reads per file " + Coverage + " is not a valid integer value\n")   
			
	Protocol = args.protocol
	if Protocol == None:
		Protocol = "linear"
	else:
		if Protocol not in ["L", "l", "linear", "C", "c", "circular"]:
			raise Exception("Unexpected error, for linear library prep protocol type 'linear' or 'L', for circular library prep protocol type 'circular' or 'C'\n") 
	
	cmd = "python {} -p {} -o {} -r {} -m {} -t {} -s {} -e {} -q {} -id {} -ev {} -rep {} -n {} -c {} -pr {}".format(pipeline_path, path, dir_path, ref_genome, mode, task, 
                                                            start_stage, end_stage, q_score, blast_id, e_value, min_num_repeats, Num_reads_per_file, Coverage, Protocol)
	print(cmd)
		
	pipeline_summary = dir_path + "/Summary.txt"
	try:
		with open(pipeline_summary, "a") as o:
			o.write("\n---- Pipeline running -----\n")
			o.write("{}\n\n".format(datetime.datetime.now()))
			o.write("Pipeline command used:\n{}\n\n".format(cmd))
			o.write("Blast parameters: mode = {}, task = {}, %id for blast = {}, e-value = {}\n".format(mode, task, blast_id, e_value))
			o.write("Base Calling Parameters: number of repeats used = {}, q-score = {}, protocol = {}\n\n".format(min_num_repeats, q_score, Protocol))
	except:
		raise Exception("Unexpected error, cannot write into file " + pipeline_summary + "\n")	#print warning instead???
	
	for stage in range(start_stage, end_stage+1):	
		if stage == 0:
			merge_fastq_files(sample_dir_path, sample_basename_pattern, number_of_N, dir_path)
		if stage == 1:
			file_type = ".fastq"
			input_fastq_files = FindFilesInDir(dir_path, file_type)
			file_type = ".gz"
			input_gz_files = FindFilesInDir(dir_path, file_type)
			if (len(input_fastq_files) > 0 and len(input_gz_files) > 0):
				raise Exception("Unexpected error, there is a mix of fastq and gz files in directory " + dir_path + "\n")
			elif (len(input_fastq_files) == 0 and len(input_gz_files) == 0):
				raise Exception("Unexpected error, there are no fastq or gz files in directory " + dir_path + "\n")
			elif (len(input_fastq_files) > 0 and len(input_gz_files) == 0):
				input_files = input_fastq_files
			elif (len(input_fastq_files) == 0 and len(input_gz_files) > 0):
				input_files = input_gz_files	
			toFastaAndSplit(dir_path, input_files, Num_reads_per_file)
		if stage == 2:
			Blast(dir_path, ref_genome, task, mode, e_value, blast_id)	
		if stage == 3:
			BaseCall(dir_path, ref_genome, min_num_repeats, q_score, mode, Protocol) 
		if stage == 4:
			Join(dir_path, path, ref_genome, Coverage)
		if stage == 5:
			Summary(dir_path, Coverage) 
		if stage == 6:
			os.system("zip " + dir_path + "/OutputFiles.zip " + dir_path + "/*.part* " + dir_path + "/*.OU")
			os.system("rm -rf " + dir_path + "/*.part* " + dir_path + "/*.OU")
			
	print("END OF PIPELINE RUN")		
    
if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-p", "--path", type=str, help="a path with a pattarn of sample name for merge. Example: dir1/dir2/sample1_", required=False)
	parser.add_argument("-N", "--num_of_N", type=int, help="number of N's to add for merge of R1 and R2 pair-end reads", required=False, default=60) 
	parser.add_argument("-n", "--num_reads", type=int, help="number of reads per split file", required=False, default=25000) 
	parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory", required=True)
	parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference fasta file", required=True)
	parser.add_argument("-m", "--blast_mode", type=str, help="mode for blast, for Seq to Ref blast type SR, sr or SeqtoRef, for Ref to Seq blast type RS, rs or ReftoSeq, default = 'SeqtoRef'", required=False, default="SeqtoRef")
	parser.add_argument("-id", "--blast_id", type=int, help="% blast id, default=85", required=False, default=85)
	parser.add_argument("-t", "--blast_task", type=str, help="task for blast, blastn/megablast/dc-megablast?, default='blastn'", required=False, default="blastn")
	parser.add_argument("-ev", "--evalue", type=float, help="E-value for blast, default=1e-7", required=False, default=1e-7)
	parser.add_argument("-s", "--start", type=int, help="start step number. default=1", required=False, default=0)
	parser.add_argument("-e", "--end", type=int, help="end step number. default=6", required=False, default=6)
	parser.add_argument("-q", "--q_score", type=int, help="Q-score cutoff, default=30", required=False, default=30)
	parser.add_argument("-rep", "--repeats", type=int, help="number of repeats, default = 2", required=False, default=2)
	parser.add_argument("-c", "--coverage", type=int, help="coverage cut-off for statistics, default = 10000", required=False, default=10000)
	parser.add_argument("-pr", "--protocol", type=str, help="Library prep protocol is linear = 'L', 'l' or 'linear', or circular = 'C', 'c' or 'circular'. Default = 'linear'", required=False, default="linear")
	args = parser.parse_args()
	main(args)


