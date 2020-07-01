#! python/python-anaconda3.2019.7

import argparse
import time
import datetime
import os
import sys
#append parent directory to path so that we can import from there:
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.runner_utils import FindFilesInDir, check_queue, Sleep, create_array
from utils.pbs_jobs import create_pbs_cmd, submit
from utils.logger import pipeline_logger

#TODO?: relative path support for input files

''' 
Pipeline generating frequency files from raw sequencing data files, either fastq.gz or fastq.
0.	Merge files. Receiving a full path to a sample with pattern (dir1/dir2/sample1_) and a directory, where merged file will be created.
1.  Split files. Receiving an input directory containing the (merged or single read) files, either fastq.gz or fastq, and a directory where split files will be created. 
	If needed, convert fastq.gz to fastq files, and split into N equally sized smaller fasta files containing 25K reads as default, or Num_of_reads_per_file as user input. 
2.  Blast. Running blast on each of the N fasta files, either sequences vs. reference (default) mode or reference against sequences (as indexed db) mode.
    Blast task of choice, either blastn (default), megablast or dc-megablast (dc-megablast available in sequence to reference mode only).
3.  BaseCalling. Running base calling on each of the N blast files and creating a frequency output file based on the requested parameters (repeats and q-score).
4.  Join. Combine N frequency files to a single frequencies file. Create linked mutations file. Create consensus.
5.  Summarize. Create basic plots and summarize data from N stats files.
6.	Warp up. zip files.
'''

def merge_fastq_files(sample_dir_path, sample_basename_pattern, number_of_N, dir_path, queue):
	alias = "MergeFiles"
	script_path = "/sternadi/home/volume1/shared/SternLab/scripts/merge_fastq_files.py"

	files_to_merge = FindFilesInDir(sample_dir_path, sample_basename_pattern)
	num_of_files_to_merge = len(files_to_merge)
	if num_of_files_to_merge == 0 or num_of_files_to_merge % 2 != 0:
		raise Exception("Unexpected error, number of files to merge with pattern " + sample_basename_pattern + " in directory " + sample_dir_path + " is zero, or does not divide by 2\n")
	num_of_expected_merged_files = int(num_of_files_to_merge/2)

	R1_files = [];	R2_files = []
	for file in files_to_merge:
		if "_R1" in file:
			R1_files.append(file)
		elif "_R2" in file:
			R2_files.append(file)
		else:
			raise Exception("File" + file + "does not contain R1 or R2 in sample name. Unable to merge file\n")

	if len(R1_files) != len(R2_files) or len(R1_files) != num_of_expected_merged_files:
		raise Exception("Number of R1 sample files " + str(len(R1_files)) + " is different than number of R2 sample files " + str(len(R2_files)) + \
						" or from the number of expected merged files " + str(num_of_expected_merged_files) + "\n")

	files_to_merge_list = []
	for i in range(num_of_expected_merged_files):
		output_file_basename = os.path.basename(files_to_merge[2*i]).replace('_001', '').replace('_R1', '_merged')
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
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_expected_merged_files, gmem=gmem, cmds=cmds, queue=queue, load_python=True)
	job_id = submit(cmdfile)
	Sleep(alias, job_id)

	time.sleep(10)
	file_type = "merged*"
	merged_files = FindFilesInDir(dir_path, file_type)
	if len(merged_files) != num_of_expected_merged_files:
		raise Exception("Unexpected error, number of merged files does not match expected number of merged files " + str(num_of_expected_merged_files) + "\n")

def toFastaAndSplit(pipeline_dir, dir_path, input_files, Num_reads_per_file, queue):
	alias = "toFastaAndSplit"
	script_path = pipeline_dir + "/ToFastaAndSplit.py"

	num_of_input_files = len(input_files)
	if num_of_input_files == 0:
		raise Exception("Unexpected error, there are no input files in input directory " + dir_path + "\n")
	elif num_of_input_files == 1:
		index = '0'
		gmem = 2
	else:
		index = '$PBS_ARRAY_INDEX-1'
		gmem = 7

	array = create_array(input_files)
	cmd1 = 'declare -a FILENAMES\n'
	cmd2 = 'FILENAMES=' + array + "\n\n"
	cmd3 = "python " + script_path + " -o " + dir_path + " -f ${FILENAMES[" + index + "]} -n " + str(Num_reads_per_file) + "\n"
	cmds = cmd1 + cmd2 + cmd3

	cmdfile = dir_path + "/FastaAndSplit.cmd"
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_input_files, gmem=gmem, cmds=cmds, queue=queue, load_python=True)
	job_id = submit(cmdfile)
	Sleep(alias, job_id)

	time.sleep(10)
	file_type = ".fasta"
	fasta_files = FindFilesInDir(dir_path, file_type)
	file_type = ".qual"
	quality_files = FindFilesInDir(dir_path, file_type)
	if len(fasta_files) != len(quality_files) or len(fasta_files) == 0 or len(quality_files) == 0:
		raise Exception("Unexpected error, number of fasta and / or quality output files does not match expected number of output files\n")

def Blast (dir_path, ref_genome, task, mode, e_value, ID_blast, queue):
	alias = "Blast"
	blast_dir = "/sternadi/home/volume1/shared/tools/ncbi-blast-2.2.30+/bin"

	file_type = "part*.fasta"
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
				"-num_alignments 1000000 -dust no -soft_masking F -perc_identity " + str(ID_blast) + " -evalue " + str(e_value) + "\n"
		cmds = cmd1 + cmd2 + cmd3 + cmd4

	else:
		raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")

	cmdfile = dir_path + "/Blast.cmd"
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_input_files, gmem=gmem, cmds=cmds, queue=queue, load_python=True)
	job_id = submit(cmdfile)
	Sleep(alias, job_id)

	time.sleep(10)
	file_type = ".blast"
	blast_files = FindFilesInDir(dir_path, file_type)
	if len(blast_files) != len(input_fasta_files):
		raise Exception("Unexpected error, number of blast output files " + str(len(blast_files)) + " does not match number of input fasta files " + str(len(input_fasta_files)) + "\n")

def BaseCall(pipeline_dir, dir_path, ref_genome, min_num_repeats, q_score, mode, Protocol, queue):
	alias = "BaseCalling"
	script_path = pipeline_dir + "/BaseCall.py"

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
	cmd3 = "python " + script_path + " -f ${FILENAMES[" + index + "]} -r " + ref_genome + " -q " + str(q_score) + " -x " + str(min_num_repeats) + " -m " + mode + " -p " + Protocol
	cmds = cmd1 + cmd2 + cmd3

	cmdfile = dir_path + "/BaseCalling.cmd"
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_input_files, gmem=gmem, cmds=cmds, queue=queue, load_python=True)
	job_id = submit(cmdfile)
	Sleep(alias, job_id)

	time.sleep(30)
	file_type = ".freqs"
	freqs_files = FindFilesInDir(dir_path, file_type)
	if len(input_blast_files) != len(freqs_files):
		raise Exception("Unexpected error, number of freqs output files " + str(len(freqs_files)) + " does not match number of input blast files " + str(len(input_blast_files)) + "\n")

def Join (pipeline_dir, dir_path, ref_genome, Coverage, queue):
	alias = "Join"
	script_path = pipeline_dir + "/Join.py"

	file_type = ".freqs"
	input_freqs_files = FindFilesInDir(dir_path, file_type)
	if len(input_freqs_files) == 0:
		raise Exception("Unexpected error, there are no freqs files in directory " + dir_path + "\n")

	cmdfile = dir_path + "/Join.cmd"
	cmds = "python " + script_path + " -o " + dir_path + " -r " + ref_genome + " -c " + str(Coverage)
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=1, gmem=2, cmds=cmds, queue=queue, load_python=True)
	job_id = submit(cmdfile)
	Sleep(alias, job_id)

	time.sleep(10)
	file_type = "merge.freqs.csv"
	merge_file = FindFilesInDir(dir_path, file_type)
	if len(merge_file) != 1:
		raise Exception("Unexpected error, merge.freqs.csv file does not exist in directory " + dir_path + " or is different than 1\n")

def Summary (pipeline_dir, dir_path, Coverage, ref_genome, queue):
	alias = "Summary"
	script_path = pipeline_dir + "/Summary.py"

	file_type = "merge.freqs.csv"
	input_csv_files = FindFilesInDir(dir_path, file_type)
	if len(input_csv_files) != 1:
		raise Exception("Unexpected error, number of merge.freqs.csv file in directory " + dir_path + " is different than 1. Cannot perform summary analysis\n")

	file_type = ".stats"
	input_stats_files = FindFilesInDir(dir_path, file_type)
	if len(input_stats_files) < 1:
		raise Exception("Unexpected error, there are no stats files in directory " + dir_path + ". Cannot perform summary analysis\n")

	file_type = ".blast"
	input_blast_files = FindFilesInDir(dir_path, file_type)
	if len(input_blast_files) < 1:
		raise Exception("Unexpected error, there are no blast files in directory " + dir_path + ". Cannot perform summary analysis\n")

	cmdfile = dir_path + "/Summary.cmd"
	cmds = "python " + script_path + " -o " + dir_path + " -c " + str(Coverage) + " -r " + ref_genome
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=1, gmem=2, cmds=cmds, queue=queue, load_python=True)
	job_id = submit(cmdfile)
	Sleep(alias, job_id)

	time.sleep(10)
	file_type = "Summary.txt"
	summary_file = FindFilesInDir(dir_path, file_type)
	if len(summary_file) != 1:
		raise Exception("Unexpected error, number of Summary.txt file in directory" + dir_path + " is different than 1\n")

def main(args):


	pipeline_dir = os.path.dirname(os.path.abspath(__file__).strip())
	pipeline_path = pipeline_dir + "/Runner.py"

	if not (os.path.isfile(pipeline_path) and os.path.splitext(pipeline_path)[1] == '.py'):
		raise Exception("Unexpected error, pipeline path " + pipeline_path + " does not exist or not a .py file\n")

	start_stage = args.start
	if start_stage != None:
		if start_stage not in [0,1,2,3,4,5,6]:
			raise Exception("Unexpected error, start_stage " + str(start_stage) + " is not a valid value. Valid values are integers between 0-6\n")

	if args.path != None:
		sample_path = args.path.strip()
		sample_dir_path = os.path.dirname(sample_path)
		sample_basename_pattern = os.path.basename(sample_path)
	else:
		sample_path =""; sample_dir_path = ""; sample_basename_pattern = ""

	dir_path = args.output_dir.strip()
	if start_stage in [0, 1]:
		if not os.path.isdir(dir_path):
			try:
				os.system("mkdir -p " + dir_path)
			except:
				raise Exception("failed to create directory " + dir_path + "\n")
	if not os.path.isdir(dir_path):
		raise Exception("Directory " + dir_path + " does not exist or is not a valid directory path\n")
	if args.log_folder is not None:
		log_folder = args.log_folder
	else:
		log_folder = dir_path
	log = pipeline_logger('PipelineRunner', log_folder)
	if start_stage == None or start_stage in [0,1]:
		if not os.path.isdir(sample_dir_path):
			raise Exception("Directory of sample to merge/split " + sample_dir_path + " does not exist or is not a valid directory path\n")
		if sample_basename_pattern == "":
			raise Exception("Sample basename pattern is missing in path " + sample_path + "\n")
		else:
			sample_basename_pattern += "*"
		file_type = sample_basename_pattern + "_R2*"
		paired_samples = FindFilesInDir(sample_dir_path, file_type)
		file_type = "merged.fastq*"
		merged_files = FindFilesInDir(dir_path, file_type)
		if len(paired_samples) > 0 and start_stage == None:
			start_stage = 0
		elif len(paired_samples) == 0 and start_stage == None:
			start_stage = 1
		elif len(paired_samples) > 0 and start_stage > 0 and len(merged_files) == 0:
			raise Exception("R2 reads for sample pattern " + os.path.basename(sample_path) + " were identified in directory " + sample_dir_path + ". Pipeline should start at stage 0\n")
		elif len(paired_samples) == 0 and start_stage == 0:
			raise Exception("R2 reads for sample pattern " + os.path.basename(sample_path) + " were not identified in directory " + sample_dir_path + ". Pipeline should start at stage 1\n")

	end_stage = args.end
	if end_stage != None:
		if end_stage not in [0,1,2,3,4,5,6]:
			raise Exception("Unexpected error, end_stage is not a valid integer value between 0-6. Valid values are integers between 0-6\n")

	if start_stage > end_stage:
		raise Exception ("Unexpected error, start stage " + str(start_stage) + " is larger than end stage " + str(end_stage) + "\n")

	number_of_N = args.num_of_N
	if start_stage == 0 and number_of_N != None:
		if number_of_N < 60:
			log.warning("Running merge reads with number_of_N smaller than 60")

	ref_genome = args.ref.strip()
	if not (os.path.isfile(ref_genome) and os.path.splitext(ref_genome)[1] == '.fasta'):
		raise Exception("Unexpected error, " + ref_genome + " reference genome file does not exist, is not a file or is not a fasta file\n")

	q_score = args.q_score
	if q_score != None:
		if q_score < 0 or q_score > 40:
			raise Exception("Unexpected error, q-score value " + str(q_score) + " is not valid, should be an integer value between 0-40\n")
		if q_score < 16:
			log.warning("Running pipeline with q-score value of " + str(q_score))

	blast_id = args.blast_id
	if blast_id != None:
		if blast_id < 85:
			log.warning("Running pipeline with blast id value of " + str(blast_id))
		if blast_id < 0 or blast_id > 100:
			raise Exception("Unexpected error, identity % for blast is not a valid value: " + str(blast_id) + " \n")

	e_value = args.evalue
	if e_value != None:
		if e_value > 1e-7:
			log.warning("Running pipeline with e_value > " + str(e_value) + "\n")

	task = args.blast_task.strip()
	if task != None:
		if task not in ["megablast", "blastn", "dc-megablast"]:
			raise Exception("Unexpected error, blast task has to be 'blastn', 'megablast' or 'dc-megablast'\n")

	mode = args.blast_mode.strip()
	if mode != None:
		if mode not in ["ReftoSeq", "RS", "rs", "sr", "SR", "SeqtoRef"]:
			raise Exception("Unexpected error, blast mode has to be either ReftoSeq, RS, rs or SeqtoRef, SR, sr\n")

	queue = args.queue.strip()
	check_queue(queue)

	min_num_repeats = args.repeats
	if min_num_repeats != None:
		if min_num_repeats < 1:
			raise Exception ("Unexpected error, min number of repeats is less than 1\n")
		if min_num_repeats < 2:
			log.warning("Running pipeline with min number of repeats less than 2")
		if min_num_repeats > 2:
			log.warning("Running pipeline with min number of repeats bigger than 2")

	Min_Num_reads_per_file = 10000
	Max_Num_reads_per_file = 40000
	Num_reads_per_file = args.num_reads
	if Num_reads_per_file != None:
		if Num_reads_per_file < Min_Num_reads_per_file:
			log.warning("Running pipeline with less than " + str(Min_Num_reads_per_file) + " reads per split file")
		if Num_reads_per_file > Max_Num_reads_per_file:
			log.warning("Running pipeline with more than " + str(Max_Num_reads_per_file) + " reads per split file")

	Min_Coverage = 1000
	Coverage = args.coverage
	if Coverage != None:
		if Coverage < Min_Coverage:
			log.warning("Running pipeline with coverage smaller than " + str(Min_Coverage))

	Protocol = args.protocol.strip()
	if Protocol not in ["L", "l", "linear", "C", "c", "circular"]:
		raise Exception("Unexpected error, for linear library prep protocol type 'linear', 'L' or 'l', for circular library prep protocol type 'circular', 'C' or 'c'\n")

	overwrite = args.overwrite.strip()
	if overwrite not in ["Y", "y", "N", "n"]:
		raise Exception("Unexpected error, overwrite value is different than Y/N\n")

	cmd = "python {} -i {} -o {} -r {} -m {} -t {} -s {} -e {} -q {} -d {} -v {} -x {} -n {} -c {} -p {} -u {} -w {}".format(pipeline_path, sample_path, dir_path, ref_genome, mode, task,
                                                            start_stage, end_stage, q_score, blast_id, e_value, min_num_repeats, Num_reads_per_file, Coverage, Protocol, queue, overwrite)
	log.info(cmd)

	pipeline_summary = dir_path + "/Summary.txt"
	try:
		with open(pipeline_summary, "a") as o:
			o.write("\n---- Pipeline running -----\n")
			o.write("{}\n\n".format(datetime.datetime.now()))
			o.write("Pipeline command used:\n{}\n\n".format(cmd))
			o.write("Blast parameters: mode = {}, task = {}, %id for blast = {}, e-value = {}\n".format(mode, task, blast_id, e_value))
			o.write("Base Calling Parameters: number of repeats used = {}, q-score = {}, protocol = {}\n\n".format(min_num_repeats, q_score, Protocol))
	except:
		raise Exception("Unexpected error, cannot write into file " + pipeline_summary + "\n")

	for stage in range(start_stage, end_stage+1):
		if stage == 0:
			file_type = "merged.fastq*"
			merged_files = FindFilesInDir(dir_path, file_type)
			if len(merged_files) > 0 and overwrite in ["N", "n"]:
				raise Exception("Unexpected error, merged files were found in directory" + dir_path + ". Cannot run merge_fastq_files. To re-write add -w Y to command line\n")
			elif len(merged_files) > 0	and overwrite in ["Y","y"]:
				os.system("rm -rf " + dir_path + "/*.merged.fastq*")
			merge_fastq_files(sample_dir_path, sample_basename_pattern, number_of_N, dir_path, queue)

		if stage == 1:
			file_type = ".fasta"
			fasta_files = FindFilesInDir(dir_path, file_type)
			file_type = ".qual"
			quality_files = FindFilesInDir(dir_path, file_type)
			if (len(fasta_files) > 0 or len(quality_files) > 0) and overwrite == "N":
				raise Exception("Unexpected error, fasta and / or quality files were found in directory " + dir_path + ". Cannot run toFastaAndSplit. To re-write add -w Y to command line\n")
			elif (len(fasta_files) > 0 or len(quality_files) > 0) and overwrite == "Y":
				os.system("rm -rf " + dir_path + "/*.fasta " + dir_path + "/*.qual")

			file_type = sample_basename_pattern + "_R2*"
			paired_samples = FindFilesInDir(sample_dir_path, file_type)
			file_type = "merged.fastq*"
			merged_files = FindFilesInDir(dir_path, file_type)
			if start_stage == 1 and len(paired_samples) == 0:	#for single-end reads use sample_dir_path to get files for split
				file_type = ".fastq*"
				log.info("Single-end read, fetching files from {}\n".format(sample_dir_path))
				sample_files = FindFilesInDir(sample_dir_path, file_type)
			elif len(merged_files) > 0:		#for paired-end reads use dir_path to get merged files for split
				log.info("Paired-end read, fetching files from {}\n".format(dir_path))
				sample_files = FindFilesInDir(dir_path, file_type)
			else:
				raise Exception("Unable to detect sample files with base name pattern " + sample_basename_pattern + " in directory\n")
			input_gz_files = []; input_fastq_files = []
			for file in sample_files:
				if os.path.splitext(file)[1] == ".gz":
					input_gz_files.append(file)
				elif os.path.splitext(file)[1] == ".fastq":
					input_fastq_files.append(file)
				else:
					raise Exception("Sample pattern " + os.path.basename(sample_path) + " is not a .gz or .fastq file")
			if len(input_fastq_files) > 0 and len(input_gz_files) > 0 and overwrite in ["N", "n"]:
				raise Exception("Unexpected error, there is a mix of fastq and gz files in directory " + dir_path + "\n")
			elif len(input_fastq_files) > 0 and len(input_gz_files) > 0 and overwrite in ["Y", "y"]:
				os.system("rm -rf " + dir_path + "/*.fastq")
				input_files = input_gz_files
			elif len(input_fastq_files) == 0 and len(input_gz_files) == 0:
				raise Exception("Unexpected error, there are no fastq or gz files in directory " + dir_path + "\n")
			elif len(input_fastq_files) > 0 and len(input_gz_files) == 0:
				input_files = input_fastq_files
			elif len(input_fastq_files) == 0 and len(input_gz_files) > 0:
				input_files = input_gz_files
			toFastaAndSplit(pipeline_dir, dir_path, input_files, Num_reads_per_file, queue)

		if stage == 2:
			file_type = ".blast"
			blast_files = FindFilesInDir(dir_path, file_type)
			if len(blast_files) > 0 and overwrite in ["N", "n"]:
				raise Exception("Unexpected error, blast files were found in directory " + dir_path + ". Cannot run Blast. To re-write add -w Y to command line\n")
			elif len(blast_files) > 0 and overwrite in ["Y","y"]:
				os.system("rm -rf " + dir_path + "/*.blast")
			Blast(dir_path, ref_genome, task, mode, e_value, blast_id, queue)

		if stage == 3:
			file_type = ".freqs"
			freqs_files = FindFilesInDir(dir_path, file_type)
			if len(freqs_files) > 0 and overwrite in ["N", "n"]:
				raise Exception("Unexpected error, freqs files were found in directory " + dir_path + ". Cannot run BaseCall. To re-write add -w Y to command line\n")
			elif len(freqs_files) > 0 and overwrite in ["Y","y"]:
				os.system("rm -rf " + dir_path + "/*.freqs* " + dir_path + "/*.good* " + dir_path + "/*.NonContributing")
			BaseCall(pipeline_dir, dir_path, ref_genome, min_num_repeats, q_score, mode, Protocol, queue)

		if stage == 4:
			file_type = "merge.freqs.csv"
			merged_freqs = FindFilesInDir(dir_path, file_type)
			if len(merged_freqs) > 0 and overwrite in ["N", "n"]:
				raise Exception("Unexpected error, merge.freqs.csv was found in directory " + dir_path + ". Cannot run Join. To re-write add -w Y to command line\n")
			elif len(merged_freqs) > 0 and overwrite in ["Y","y"]:
				os.system("rm -rf " + dir_path + "/*.merge.freqs.csv")
			Join(pipeline_dir, dir_path, ref_genome, Coverage, queue)

		if stage == 5:
			Summary(pipeline_dir, dir_path, Coverage, ref_genome, queue)

		if stage == 6:
			os.system("zip " + dir_path + "/OutputFiles.zip " + dir_path + "/*.part* " + dir_path + "/*.OU")
			os.system("rm -rf " + dir_path + "/*.part* " + dir_path + "/*.OU")

	log.info("END OF PIPELINE RUN")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--path", type=str, help="a full path with a pattern of sample name for merge and/or split. Example: dir1/dir2/sample1_", required=False)
	parser.add_argument("-N", "--num_of_N", type=int, help="number of N's to add for merge of R1 and R2 pair-end reads", required=False, default=60)
	parser.add_argument("-n", "--num_reads", type=int, help="number of reads per split file, default=25000", required=False, default=25000)
	parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory", required=True)
	parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference fasta file", required=True)
	parser.add_argument("-m", "--blast_mode", type=str, help="mode for blast, for Seq to Ref blast type SR, sr or SeqtoRef, for Ref to Seq blast type RS, rs or ReftoSeq, default = 'SeqtoRef'",
						required=False, default="SeqtoRef")
	parser.add_argument("-d", "--blast_id", type=int, help="% blast id, default=85", required=False, default=85)
	parser.add_argument("-t", "--blast_task", type=str, help="task for blast, blastn/megablast/dc-megablast?, default='blastn'", required=False, default="blastn")
	parser.add_argument("-v", "--evalue", type=float, help="E-value for blast, default=1e-7", required=False, default=1e-7)
	parser.add_argument("-s", "--start", type=int, help="start stage number", required=False)
	parser.add_argument("-e", "--end", type=int, help="end stage number, default=5", required=False, default=5)
	parser.add_argument("-q", "--q_score", type=int, help="Q-score cutoff, default=30", required=False, default=30)
	parser.add_argument("-x", "--repeats", type=int, help="number of repeats, default=2", required=False, default=2)
	parser.add_argument("-c", "--coverage", type=int, help="coverage cut-off for statistics, default=10000", required=False, default=10000)
	parser.add_argument("-p", "--protocol", type=str, help="Library prep protocol is linear = 'L', 'l' or 'linear', or circular = 'C', 'c' or 'circular'. Default='linear'",
						required=False, default="linear")
	parser.add_argument("-u", "--queue", type=str, help="queue to run pipeline, default='tzachi@power9'", required=False, default="tzachi@power9")
	parser.add_argument("-w", "--overwrite", type=str, help="overwrite? Y/N, default='N'", required=False, default="N")
	parser.add_argument("-L", "--log_folder", type=str, help="Folder path to write .log file in, default is None", default=None)
	args = parser.parse_args()
	main(args)


