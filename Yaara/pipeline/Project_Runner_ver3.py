#! python/python-anaconda3.2019.7

import argparse
import time
import glob
import os

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
		print("cmd file was not submitted")	

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
		process = os.popen("qstat -t " + job_id + " | wc -l").read() #check
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

def run_project(pipeline_path, input_dir, dir_path, ref_genome, mode, task, start_stage, end_stage, q_score, blast_id, e_value, min_num_repeats, Num_reads_per_file, Coverage, Protocol):
	alias = "RunProject"

	file_type = "_R1*"
	samples_files = FindFilesInDir(input_dir, file_type)   
	if len(samples_files) == 0:
		raise Exception("Unexpected error, there are no input files in input directory " + input_dir + "\n")
	
	samples_list = []
	for sample in samples_files:
		if not '_L00' in sample:
			raise Exception("Unexpected error. '_L00' is missing in sample name, cannot identify sample name pattern\n")
		else:
			sample_name = os.path.basename(sample).split('_L00')[0]
		sample_pattern = input_dir + "/" + sample_name + "_"
		if not sample_pattern in samples_list:
			sample_dir_path = dir_path + "/" + sample_name
			if not os.path.isdir(sample_dir_path):
				try:
					os.system("mkdir -p " + sample_dir_path)
				except:
					raise Exception("failed to create input directory " + sample_dir_path + "\n")
				if not os.path.isdir(sample_dir_path):
					raise Exception("Sample directory " + sample_dir_path + " does not exist or is not a valid directory path\n")	
			samples_list.extend((sample_pattern, sample_dir_path))

	if len(samples_list) % 2 != 0:
		raise Exception("Unexpected error, number of samples for pipeline does not match number of output directory created\n")

	file_type = "_R2*"
	paired_samples = FindFilesInDir(input_dir, file_type)

	if start_stage == None:
		if len(paired_samples) > 0:		#paired-end reads
			start_stage = 0
		else: 	#single-end reads
			start_stage = 1
	if start_stage > end_stage:
		raise Exception ("Unexpected error, start stage " + str(start_stage) + " is larger than end stage " + str(end_stage) + "\n")

	num_of_samples = len(samples_files)
	if num_of_samples == 1:
		p = '0'
		o = '1'
		gmem = 2
	else:
		p = '$((PBS_ARRAY_INDEX*2))-2'
		o = '$((PBS_ARRAY_INDEX*2))-1'
		gmem = 7

	array = create_array(samples_list)
	cmd1 = 'declare -a SAMPLENAMES\n'
	cmd2 = 'SAMPLENAMES=' + array + "\n\n"
	cmd3 = "python " + pipeline_path + " -p ${SAMPLENAMES[" + p + "]} -o ${SAMPLENAMES[" + o + "]} -r " + ref_genome + " -m " + mode + " -t " + task + " -s " + str(start_stage) + " -e " + str(end_stage) + " -q " + str(q_score) + \
			" -id " + str(blast_id) + " -ev " + str(e_value) + " -rep " + str(min_num_repeats) + " -n " + str(Num_reads_per_file) + " -c " + str(Coverage) + " -pr " + Protocol + "\n"
	cmds = cmd1 + cmd2 + cmd3
	
	cmdfile = dir_path + "/pipeline_project_runner.cmd"	
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_samples, gmem=gmem, cmds=cmds, load_python=True)
	job_id = submit(cmdfile)
	print(job_id)
	Sleep(alias, job_id)
	
def main(args):

	pipeline_dir = os.path.dirname(os.path.abspath(__file__))
	pipeline_path = pipeline_dir + "/runner.py"

	if not os.path.isfile(pipeline_path):
		raise Exception("Unexpected error, " + pipeline_path + " does not exist, is not a file or or is not a blast file\n")
	
	start_stage = args.start
	if start_stage != None:
		try:
			start_stage = int(start_stage) 
		except:
			raise Exception("Unexpected error, start_stage " + str(start_stage) + " is not a valid value\n")
		if not start_stage in [0,1,2,3,4,5,6]:
			raise Exception("Unexpected error, start_stage " + str(start_stage) + " is not a valid value\n")
    
	end_stage = args.end
	if end_stage != None:
		try:
			end_stage = int(end_stage) 
		except:
			raise Exception("Unexpected error, end_stage " + str(end_stage) + " is not a valid integer value between 0-6")
		if not end_stage in [0,1,2,3,4,5,6]:
			raise Exception("Unexpected error, end_stage is not a valid integer value between 0-6\n")
    
	if start_stage != None and end_stage != None and start_stage > end_stage:
		raise Exception ("Unexpected error, start stage " + str(start_stage) + " is larger than end stage " + str(end_stage) + "\n")
			
	input_dir = args.input_dir
	if not os.path.isdir(input_dir):
		raise Exception("Directory for input files " + input_dir + " does not exist or is not a valid directory path\n")
	
	number_of_N = args.num_of_N
	if number_of_N != None:
		try:
			number_of_N = int(number_of_N) 
		except:
			raise Exception("Unexpected error, number of Ns per merge file " + number_of_N + " is not a valid integer value\n")

	dir_path = args.output_dir
	if not os.path.isdir(dir_path):
		try:
			os.system("mkdir -p " + dir_path)
		except:
			raise Exception("failed to create input directory " + dir_path + "\n")
		if not os.path.isdir(dir_path):
			raise Exception("Directory " + dir_path + " does not exist or is not a valid directory path\n")
	
	ref_genome = args.ref
	if not (os.path.isfile(ref_genome) and os.path.splitext(ref_genome)[1] == '.fasta'):
		raise Exception("Unexpected error, reference genome file does not exist, is not a file or is not a fasta file\n")
    	
	q_score = args.q_score
	if q_score != None:
		try:
			q_score = int(q_score) 
		except:
			raise Exception("Unexpected error, q_score should be an integer value between 0-40\n")
		if q_score < 0 or q_score > 40: 
			raise Exception("Unexpected error, q-score value " + str(q_score) + " is not valid, should be an integer value between 0-40\n")
		if q_score < 16:
			print("\nWarning, running pipeline with q-score value of " + str(q_score) + "\n")
	
	blast_id = args.blast_id
	if blast_id != None:
		try:
			blast_id = int(blast_id) 
		except:
			raise Exception("Unexpected error, identity % for blast should be an integer value\n")
		if blast_id < 0 or blast_id > 100:
			raise Exception("Unexpected error, identity % for blast is not a valid value: " + str(blast_id) + " \n")
	
	e_value = args.evalue
	if e_value != None:
		try:
			e_value = float(e_value) 
		except:
			raise Exception("Unexpected error, e-value " + e_value + " is not a valid value\n")
		if e_value < 1e-7:
			print("\nWarning, running pipeline with e_value < " + str(e_value) + "\n")
	
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
			raise Exception("Unexpected error, min number of repeats " + str(min_num_repeats) + " is not a valid integer value\n")
		if min_num_repeats < 1: 
			raise Exception ("Unexpected error, min number of repeats is less than 1\n")
		if min_num_repeats < 2:  
			print("\nWarning. Running pipeline with min number of repeats less than 2\n")
		if min_num_repeats > 2:  
			print("\nWarning. Running pipeline with min number of repeats bigger than 2\n")
			
	Num_reads_per_file = args.num_reads
	if Num_reads_per_file != None:	
		try:
			Num_reads_per_file = int(Num_reads_per_file) 
		except:
			raise Exception("Unexpected error, number of reads per file " + str(Num_reads_per_file) + " is not a valid integer value\n")
        
	Coverage = args.coverage
	if Coverage != None:
		try:
			Coverage = int(Coverage) 
		except:
			raise Exception("Unexpected error, number of reads per file " + str(Coverage) + " is not a valid integer value\n")

	Protocol = args.protocol
	if Protocol == None:
		Protocol = "linear"
	else:
		if Protocol not in ["L", "l", "linear", "C", "c", "circular"]:
			raise Exception("Unexpected error, for linear library prep protocol type 'linear' or 'L', for circular library prep protocol type 'circular' or 'C'\n") 			
	
	cmd = "python {} -i {} -o {} -r {} -m {} -t {} -s {} -e {} -q {} -id {} -ev {} -rep {} -n {} -c {} -pr {}".format(pipeline_path, input_dir, dir_path, ref_genome, mode, task, 
                                                            start_stage, end_stage, q_score, blast_id, e_value, min_num_repeats, Num_reads_per_file, Coverage, Protocol)
	print(cmd)
	
	run_project(pipeline_path, input_dir, dir_path, ref_genome, mode, task, start_stage, end_stage, q_score, blast_id, e_value, min_num_repeats, Num_reads_per_file, Coverage, Protocol)
	
	print("END OF RUN PROJECT")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-N", "--num_of_N", type=int, help="number of N's to add for merge of R1 and R2 pair-end reads", required=False, default=60) 
	parser.add_argument("-i", "--input_dir", type=str, help="a path to an input directory containing sequencing samples for pipeline analysis", required=True)
	parser.add_argument("-n", "--num_reads", type=int, help="number of reads per split file", required=False, default=25000) 
	parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory", required=True)
	parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference fasta file", required=True)
	parser.add_argument("-m", "--blast_mode", type=str, help="mode for blast, for Seq to Ref blast type SR, sr or SeqtoRef, for Ref to Seq blast type RS, rs or ReftoSeq, default = 'SeqtoRef'", required=False, default="SeqtoRef")
	parser.add_argument("-id", "--blast_id", type=int, help="% blast id, default=85", required=False, default=85)
	parser.add_argument("-t", "--blast_task", type=str, help="task for blast, blastn/megablast/dc-megablast?, default='blastn'", required=False, default="blastn")
	parser.add_argument("-ev", "--evalue", type=float, help="E-value for blast, default=1e-7", required=False, default=1e-7)
	parser.add_argument("-s", "--start", type=int, help="start step number", required=False)
	parser.add_argument("-e", "--end", type=int, help="end step number, default=6", required=False, default=6)
	parser.add_argument("-q", "--q_score", type=int, help="Q-score cutoff, default=30", required=False, default=30)
	parser.add_argument("-rep", "--repeats", type=int, help="number of repeats, default=2", required=False, default=2)
	parser.add_argument("-c", "--coverage", type=int, help="coverage cut-off for statistics, default=10000", required=False, default=10000)
	parser.add_argument("-pr", "--protocol", type=str, help="Library prep protocol is linear = 'L', 'l' or 'linear', or circular = 'C', 'c' or 'circular'. Default='linear'", required=False, default="linear")
	args = parser.parse_args()
	main(args)