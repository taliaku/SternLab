#! python/python-anaconda3.2019.7

import argparse
import time
import glob
import os

def create_pbs_cmd(cmdfile, alias, jnum, gmem, cmds, queue, load_python=True):
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
			o.write("module load python/anaconda_python-3.6.1\n")       
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

def run_multi_projects(pipeline_path, cmds_file, queue):
	alias = "RunMultiProjects"
	
	num_of_cmd = int(os.popen("awk 'END {print NR}' " + cmds_file).read())
	if num_of_cmd == 0:
		raise Exception ("Unexpected error, file " + cmds_file + " is empty, or number of lines in file does not divide by 2\n")
	elif num_of_cmd == 1:
		NR_value = '1'
		gmem = 2
	else:
		NR_value = '$PBS_ARRAY_INDEX'
		gmem = 7
	
	dir_path = os.path.dirname(cmds_file)
	cmdfile = dir_path + "/RunMultiProjects.cmd"	
	cmd1 = 'CMD=$(awk "NR==' + NR_value + '" ' + cmds_file + ')\n'
	cmd2 = "python " + pipeline_path + " $CMD "
	cmds = cmd1 + cmd2
	create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_cmd, gmem=gmem, cmds=cmds, queue=queue, load_python=True)
	job_id = submit(cmdfile)
	print(job_id)
	Sleep(alias, job_id)
	time.sleep(10)
	
def main(args):
	
	#pipeline_path = args.pipeline_runner
	pipeline_dir = os.path.dirname(os.path.abspath(__file__).strip())
	pipeline_path = pipeline_dir + "/Runner.py"

	if not os.path.isfile(pipeline_path):
		raise Exception("Unexpected error, " + pipeline_path + " does not exist, is not a file or or is not a blast file\n")
	
	cmds_file = args.cmds_file
	if not os.path.isfile(cmds_file):
		raise Exception("Unexpected error, " + cmds_file + " cmds file does not exist or is not a file\n")

	queue = args.queue
	if queue != None:
		if queue not in ["inf", "hugemem", "pup-interactive", "parallel", "adis", "adis-long"]:
			raise Exception("Unexpected error, queue has to be either 'inf', 'hugemem', 'pup-interactive', 'parallel', 'adis', 'adis-long'\n")

	run_multi_projects(pipeline_path, cmds_file, queue)
	
	print("END OF RUN MULTI PROJECTS")

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("-f", "--cmds_file", type=str, help="a path to a file containing a list of cmds to run. Different variables for each cmd are excepted", required=True)
	parser.add_argument("-u", "--queue", type=str, help="queue to run pipeline, default='adis'", required=False, default="adis")
	args = parser.parse_args()
	main(args)
