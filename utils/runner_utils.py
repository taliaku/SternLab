#! python/python-anaconda3.2019.7

import os
import sys
import glob
import time
from utils.logger import pipeline_logger

def check_queue(queue):
	allowed_queues = ["inf", "hugemem", "pup-interactive", "parallel", "adis", "adis-long", "tzachi@power9", 'adistzachi'] 
	if queue not in allowed_queues:
		raise Exception(f"Sorry but queue must be one of {allowed_queues}, not '{queue}'")

def create_pbs_cmd(cmdfile, alias, jnum, gmem, cmds, queue, load_python=True, run_after_job=None):
	with open(cmdfile, 'w') as o:
		o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
		o.write("#PBS -q %s\n" % queue)
		o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
		o.write("#PBS -N "+ alias+"\n")
		o.write("#PBS -o %s\n" % "/".join(cmdfile.split("/")[:-1]))
		o.write("#PBS -e %s\n" % "/".join(cmdfile.split("/")[:-1]))
		if run_after_job != None:
			o.write(f"#PBS -W depend=afterok: {run_after_job}\n\n")
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
	log = pipeline_logger()
	cmd = "/opt/pbs/bin/qsub " + cmdfile
	result = os.popen(cmd).read()
	if 'power' in result:
		return result.split(".")[0]
	else:
		log.error(f"{cmdfile} was not submitted")	

def Sleep (alias, job_id, sleep_max=1200000, sleep_quantum=10, queue='adistzachi'):
	log = pipeline_logger()
	log.info(f"Starting {alias} with job id: {job_id}")
	start_time = time.time()
	i = 0
	if job_id[-2:]=='[]': # grep doesn't like these..
		job_id = job_id[:-2]
	qstat_command = f"qstat @power9 | grep {job_id}" #TODO: do we really need to specify power9?
	qstat_result = os.popen(qstat_command).read()
	while len(qstat_result) > 0 and i <= sleep_max: 
		for second in range(0, sleep_quantum):
			elapsed_time = time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time))
			sys.stdout.write("\r")
			sys.stdout.write(f'Elapsed time: {elapsed_time}')
			sys.stdout.flush()
			time.sleep(1)
		i += sleep_quantum
		qstat_result = os.popen(qstat_command).read()
	if len(qstat_result) > 0: 
		raise Exception(alias + " stage was not completed. Max sleep time reached\n")
	sys.stdout.write("\n")
	log.info(f"{alias} Done.")

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