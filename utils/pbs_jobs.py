#! /usr/local/python_anaconda/bin/python3.4

import os
import subprocess
from time import sleep
import getpass
import datetime

USER_FOLDER_DICT = {"taliakustin": "/sternadi/home/volume1/taliakustin/temp",
                  "daniellem1":"/sternadi/home/volume1/daniellem1/temp",
                  "okushnir": "/sternadi/home/volume3/okushnir/running",
                  "omertirosh": "/sternadi/home/volume3/omer/logs",
                  'noamharel':'/sternadi/home/volume2/noam/logs',
                  'ita': '/sternadi/home/volume3/ita/logs',
                  'farage': '/sternadi/home/volume3/farage/temp'}


def create_pbs_cmd(cmdfile, alias, queue="adistzachi@power9", gmem=2, ncpus=1, ngpus=1, cmds="", dir = "",
                   load_python=True, jnum=False, run_after_job=None, nodes=1):
    with open(cmdfile, 'w') as o:
        o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
        o.write("#PBS -q %s\n" % queue)
        o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
        o.write("#PBS -N "+alias+"\n")
        #TODO: did dropping this "if" do something horrible?
        #if job_name in cmdfile and datetime.datetime.today().strftime('%Y-%m') in cmdfile:
        o.write("#PBS -o %s\n" % "/".join(cmdfile.split("/")[:-1]))
        o.write("#PBS -e %s\n" % "/".join(cmdfile.split("/")[:-1]))

        # running on GPUs 
        if queue == 'gpu':
            o.write(f"#PBS -l select=ngpus={ngpus}\n")
        else:    
            o.write(f"#PBS -l select={nodes}:ncpus={ncpus}:mem={gmem}gb\n")

        if jnum:
           if jnum != 1:
               o.write("#PBS -J 1-"+str(jnum)+"\n\n")
        if run_after_job != None and 'adi' in queue:
            o.write("#PBS -W depend=afterok:" + str(run_after_job)+ ".power9.tau.ac.il\n\n")
        if run_after_job != None and 'dudu' in queue:
            o.write("#PBS -W depend=afterok:" + str(run_after_job)+ ".power8.tau.ac.il\n\n")
    
        if dir != "":
            o.write("ls -land %s\n" % dir)
        o.write("id\n")
        o.write("date\n")
        o.write("hostname\n")
        if load_python:
            o.write("module load python/python-anaconda3.2019.10\n")
        
        o.write("\n")
        o.write('echo "%s"' % cmds)
        o.write("\n")
        o.write(cmds)
        o.write("\n")
        o.write("date\n")
    o.close()


def create_array_pbs_cmd(cmdfile, jnum, alias, gmem=7, cmds="", dir="", load_python=False, queue="adis", run_after_job=None):
    with open(cmdfile, 'w') as o:
        o.write("#!/bin/bash\n#PBS -S /bin/bash\n#PBS -j oe\n#PBS -r y\n")
        o.write("#PBS -q %s\n" % queue)
        o.write("#PBS -v PBS_O_SHELL=bash,PBS_ENVIRONMENT=PBS_BATCH \n")
        o.write("#PBS -N " + alias + "\n")
        if alias in cmdfile and datetime.datetime.today().strftime('%Y-%m') in cmdfile:
            o.write("#PBS -o %s\n" % "/".join(cmdfile.split("/")[:-1]))
        if (gmem):
            mem = gmem * 1000
            o.write("#PBS -l mem=" + str(mem) + "mb\n")
        if jnum == 1:
            raise Exception('Jnum parameter should be larger then 1\n')
        else:
            o.write("#PBS -J 1-" + str(jnum) + "\n\n")
            # #o.write("#PBS -J 3-4 \n")
        if run_after_job != None:
            o.write("#PBS -W depend=afterok:" + str(run_after_job)+ ".power9.tau.ac.il\n\n")
        if dir != "":
            o.write("ls -land %s\n" % dir)
        o.write("id\n")
        o.write("date\n")
        o.write("hostname\n")
        if load_python:
            o.write("module load python/python-anaconda3.2019.7\n")

        o.write("\n")
        o.write(cmds)
        o.write("\n")
        o.write("date\n")
        o.write("echo DONE\n")
    o.close()


def submit(cmdfile):
    cmd = "/opt/pbs/bin/qsub " + cmdfile
    result = os.popen(cmd).read()
    return result.split(".")[0]


def check_pbs(job_id):
    """
    :param job_id: The PBS job id
    :return: "Done!", when the job is done
    """
    status = "Running..."
    try:
        process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
        while process != "":
            process = subprocess.check_output("qstat | grep " + str(job_id), shell=True)
            sleep(0.05)
        print("")
    except (subprocess.CalledProcessError):
        process = ""
    if process == "":
        status = "Done"
    return status


def assign_cmdfile_path(cmdname, alias):
    username = getpass.getuser()
    if username in USER_FOLDER_DICT.keys():
        tmp_dir = USER_FOLDER_DICT[username]
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        date = datetime.datetime.today().strftime('%Y-%m')
        tmp_dir = tmp_dir + "/%s" % date
        if not os.path.exists(tmp_dir):
            os.system("mkdir %s" % tmp_dir)
        tmp_dir = tmp_dir + "/%s" % alias
        if not os.path.exists(tmp_dir):
            os.mkdir(tmp_dir)
        cmdname = tmp_dir + "/" + cmdname
    return cmdname


def check_killed_daily():
    username = getpass.getuser()
    if username not in USER_FOLDER_DICT.keys():
        return
    tmp_dir = USER_FOLDER_DICT[username]
    if not os.path.exists(tmp_dir):
        return
    date = datetime.datetime.today().strftime('%Y-%m')
    tmp_dir = os.path.join(tmp_dir, date)
    os.chdir(tmp_dir)
    process = subprocess.getoutput('find ./* -maxdepth 1 -mtime -1 | grep OU | xargs grep --include "*OU" kill')    # get all files that had SIG KILL that day
    return process



