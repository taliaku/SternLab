#! python/python-anaconda3.2019.7

import argparse
import time
import glob
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))) 
from utils.runner_utils import FindFilesInDir, check_queue, create_pbs_cmd, submit, Sleep, create_array		

def run_project(pipeline_path, input_dir, dir_path, ref_genome, mode, task, start_stage, end_stage, q_score, blast_id,
                e_value, min_num_repeats, Num_reads_per_file, Coverage, Protocol, queue, overwrite, sample_basename_pattern):
    alias = "RunProject"

    file_type = sample_basename_pattern + "*"
    project_samples = FindFilesInDir(input_dir, file_type)
    if len(project_samples) == 0:
        raise Exception("Unexpected error, input directory " + input_dir + " is empty\n")

    file_type = os.path.basename(project_samples[0]).split('L001')[0] + "*" # pattern extraction logic- 'pattern_L00*_R*.fastq'
    Lanes = len(FindFilesInDir(input_dir, file_type))

    samples_list = []
    for sample in project_samples:
        if not '_L00' in sample:
            raise Exception("Unexpected error. '_L00' is missing in sample name, cannot identify sample name pattern\n")
        else:
            sample_name = os.path.basename(sample).split('_L00')[0]

        sample_pattern = sample + "/" + sample_name + "_"
        if "FASTQ_Generation_" in sample:
            sample_pattern = sample + "/" + sample_name.replace("_","-") + "_"

        if not sample_pattern in samples_list:
            sample_dir_path = dir_path + "/" + sample_name
            if not os.path.isdir(sample_dir_path):
                try:
                    os.system("mkdir -p " + sample_dir_path)
                except:
                    raise Exception("failed to create input directory " + sample_dir_path + "\n")
                if not os.path.isdir(sample_dir_path):
                    raise Exception("Sample directory " + sample_dir_path + " does not exist or is not a valid directory path\n")

            if Lanes == 1:
               samples_list.extend((sample_pattern, sample_dir_path))
            elif Lanes > 1:
                input_samples_dir = sample_dir_path + "/input_files"
                try:
                    os.system("mkdir -p " + input_samples_dir)
                except:
                    raise Exception("failed to create input directory " + input_samples_dir + "\n")
                if not os.path.isdir(input_samples_dir):
                    raise Exception("Sample directory " + input_samples_dir + " does not exist or is not a valid directory path\n")

        if Lanes > 1:
            file_type = "*_R*"
            Sample_reads = FindFilesInDir(sample, file_type)
            for read in Sample_reads:
                os.system("cp " + read + " " + input_samples_dir + "/" + os.path.basename(read))
            sample_pattern = input_samples_dir + "/"
            if not sample_pattern in samples_list:
                samples_list.extend((sample_pattern, sample_dir_path))

    if Lanes > 1:
        for i in range(0, len(samples_list)-1, 2):
            file_type = "*_R*"
            Sample_reads = FindFilesInDir(os.path.dirname(samples_list[i]), file_type)
            sample_pattern = os.path.commonprefix(Sample_reads)
            samples_list[i] = sample_pattern

    if len(samples_list) % 2 != 0:
        raise Exception("Unexpected error, number of samples for pipeline does not match number of output directory created\n")

    file_type = "*_R2*"
    R2_samples_files = len(FindFilesInDir(project_samples[0], file_type))
    if start_stage == None:
        if R2_samples_files > 0:  # paired-end reads
            start_stage = 0
        else:  # single-end reads
            start_stage = 1

    if start_stage > end_stage:
        raise Exception("Unexpected error, start stage " + str(start_stage) + " is larger than end stage " + str(end_stage) + "\n")

    num_of_samples = int(len(samples_list)/2)
    if num_of_samples == 1:
        p = '0'
        o = '1'
        gmem = 2
    else:
        p = '$((PBS_ARRAY_INDEX*2))-2'
        o = '$((PBS_ARRAY_INDEX*2))-1'
        gmem = 7

    array = create_array(samples_list)
    pipeline_dir = pipeline_path[0:pipeline_path.rfind('/')]
    cmd1 = 'declare -a SAMPLENAMES\n'
    cmd2 = 'SAMPLENAMES=' + array + "\n\n"
    cmd3 = "python " + pipeline_path + " -i ${SAMPLENAMES[" + p + "]} -o ${SAMPLENAMES[" + o + "]} -r " + ref_genome + " -m " + mode + " -t " + task + " -s " + str(start_stage) + \
           " -e " + str(end_stage) + " -q " + str(q_score) + " -d " + str(blast_id) + " -v " + str(e_value) + " -x " + str(min_num_repeats) + " -n " + str(Num_reads_per_file) + \
           " -c " + str(Coverage) + " -p " + Protocol + " -u " + queue + " -w " + overwrite + "\n"
    cmd4 = f"python {pipeline_dir}/AggregateSummaries.py -i {dir_path} -o {dir_path}/AggregatedSummary.csv"
    cmds = cmd1 + cmd2 + cmd3 + cmd4
    cmdfile = dir_path + "/pipeline_project_runner.cmd"
    create_pbs_cmd(cmdfile=cmdfile, alias=alias, jnum=num_of_samples, gmem=gmem, cmds=cmds, queue=queue, load_python=True)
    job_id = submit(cmdfile)
    print(job_id)
    Sleep(alias, job_id)

def main(args):
    pipeline_dir = os.path.dirname(os.path.abspath(__file__).strip())
    pipeline_path = pipeline_dir + "/Runner.py"

    if not (os.path.isfile(pipeline_path) and os.path.splitext(pipeline_path)[1] == '.py'):
        raise Exception("Unexpected error, pipeline path " + pipeline_path + " does not exist or not a .py file\n")

    start_stage = args.start
    if start_stage != None:
        if not start_stage in [0, 1, 2, 3, 4, 5, 6]:
            raise Exception("Unexpected error, start_stage " + str(start_stage) + " is not a valid value\n")

    end_stage = args.end
    if end_stage != None:
        if not end_stage in [0, 1, 2, 3, 4, 5, 6]:
            raise Exception("Unexpected error, end_stage is not a valid integer value between 0-6\n")

    if start_stage != None and end_stage != None and start_stage > end_stage:
        raise Exception("Unexpected error, start stage " + str(start_stage) + " is larger than end stage " + str(end_stage) + "\n")

    input_dir = os.path.dirname(args.input_dir.strip())
    sample_basename_pattern = os.path.basename(args.input_dir.strip())
    if not os.path.isdir(input_dir):
        raise Exception("Directory for input files " + input_dir + " does not exist or is not a valid directory path\n")
    if sample_basename_pattern != "":
        sample_basename_pattern = sample_basename_pattern + "*"

    number_of_N = args.num_of_N
    if number_of_N != None:
        if number_of_N < 60:
            print("\nWarning. Running merge reads with number_of_N smaller than 60\n")

    dir_path = args.output_dir.strip()
    if not os.path.isdir(dir_path):
        try:
            os.system("mkdir -p " + dir_path)
        except:
            raise Exception("failed to create input directory " + dir_path + "\n")
        if not os.path.isdir(dir_path):
            raise Exception("Directory " + dir_path + " does not exist or is not a valid directory path\n")

    ref_genome = args.ref.strip()
    if not (os.path.isfile(ref_genome) and os.path.splitext(ref_genome)[1] == '.fasta'):
        raise Exception("Unexpected error, reference genome file does not exist, is not a file or is not a fasta file\n")

    q_score = args.q_score
    if q_score != None:
        if q_score < 0 or q_score > 40:
            raise Exception("Unexpected error, q-score value " + str(q_score) + " is not valid, should be an integer value between 0-40\n")
        if q_score < 16:
            print("\nWarning, running pipeline with q-score value of " + str(q_score) + "\n")

    blast_id = args.blast_id
    if blast_id != None:
        if blast_id < 85:
            print("\nWarning, running pipeline with blast id value of " + str(blast_id) + "\n")
        if blast_id < 0 or blast_id > 100:
            raise Exception("Unexpected error, identity % for blast is not a valid value: " + str(blast_id) + " \n")

    e_value = args.evalue
    if e_value != None:
        if e_value > 1e-7:
            print("\nWarning, running pipeline with e_value > " + str(e_value) + "\n")

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
            raise Exception("Unexpected error, min number of repeats is less than 1\n")
        if min_num_repeats < 2:
            print("\nWarning. Running pipeline with min number of repeats less than 2\n")
        if min_num_repeats > 2:
            print("\nWarning. Running pipeline with min number of repeats bigger than 2\n")

    Min_Num_reads_per_file = 10000
    Max_Num_reads_per_file = 40000
    Num_reads_per_file = args.num_reads
    if Num_reads_per_file != None:
        if Num_reads_per_file < Min_Num_reads_per_file:
            print("\nWarning. Running pipeline with less than " + str(Min_Num_reads_per_file) + " reads per split file\n")
        if Num_reads_per_file > Max_Num_reads_per_file:
            print("\nWarning. Running pipeline with more than " + str(Max_Num_reads_per_file) + " reads per split file\n")

    Min_Coverage = 1000
    Coverage = args.coverage
    if Coverage != None:
        if Coverage < Min_Coverage:
            print("\nWarning. Running pipeline with coverage smaller than " + str(Min_Coverage) + "\n")

    Protocol = args.protocol.strip()
    if Protocol not in ["L", "l", "linear", "C", "c", "circular"]:
        raise Exception("Unexpected error, for linear library prep protocol type 'linear' or 'L', for circular library prep protocol type 'circular' or 'C'\n")

    overwrite = args.overwrite.strip()
    if overwrite not in ["Y", "y", "N", "n"]:
        raise Exception("Unexpected error, overwrite value is different than Y/N\n")

    cmd = "python {} -i {} -o {} -r {} -m {} -t {} -s {} -e {} -q {} -d {} -v {} -x {} -n {} -c {} -p {} -u {} -w {}".format(
        pipeline_path, input_dir, dir_path, ref_genome, mode, task, start_stage, end_stage, q_score, blast_id, e_value, min_num_repeats, Num_reads_per_file, Coverage, Protocol, queue, overwrite)
    print(cmd)

    run_project(pipeline_path, input_dir, dir_path, ref_genome, mode, task, start_stage, end_stage, q_score, blast_id, e_value,
                min_num_repeats, Num_reads_per_file, Coverage, Protocol, queue, overwrite, sample_basename_pattern)

    print("END OF RUN PROJECT")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-N", "--num_of_N", type=int, help="number of N's to add for merge of R1 and R2 pair-end reads", required=False, default=60)
    parser.add_argument("-i", "--input_dir", type=str, help="a path to an input directory containing sequencing samples for pipeline analysis, example: dir1/dir2/. \
                        Sample pattern is allowed, example: dir1/dir2/sample_1_", required=True)
    parser.add_argument("-n", "--num_reads", type=int, help="number of reads per split file", required=False, default=25000)
    parser.add_argument("-o", "--output_dir", type=str, help="a path to an output directory", required=True)
    parser.add_argument("-r", "--ref", type=str, help="a path to a genome reference fasta file", required=True)
    parser.add_argument("-m", "--blast_mode", type=str, help="mode for blast, for Seq to Ref blast type SR, sr or SeqtoRef, for Ref to Seq blast type RS, rs or ReftoSeq, default = 'SeqtoRef'", required=False, default="SeqtoRef")
    parser.add_argument("-d", "--blast_id", type=int, help="% blast id, default=85", required=False, default=85)
    parser.add_argument("-t", "--blast_task", type=str, help="task for blast, blastn/megablast/dc-megablast?, default='blastn'", required=False, default="blastn")
    parser.add_argument("-v", "--evalue", type=float, help="E-value for blast, default=1e-7", required=False, default=1e-7)
    parser.add_argument("-s", "--start", type=int, help="start stage number", required=False)
    parser.add_argument("-e", "--end", type=int, help="end stage number, default=5", required=False, default=5)
    parser.add_argument("-q", "--q_score", type=int, help="Q-score cutoff, default=30", required=False, default=30)
    parser.add_argument("-x", "--repeats", type=int, help="number of repeats, default=2", required=False, default=2)
    parser.add_argument("-c", "--coverage", type=int, help="coverage cut-off for statistics, default=10000", required=False, default=10000)
    parser.add_argument("-u", "--queue", type=str, help="queue to run pipeline, default='tzachi@power9'", required=False, default="tzachi@power9")
    parser.add_argument("-p", "--protocol", type=str, help="Library prep protocol is linear = 'L', 'l' or 'linear', or circular = 'C', 'c' or 'circular'. Default='linear'",
                        required=False, default="linear")
    parser.add_argument("-w", "--overwrite", type=str, help="overwrite? Y/N, default='N'", required=False, default="N")

    args = parser.parse_args()
    main(args)