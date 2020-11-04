
import os,sys,inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0,parentdir)
from utils import pbs_jobs
import os
from file_utilities import check_filename, check_dirname, make_dir
from pbs_runners import mafft_runner, phyml_runner, script_runner
import argparse

OUTPUT_DIR = '/sternadi/nobackup/volume1/covid/israel_artic_pipeline/'
UNITED_OUTPUT_DIR = '/sternadi/nobackup/volume1/covid/israel_artic_pipeline/all_results/'
REFERENCE_FILE = '/sternadi/home/volume2/noam/covid/MN908947.fasta'


def covid_artic_pipeline(input_dir, output_dir, alias="ARTIC_pipeline", queue="adistzachi"):
    input_dir = check_dirname(input_dir + '/Raw_data/')
    output_dir = check_dirname(output_dir, Truedir = False)
    # create output dirs
    make_dir(output_dir)
    make_dir(f"{output_dir}/ptrimmer_cleanup/")
    make_dir(f"{output_dir}/python_pipeline/")

    cmdfile = pbs_jobs.assign_cmdfile_path('ARTIC_pipeline_p1', alias);
    gmem = 10;

    # ptrimmer
    cmds = f"for d in {input_dir}/*R1*; " \
           f"do dir_name=$(basename $d); " \
           f"dir_name=$(echo $dir_name | cut -d'_' -f1); " \
           f"mkdir {output_dir}/ptrimmer_cleanup/$dir_name; " \
           f"echo $d >> {output_dir}/ptrimmer_cleanup/ptrimmer.log; " \
           f'/sternadi/home/volume1/shared/tools/pTrimmer/pTrimmer-1.3.1 -s pair -a /sternadi/home/volume2/noam/covid/artic_amplicons/ptrimmer_primers.txt --read1 $d --read2 "${{d//_R1_/_R2_}}" -o {output_dir}/ptrimmer_cleanup/$dir_name &>> {output_dir}/ptrimmer_cleanup/ptrimmer.log; ' \
           f"done\n"
    # organize file names for pipeline
    cmds += f'for file in {output_dir}/ptrimmer_cleanup/*/*.fq; do mv "$file" "${{file%.*}}.fastq"; done\n'
    cmds += f'for file in {output_dir}/ptrimmer_cleanup/*; do mv "$file" "${{file%.*}}_L001"; done\n'
    # run pipeline
    cmds += f'python /sternadi/home/volume2/noam/SternLab/Python_pipeline/Project_Runner.py -o {output_dir}/python_pipeline -i {output_dir}/ptrimmer_cleanup/ -r {REFERENCE_FILE} -x 1 -c 0 -v 1e-9'
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, queue=queue, gmem=gmem, cmds=cmds)
    job_id_p1 = pbs_jobs.submit(cmdfile)

    # after pipelines are finished: create concensuses, realign together with Wuhan reference and create tree and mutation csv, create_big_freqs_file,
    cmdfile = pbs_jobs.assign_cmdfile_path('ARTIC_pipeline_p2', alias)
    # create concensuses
    make_dir(f"{output_dir}/python_pipeline/freqs/")
    cmds = f"cp {output_dir}/python_pipeline/*/*.freqs.csv {output_dir}/python_pipeline/freqs/\n"
    cmds += f"python /sternadi/home/volume2/noam/SternLab/covid19/get_consensus_covid19.py -d {output_dir}/python_pipeline/freqs/ -o $(basename {output_dir})\n"
    # concat together with Wuhan reference
    cmds += f"cat /sternadi/home/volume2/noam/covid/MN908947.fasta {output_dir}/python_pipeline/freqs/$(basename {output_dir})_consensus_all.fasta > {output_dir}/consensuses_w_wuhan.fasta\n"
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias+'_p2', queue=queue, gmem=gmem, cmds=cmds, run_after_job=job_id_p1)
    job_id_p2 = pbs_jobs.submit(cmdfile)
    # align and create tree
    os.system(f'touch {output_dir}/consensuses_w_wuhan.fasta {output_dir}/consensuses_w_wuhan.aln {output_dir}/consensuses_w_wuhan.aln.phy')
    job_id_p3 = mafft_runner(f"{output_dir}/consensuses_w_wuhan.fasta", run_after_job=job_id_p2)
    job_id_p4 = script_runner(f"python /sternadi/home/volume2/noam/SternLab/scripts/fasta_to_phylip.py -i {output_dir}/consensuses_w_wuhan.aln -o {output_dir}/consensuses_w_wuhan.aln.phy", run_after_job=job_id_p3)
    job_id_p5 = phyml_runner(f"{output_dir}/consensuses_w_wuhan.aln.phy", run_after_job=job_id_p4)
    # create mutations csv
    job_id_p6 = script_runner(f"python /sternadi/home/volume2/noam/SternLab/scripts/compare_aligned_msa.py -i {output_dir}/consensuses_w_wuhan.aln -r 'MN908947.3' -o {output_dir}/consensuses_w_wuhan.mutations.csv -e y", alias='mutations_from_msa', run_after_job=job_id_p3)
    # create big freqs file
    job_id_p7 = script_runner(f"python /sternadi/home/volume2/noam/SternLab/scripts/unite_all_freq_files.py -d '{output_dir}/python_pipeline/freqs/' -o '{output_dir}/python_pipeline/freqs/all_freqs.csv'", alias='unite_freqs', run_after_job=job_id_p3)
    return (job_id_p1, job_id_p2, job_id_p3, job_id_p4, job_id_p5, job_id_p6, job_id_p7)

def add_to_previous_results(run_after_job, new_dir, all_results_dir=UNITED_OUTPUT_DIR):
    # unite concensuses
    job_id_p1 = script_runner(f'cat {REFERENCE_FILE} {all_results_dir}/../*/python_pipeline/freqs/*consensus_all.fasta > {all_results_dir}israel_sequences.fasta', run_after_job=run_after_job)
    # align and make tree
    job_id_p2 = mafft_runner(f"{all_results_dir}israel_sequences.fasta", run_after_job=job_id_p1)
    job_id_p3 = script_runner(f"python /sternadi/home/volume2/noam/SternLab/scripts/fasta_to_phylip.py -i {all_results_dir}israel_sequences.aln -o {all_results_dir}israel_sequences.aln.phy", run_after_job=job_id_p2)
    job_id_p4 = phyml_runner(f"{all_results_dir}israel_sequences.aln.phy", run_after_job=job_id_p3)
    #### TO DO - color tree according to directories
    # mutations csv
    job_id_p5 = script_runner(f"python /sternadi/home/volume2/noam/SternLab/scripts/compare_aligned_msa.py -i {all_results_dir}israel_sequences.aln -r 'MN908947.3' -o {all_results_dir}israel_sequences.mutations.csv -e y", alias='mutations_from_msa', run_after_job=job_id_p2)
    # concat freqs
    job_id_p6 = script_runner(f'python /sternadi/home/volume2/noam/SternLab/scripts/concat_dfs.py {all_results_dir}/../*/python_pipeline/freqs/all_freqs.csv -o {all_results_dir}israel_freqs.csv', run_after_job=run_after_job)
    return (job_id_p1, job_id_p2, job_id_p3, job_id_p4, job_id_p5, job_id_p6)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", type=str, help="input directory as downloaded from Technion, contains gzipped fastqs",
                        required=True)
    args = parser.parse_args()
    if not vars(args):
        parser.print_help()
        parser.exit(1)
    
    input_dir = check_dirname(args.input_dir)
    if input_dir.endswith('/'):
        input_dir = input_dir[:-1]
    output_dir = OUTPUT_DIR + os.path.basename(input_dir)
    make_dir(check_dirname(output_dir, Truedir=False))
    print(f'reading from: {output_dir}')
    print(f'writing to : {input_dir}')
    # run analysis on new data
    job_ids = covid_artic_pipeline(input_dir, output_dir)
    # add to all israel analysis files
    add_to_previous_results(job_ids[-1], output_dir)