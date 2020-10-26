
import sys
sys.path.append("..")
from utils import pbs_jobs
import os
from file_utilities import check_filename, check_dirname, make_dir

def covid_artic_pipeline(input_dir, output_dir, alias="ARTIC_pipeline", queue="adistzachi"):
    input_dir = check_dirname(input_dir)
    output_dir = check_dirname(output_dir, Truedir = False)
    # create output dirs
    make_dir(output_dir)
    make_dir(f"{output_dir}/ptrimmer_cleanup/")
    make_dir(f"{output_dir}/python_pipeline_x1_c0_after_ptrimmer/")

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
    cmds += f'python /sternadi/home/volume2/noam/SternLab/Python_pipeline/Project_Runner.py -o {output_dir}/python_pipeline_x1_c0_after_ptrimmer -i {output_dir}/ptrimmer_cleanup/ -r /sternadi/home/volume2/noam/covid/MN908947.fasta -x 1 -c 0'
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias, queue=queue, gmem=gmem, cmds=cmds)
    job_id_p1 = pbs_jobs.submit(cmdfile)

    # after pipelines are finished: create concensuses, realign together with Wuhan reference and create tree and mutation csv, create_big_freqs_file,
    cmdfile = pbs_jobs.assign_cmdfile_path('ARTIC_pipeline_p2', alias);
    # create concensuses
    make_dir(f"{output_dir}/python_pipeline_x1_c0_after_ptrimmer/freqs/")
    cmds = f"cp {output_dir}/python_pipeline_x1_c0_after_ptrimmer/*/*.freqs.csv {output_dir}/python_pipeline_x1_c0_after_ptrimmer/freqs/\n"
    cmds += f"python /sternadi/home/volume2/noam/SternLab/covid19/get_consensus_covid19.py -d {output_dir}/python_pipeline_x1_c0_after_ptrimmer/freqs/ -o $(basename {output_dir})\n"
    # concat together with Wuhan reference
    cmds += f"cat /sternadi/home/volume2/noam/covid/MN908947.fasta {output_dir}/python_pipeline_x1_c0_after_ptrimmer/freqs/$(basename {output_dir})_consensus_all.fasta > {output_dir}/consensuses_w_wuhan.fasta\n"
    pbs_jobs.create_pbs_cmd(cmdfile, alias=alias+'_p2', queue=queue, gmem=gmem, cmds=cmds, run_after_job=job_id_p1)
    job_id_p2 = pbs_jobs.submit(cmdfile)
    # align and create tree
    os.system(f'touch {output_dir}/consensuses_w_wuhan.fasta {output_dir}/consensuses_w_wuhan.aln')
    job_id_p3 = mafft_runner(f"{output_dir}/consensuses_w_wuhan.fasta", run_after_job=job_id_p2)
    job_id_p4 = phyml_runner(f"{output_dir}/consensuses_w_wuhan.aln", phylip=False)# run_after_job=job_id_p3)
    # create mutations csv
    job_id_p5 = script_runner(f"python /sternadi/home/volume2/noam/SternLab/scripts/compare_aligned_msa.py -i {output_dir}/consensuses_w_wuhan.aln -r 'MN908947.3' -o {output_dir}/consensuses_w_wuhan.mutations.csv", alias='mutations_from_msa', run_after_job='job_id_p4')
    # create big freqs file
    job_id_p6 = script_runner(f"python /sternadi/home/volume2/noam/SternLab/scripts/unite_all_freq_files.py -d f'{output_dir}/python_pipeline_x1_c0_after_ptrimmer/freqs/' -o f'{output_dir}/python_pipeline_x1_c0_after_ptrimmer/freqs/all_freqs.csv'", alias='unite_freqs', run_after_job='job_id_p4')

    #### to do: unite for all libraries, add hard coded paths for where the files are

    return (job_id_p1, job_id_p2, job_id_p3, job_id_p4, job_id_p5, job_id_p6)