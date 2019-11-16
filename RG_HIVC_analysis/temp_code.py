from pbs_runners import script_runner, array_script_runner

# Job submission without PBS files:

# array_script_runner
# def run_array_script_runner():
#     freqs_filename= sample.split('_')[0] + '.freqs'
#     output_filename= sample + '.fasta'
#     array_script_runner('python /sternadi/home/volume3/omer/SternLab/NGS_analysis/make_reference_from_consensus.py '
#                         '-f sternadi/home/volume1/shared/data/ref_genomes/HXB2.fasta '
#                         '-p /sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run2/$sample/$freqs_filename '
#                         '-o /sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/refs/$output_filename '
#                         '-i $PBS_ARRAY_INDEX', jnum=76, alias='make_ref_from_con_{}'.format(sample), load_python=True)


# script_runner
def run_script_runner():
    path = '/sternadi/datasets/volume2/HIV_ravi_gupta_processed'
    # original_samples_leftovers = ['105090_S50','105094_S45','105257_S39','130945_S2','504181_S25','504182_S26','504184_S28','504185_S29','504186_S30','504188_S32','504189_S33','504190_S34','504191_S35','504192_S36','504193_S37','504194_S38','504195_S39','504196_S40','504197_S41','504198_S42','504199_S43','504200_S44','504201_S45','504202_S46','504203_S47','504204_S48','504205_S49','504206_S50','504207_S51','504208_S52','504209_S53','504210_S54','504211_S55','504212_S56','504214_S58','504215_S59','504217_S61','504218_S62','504220_S64','504221_S65','504223_S67','504224_S68','504225_S69','504226_S70','504227_S71','504228_S72','504230_S74','504231_S75','504233_S77','504234_S78','504235_S79','79504_S23']
    original_samples_leftovers = ['105090_S50']
    for sample in original_samples_leftovers:
        print('sample is {sample}, path is {path}'.format(sample= sample, path= path))
        script_runner('python /sternadi/home/volume3/omer/SternLab/scripts/merge_fastq_files.py '
                      '-f {path}/{sample}/{sample}_R1.fastq.gz '
                      '-e {path}/{sample}/{sample}_R2.fastq.gz '
                      '-o {path}/{sample}/{sample}.merged.fastq.gz '
                      '-r 60'.format(sample= sample, path= path),
                      alias='merge_paired_end_{}'.format(sample),
                      load_python=True)


if __name__ == "__main__":
    run_script_runner()

