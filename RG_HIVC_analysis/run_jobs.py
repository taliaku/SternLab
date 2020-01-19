# from RG_HIVC_analysis.constants import orig_high_quality_patients
from pbs_runners import script_runner


## job submission without PBS files

# script_runner example
def run_merge_script_runner():
    original_samples_leftovers = ['105094_S45','105257_S39','130945_S2','504181_S25','504182_S26','504184_S28','504185_S29','504186_S30','504188_S32','504189_S33','504190_S34','504191_S35','504192_S36','504193_S37','504194_S38','504195_S39','504196_S40','504197_S41','504198_S42','504199_S43','504200_S44','504201_S45','504202_S46','504203_S47','504204_S48','504205_S49','504206_S50','504207_S51','504208_S52','504209_S53','504210_S54','504211_S55','504212_S56','504214_S58','504215_S59','504217_S61','504218_S62','504220_S64','504221_S65','504223_S67','504224_S68','504225_S69','504226_S70','504227_S71','504228_S72','504230_S74','504231_S75','504233_S77','504234_S78','504235_S79','79504_S23']
    original_input = '/sternadi/datasets/volume2/HIV_ravi_gupta_processed'
    original_output = '/sternadi/datasets/volume2/HIV_ravi_gupta_merged'

    fl_samples = ['78292_S1', 'X101530_S9', 'X102662_S8', 'X105350_S13', 'X105354_S14', 'X108262_S21', 'X108456_S11', 'X111406_S10', 'X111489_S12', 'X117807_S26', 'X122054_S15', 'X122061_S16', 'X122087_S18', 'X122107_S17', 'X122107_S27', 'X128277_S27', 'X130823_S22', 'X132700_S84', 'X135623_S20', 'X160392_S24', 'X160433_S25', 'X161610_S28', 'X165276_S32', 'X165290_S34', 'X165309_S33', 'X178870_S29', 'X83348_S4', 'X83744_S2', 'X87707_S25', 'X88469_S6', 'X96026_S5', 'X96051_S7']
    fl_input = '/sternadi/datasets/volume2/HIVC_Ravi_Gupta_control/1st_line_VF'
    fl_output = '/sternadi/datasets/volume2/HIV_ravi_gupta_CONTROL_merged/1st_line_VF'

    baseline_samples = ['78486_S9', '79578', '79649_S49', '81003_S48', '81088_S63', 'X101706', 'X105130_S8', 'X105289_S56', 'X105416_S169', 'X105500_S302', 'X112003_S13', 'X117547_S16', 'X117597_S55', 'X117607_S11', 'X117659_S88', 'X117890_S46', 'X119113_S320', 'X121017_S309', 'X122207_S32', 'X135606_S69', 'X160394_S25', 'X83676_S242', 'X84355_S11', 'X84469_S265', 'X84505_S40', 'X84641_S47', 'X84772_S12', 'X84990_S248', 'X88273_S273', 'X88373_S267', 'X94731_S49', 'X96465_S54']
    baseline_input = '/sternadi/datasets/volume2/HIVC_Ravi_Gupta_control/Baseline'
    baseline_output = '/sternadi/datasets/volume2/HIV_ravi_gupta_CONTROL_merged/Baseline'


    for sample in baseline_samples:
        print('sample is {sample}, input: {input}, output: {output}'.format(sample= sample, input= baseline_input, output= baseline_output))
        script_runner('mkdir -p {output}/{sample}; '
                      'python /sternadi/home/volume3/omer/SternLab/scripts/merge_fastq_files.py '
                      # '-f {input}/{sample}/{sample}_R1.fastq.gz '
                      # '-e {input}/{sample}/{sample}_R2.fastq.gz '
                      '-f {input}/{sample}/{sample}_L001_R1_001.fastq.gz '
                      '-e {input}/{sample}/{sample}_L001_R2_001.fastq.gz '
                      '-o {output}/{sample}/{sample}.merged.fastq.gz '
                      '-r 60'.format(sample= sample, input= baseline_input, output= baseline_output),
                      alias='merge_paired_end_{}'.format(sample),
                      load_python=True)


# array_script_runner
# def run_array_script_runner():
#     freqs_filename= sample.split('_')[0] + '.freqs'
#     output_filename= sample + '.fasta'
#     array_script_runner('python /sternadi/home/volume3/omer/SternLab/NGS_analysis/make_reference_from_consensus.py '
#                         '-f sternadi/home/volume1/shared/data/ref_genomes/HXB2.fasta '
#                         '-p /sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/run2/$sample/$freqs_filename '
#                         '-o /sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/refs/$output_filename '
#                         '-i $PBS_ARRAY_INDEX', jnum=76, alias='make_ref_from_con_{}'.format(sample), load_python=True)


if __name__ == "__main__":
    run_merge_script_runner()



