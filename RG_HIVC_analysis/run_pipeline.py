from pbs_runners import script_runner, array_script_runner

# from RG_HIVC_analysis.constants import orig_samples

orig_samples= ['100888_S14','130945_S2','504185_S29','504190_S34','504194_S38','504198_S42','504202_S46','504206_S50','504210_S54','504215_S59','504221_S65','504226_S70','504231_S75','79504_S23','83476_S27','84864_S47','TASPX100926_S26','X145364_S75','X83354_S92','105090_S50','504181_S25','504186_S30','504191_S35','504195_S39','504199_S43','504203_S47','504207_S51','504211_S55','504217_S61','504223_S67','504227_S71','504233_S77','81065_S64','84298_S2','87415_S75','TASPX119494_S74','X145656_S76','X84335_S20','105094_S45','504182_S26','504188_S32','504192_S36','504196_S40','504200_S44','504204_S48','504208_S52','504212_S56','504218_S62','504224_S68','504228_S72','504234_S78','83351_S15','84434_S32','TASPX100702_S59','X100748_S73','X160138_S81','X84434_S3','105257_S39','504184_S28','504189_S33','504193_S37','504197_S41','504201_S45','504205_S49','504209_S53','504214_S58','504220_S64','504225_S69','504230_S74','504235_S79','83456_S87','84785_S56','TASPX100711_S29','X145364-R_S95','X83322_S89','X84994_S90']
control_samples = ['78292_S1', 'X105350_S13', 'X108456_S11', 'X117807_S26', 'X122087_S18', 'X128277_S27', 'X135623_S20', 'X161610_S28', 'X165309_S33', 'X83744_S2', 'X96026_S5', 'X101530_S9', 'X105354_S14', 'X111406_S10', 'X122054_S15', 'X122107_S17', 'X130823_S22', 'X160392_S24', 'X165276_S32', 'X178870_S29', 'X87707_S25', 'X96051_S7', 'X102662_S8', 'X108262_S21', 'X111489_S12', 'X122061_S16', 'X122107_S27', 'X132700_S84', 'X160433_S25', 'X165290_S34', 'X83348_S4', 'X88469_S6', '78486_S9', '79649_S49', '81088_S63', 'X105130_S8', 'X105416_S169', 'X112003_S13', 'X117597_S55', 'X117659_S88', 'X119113_S320', 'X122207_S32', 'X160394_S25', 'X84355_S11', 'X84505_S40', 'X84772_S12', 'X88273_S273', 'X94731_S49', '79578', '81003_S48', 'X101706', 'X105289_S56', 'X105500_S302', 'X117547_S16', 'X117607_S11', 'X117890_S46', 'X121017_S309', 'X135606_S69', 'X83676_S242', 'X84469_S265', 'X84641_S47', 'X84990_S248', 'X88373_S267', 'X96465_S54']


def make_ref_from_con():
    samples = orig_samples
    freq_files_dir = "/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/runs/high_qual/iter3/"
    prev_refs_dir  = "/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/runs/high_qual/iter2/cons/"
    new_cons_dir =  freq_files_dir +  "cons/"

    min_coverage_for_consensus = 500

    for sample in samples:
        sample_short = sample.split('_')[0]
        freqs_filename = sample_short + '.freqs'
        sample_dir_name = sample + '/'
        cons_filename = sample + ".fasta"

        input_freq_per_sample = freq_files_dir + sample_dir_name + freqs_filename
        prev_ref_per_sample = prev_refs_dir + cons_filename
        output_con_per_sample = new_cons_dir + cons_filename


        print('running make_ref_from_con.py, sample: {sample}'.format(sample= sample))
        script_runner('mkdir -p {new_cons_dir}; '
                      'python /sternadi/home/volume3/omer/SternLab/NGS_analysis/make_reference_from_consensus.py '
                      '-c {min_coverage_for_consensus} '
                      '-f {prev_ref_per_sample} '
                      '-p {input_freq_per_sample} '
                      '-o {output_con_per_sample}'.format(new_cons_dir= new_cons_dir, min_coverage_for_consensus= min_coverage_for_consensus, prev_ref_per_sample= prev_ref_per_sample, input_freq_per_sample=input_freq_per_sample, output_con_per_sample= output_con_per_sample),
                      alias='mrfc_{}'.format(sample),
                      load_python=True)


def run_pipeline_seperated_jobs():
    inputs_dir = "/sternadi/datasets/volume2/HIV_ravi_gupta_merged/"
    # inputs_dir = "/sternadi/datasets/volume2/HIV_ravi_gupta_control_merged/" #control

    samples = orig_samples
    refs_dir =   "/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/runs/high_qual/iter2/cons/"
    output_dir = "/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/runs/high_qual/iter3/"
    qual = 2 # 2-high, 1- low

    for sample in samples:
        input_dir_per_sample = inputs_dir + sample
        output_dir_per_sample = output_dir + sample
        ref = refs_dir + sample + ".fasta"

        print('running pipeline, sample: {sample}'.format(sample= sample))
        script_runner('mkdir -p {output_dir_per_sample}; '
                      'python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py '
                      '-i {input_dir_per_sample} '
                      '-o {output_dir_per_sample} '
                      '-r {ref} '
                      '-NGS_or_Cirseq 1 '
                      '-q 30 '
                      '-rep {qual} '
                      '-b 40 '
                      '-t z '
                      '-ev 1e-2'.format(input_dir_per_sample= input_dir_per_sample, output_dir_per_sample= output_dir_per_sample, ref= ref, qual=qual),
                      alias='pipeline_{}'.format(sample),
                      load_python=True)

def run_pipeline_batch():
    inputs_dir = "/sternadi/datasets/volume2/HIV_ravi_gupta_merged/"
    # inputs_dir = "/sternadi/datasets/volume2/HIV_ravi_gupta_control_merged/" #control

    refs_dir =   "/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/runs/high_qual/iter2/cons/"
    output_dir = "/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/runs/high_qual/iter3/"
    qual = 2 # 2-high, 1- low

    input_dir_per_sample = inputs_dir + '$sample'
    output_dir_per_sample = output_dir + '$sample'
    ref = refs_dir + '$sample' + ".fasta"


    array_script_runner('declare -a SAMPLES_ARRAY = (100888_S14  130945_S2  504185_S29 504190_S34  504194_S38  504198_S42 504202_S46  504206_S50  504210_S54  504215_S59  504221_S65  504226_S70  504231_S75  79504_S23  83476_S27  84864_S47  TASPX100926_S26  X145364_S75  X83354_S92 105090_S50  504181_S25  504186_S30  504191_S35  504195_S39  504199_S43  504203_S47  504207_S51  504211_S55  504217_S61  504223_S67  504227_S71  504233_S77  81065_S64  84298_S2   87415_S75  TASPX119494_S74  X145656_S76  X84335_S20 105094_S45  504182_S26  504188_S32  504192_S36  504196_S40  504200_S44  504204_S48  504208_S52  504212_S56  504218_S62  504224_S68  504228_S72  504234_S78  83351_S15  84434_S32  TASPX100702_S59  X100748_S73  X160138_S81  X84434_S3 105257_S39  504184_S28  504189_S33  504193_S37  504197_S41  504201_S45  504205_S49  504209_S53  504214_S58  504220_S64  504225_S69  504230_S74  504235_S79  83456_S87 84785_S56  TASPX100711_S29  X145364-R_S95 X83322_S89   X84994_S90); '
                        'sample=${{SAMPLES_ARRAY[$PBS_ARRAY_INDEX-1]}}; '
                        'mkdir -p {output_dir_per_sample}; '
                        'python /sternadi/home/volume1/shared/SternLab/pipeline_runner.py '
                        '-i {input_dir_per_sample} '
                        '-o {output_dir_per_sample} '
                        '-r {ref} '
                        '-NGS_or_Cirseq 1 '
                        '-q 30 '
                        '-rep {qual} '
                        '-b 40 '
                        '-t z '
                        '-ev 1e-2'.format(input_dir_per_sample= input_dir_per_sample, output_dir_per_sample= output_dir_per_sample, ref= ref, qual=qual),
                        jnum=38,
                        alias='pl_orig_high_iter3',
                        load_python=True)

if __name__ == "__main__":
    # make_ref_from_con()
    run_pipeline_batch()



