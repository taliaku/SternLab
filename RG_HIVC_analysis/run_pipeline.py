from RG_HIVC_analysis.constants import orig_samples

from pbs_runners import script_runner


def make_ref_from_con():
    samples = orig_samples
    min_coverage_for_consensus = 500

    run_dir = "/sternadi/home/volume1/shared/analysis/HIV_ravi_gupta/runs/high_qual/"
    currect_iteration_dir = run_dir + "iter2/"
    cons_output_dir =  currect_iteration_dir +  "cons/"

    for sample in samples:
        sample_short = sample.split('_')[0]
        freqs_filename = sample_short + '.freqs'
        cons_filename = sample + ".fasta"

        input_freq_fname = currect_iteration_dir + sample + '/' + freqs_filename
        output_cons_fname = cons_output_dir + cons_filename

        prev_iteration_dir = run_dir + "iter1/"
        prev_iter_ref_fname = prev_iteration_dir + 'cons/' + cons_filename

        print('running make_ref_from_con.py, sample: {sample}'.format(sample= sample))
        script_runner('mkdir -p {cons_output_dir}; '
                      'python /sternadi/home/volume3/omer/SternLab/NGS_analysis/make_reference_from_consensus.py '
                      '-c {min_coverage_for_consensus} '
                      '-f {prev_iter_ref_fname} '
                      '-p {input_freq_fname} '
                      '-o {output_cons_fname}'.format(cons_output_dir= cons_output_dir, min_coverage_for_consensus= min_coverage_for_consensus, prev_iter_ref_fname= prev_iter_ref_fname, input_freq_fname=input_freq_fname, output_cons_fname= output_cons_fname),
                      alias='mrfc_{}'.format(sample),
                      load_python=True)


# def run_pipeline():


if __name__ == "__main__":
    make_ref_from_con()



