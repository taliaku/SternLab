import argparse
import os
import subprocess
import sys
import pandas as pd
from matplotlib import pyplot as plt
from Join import wrangle_freqs_df


STERNLAB_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(STERNLAB_PATH)
from utils.logger import pipeline_logger
from utils.pbs_jobs import create_pbs_cmd
from utils.runner_utils import submit_wait_and_log

#TODO: make this run!

def _create_python_output_folder(output_folder):
    python_output_path = os.path.join(output_folder, 'python_output')
    if not os.path.exists(python_output_path):
        os.mkdir(python_output_path)
    return python_output_path


def _create_perl_output_folder(output_folder):
    perl_output_path = os.path.join(output_folder, 'perl_output')
    if not os.path.exists(perl_output_path):
        os.mkdir(perl_output_path)
        perl_tmp_folder = os.path.join(perl_output_path, 'tmp')
        if not os.path.exists(perl_tmp_folder):
            os.mkdir(perl_tmp_folder) # seems like this is required by the pipeline..
    return perl_output_path


def _get_python_runner_flags(input_data_folder, output_folder):
    ret = {}
    input_dir_name = os.path.basename(os.path.normpath(input_data_folder))
    ret['i'] = os.path.join(input_data_folder, input_dir_name[:2])
    ret['o'] = _create_python_output_folder(output_folder)
    ret['runner_path'] = os.path.join(STERNLAB_PATH, 'Python_pipeline', 'Runner.py')
    return ret


def create_runners_cmdfile(input_data_folder, output_folder, reference_file, alias):
    perl_output_path = _create_perl_output_folder(output_folder)
    perl_runner_path = os.path.join(STERNLAB_PATH, 'pipeline_runner.py')
    python_runner_flags = _get_python_runner_flags(input_data_folder=input_data_folder, output_folder=output_folder)
    perl_runner_cmd = f"python {perl_runner_path} -i {python_runner_flags['o']} -o {perl_output_path} -r {reference_file} " \
                      f"-NGS_or_Cirseq 1"
    """ 
    the input for perl_runner_cmd is the output of python_runner_cmd because the python pipeline first created
    the fastq files which both pipelines use.
    """
    python_runner_cmd = f"python {python_runner_flags['runner_path']} -i {python_runner_flags['i']} " \
                        f"-o {python_runner_flags['o']} -r {reference_file} -m RS -L {output_folder} -x 1 -s 1"
    cmds = perl_runner_cmd + "\n" + python_runner_cmd
    cmd_file_path = os.path.join(output_folder, 'compare_pipelines.cmd')
    create_pbs_cmd(cmdfile=cmd_file_path, alias=alias, cmds=cmds)
    return cmd_file_path


def get_single_freq_file_path(path, freq_file_suffix):
    freq_files = [f for f in os.listdir(path) if f.find(freq_file_suffix) != -1]
    if len(freq_files) == 0:
        raise Exception(f"Could not find file containing {freq_file_suffix} file in {path} !")
    elif len(freq_files) > 1:
        raise Exception(f"Found more than one file containing '{freq_file_suffix}' in {path}..!")
    return os.path.join(path, freq_files[0])


def get_python_freqs(python_output_path):
    freq_file_path = get_single_freq_file_path(python_output_path, "merge.freqs.csv")
    return pd.read_csv(freq_file_path).set_index(['ref_position', 'base'], drop=True)


def get_perl_freqs(perl_output_path):
    """
    Wrangles perl freqs file to fit with Python freqs file.
    """
    perl_freqs_path = get_single_freq_file_path(perl_output_path, '.freqs')
    pe_df = pd.read_csv(perl_freqs_path, sep='\t', index_col=[0, 1, 2], usecols=[0, 1, 2, 3, 4])
    pe_df.reset_index(inplace=True)
    pe_df.columns = ['ref_position', 'base', 'freq', 'ref_base', 'coverage']
    pe_df.set_index(['ref_position', 'ref_base', 'base'], drop=True, inplace=True)
    pe_df['base_counter'] = pe_df['coverage'] * pe_df['freq']
    pe_df.drop(['freq'], axis=1, inplace=True)
    pe_df = wrangle_freqs_df(pe_df)
    return pe_df


def get_freqs_data(output_folder):
    python_output_path = _create_python_output_folder(output_folder)
    py_df = get_python_freqs(python_output_path)
    perl_output_path = _create_perl_output_folder(output_folder)
    pe_df = get_perl_freqs(perl_output_path)
    data = {'py': py_df, 'pe': pe_df}
    return data


def plot_coverage(py_df, pe_df, output_folder):
    plt.figure(figsize=(20,10))
    plt.bar(pe_df.index.get_level_values('ref_position'), pe_df['coverage'], label='perl pipeline')
    plt.bar(py_df.index.get_level_values('ref_position'), py_df['coverage'], label='python pipeline')
    plt.xlabel('ref_position')
    plt.ylabel('coverage')
    plt.title('Pipeline Coverage: perl vs python')
    plt.legend()
    plt.savefig(os.path.join(output_folder, 'coverage.png'))


def _get_mismatching_bases(zero_rank_data):
    joined = zero_rank_data['pe'].join(zero_rank_data['py'], rsuffix='_py', lsuffix='_pe')
    joined.dropna(subset=['ref_base_py', 'ref_base_pe'], inplace=True)
    mismatching_bases = joined[joined.base_py != joined.base_pe]
    return mismatching_bases


def create_analyze_data_cmdfile(output_folder, alias):
    cmd_file_path = os.path.join(output_folder, 'analyze_data.cmd')
    this_module = os.path.basename(os.path.normpath(os.path.abspath(__file__)))[:-3]
    output_folder_string = '"' + output_folder + '"'
    cmd = f"cd {os.path.join(STERNLAB_PATH, 'Python_pipeline')}; python -c " \
          f"'from {this_module} import analyze_data; analyze_data({output_folder_string})'"
    create_pbs_cmd(cmdfile=cmd_file_path, alias=alias, cmds=cmd)
    return cmd_file_path


def analyze_data(output_folder):
    """
    This function is called by PBS via the cmdfile created by create_analyze_data_cmdfile
    """
    data = get_freqs_data(output_folder)
    zero_rank_data = {key: df[df['rank'] == 0] for key, df in data.items()}
    plot_coverage(zero_rank_data['py'], zero_rank_data['pe'], output_folder)
    joined = data['pe'].join(data['py'], rsuffix='_py', lsuffix='_pe', how='outer')
    joined.to_csv(os.path.join(output_folder, 'data.csv'))


def merge_fastq_files(input_data_folder, output_folder, reference_file):
    """ Merge fastq files using the python_runner """
    python_runner_flags = _get_python_runner_flags(input_data_folder=input_data_folder, output_folder=output_folder)
    merge_fastq_cmd = f"python {python_runner_flags['runner_path']} -i {python_runner_flags['i']} " \
                      f"-o {python_runner_flags['o']} -r {reference_file} -m RS -L {output_folder} -x 1 -s 0 -e 0"
    subprocess.run(merge_fastq_cmd.split(), stdout=subprocess.PIPE)


def main(args):
    input_data_folder = args.input_data_folder
    output_folder = args.output_folder
    reference_file = args.reference_file
    just_analyze = args.just_analyze
    alias = 'ComparePipelines'
    log = pipeline_logger(alias, output_folder)
    if not just_analyze:
        log.info(f"Comparing pipelines on data from {input_data_folder} and outputting to {output_folder}")
        merge_fastq_files(input_data_folder=input_data_folder, output_folder=output_folder, reference_file=reference_file)
        compare_cmd_path = create_runners_cmdfile(input_data_folder, output_folder, reference_file, alias)
        submit_wait_and_log(compare_cmd_path, log, alias)
    log.info(f"Analyzing data...")
    alias = 'ComparePipelines-AnalyzeData'
    analyze_cmd_path = create_analyze_data_cmdfile(output_folder, alias)
    submit_wait_and_log(analyze_cmd_path, log, alias)
    log.info("Done.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_data_folder",
                        help="Fasta files to run both pipelines on",
                        required=True)
    parser.add_argument("-o", "--output_folder",
                        help="Where you want the output files to go",
                        required=True)
    parser.add_argument("-r", "--reference_file",
                        required=True)
    parser.add_argument("-j", "--just_analyze", default=False,
                        help='True will skip the pipelines and just analyze the output. default is False.')
    args = parser.parse_args()
    main(args)
