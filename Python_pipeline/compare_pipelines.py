import argparse
import os
import sys
import pandas as pd
from matplotlib import pyplot as plt
from Join import wrangle_freqs_df
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from utils.logger import pipeline_logger
from utils.pbs_jobs import create_pbs_cmd
from utils.runner_utils import submit_wait_and_log

#TODO: make this run!

def _get_python_output_path(output_folder):
    return os.path.join(output_folder, 'python_output')


def _get_perl_output_path(output_folder):
    return os.path.join(output_folder, 'perl_output')


def create_runners_cmdfile(input_data_folder, output_folder, reference_file, alias):
    perl_output_path = _get_perl_output_path(output_folder)
    python_output_path = _get_python_output_path(output_folder)
    perl_runner_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'pipeline_runner.py')
    perl_runner_cmd = f"python {perl_runner_path} -i {input_data_folder} -o {perl_output_path} -r {reference_file} -NGS_or_Cirseq 1"
    input_dir_name = os.path.basename(os.path.normpath(input_data_folder))  # python pipeline needs this..
    # TODO: call python pipeline in a way that makes sense
    python_runner_cmd = f"python Python_pipeline/Runner.py -i {os.path.join(input_data_folder, input_dir_name[:2])} " \
                        f"-o {python_output_path} -r {reference_file} -m RS"
    cmds = perl_runner_cmd + python_runner_cmd
    cmd_file_path = os.path.join(output_folder, 'compare_pipelines.cmd')
    create_pbs_cmd(cmdfile=cmd_file_path, alias=alias, cmds=cmds)
    return cmd_file_path


def get_single_freq_file_path(path, freq_file_suffix):
    freq_files = [f for f in os.listdir(path) if f.find(freq_file_suffix)]
    if len(freq_files) == 0:
        raise (ValueError, f"Could not find file containing {freq_file_suffix} file in {path} !")
    elif len(freq_files) > 1:
        raise (ValueError, f"Found more than one file containing '{freq_file_suffix}' in {path}..!")
    return os.path.join(path, freq_files[0])


def get_python_freqs(python_output_path):
    freq_file_path = get_single_freq_file_path(python_output_path, "merge.freqs.csv")
    return pd.read_csv(freq_file_path).set_index('ref_position', drop=True)


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
    python_output_path = _get_python_output_path(output_folder)
    py_df = get_python_freqs(python_output_path)
    perl_output_path = _get_perl_output_path(output_folder)
    pe_df = get_perl_freqs(perl_output_path)
    data = {'py': py_df, 'pe': pe_df}
    return data


def plot_coverage(py_df, pe_df, output_folder):
    plt.figure(figsize=(20,10))
    plt.bar(pe_df.index, pe_df['coverage'], label='perl pipeline')
    plt.bar(py_df.index, py_df['coverage'], label='python pipeline')
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
    cmd = f"python -c 'from {this_module} import analyze_data; analyze_data({output_folder_string})'"
    create_pbs_cmd(cmdfile=cmd_file_path, alias=alias, cmds=cmd)
    return cmd_file_path


def analyze_data(output_folder):
    """
    This function is called by PBS via the cmdfile created by create_analyze_data_cmdfile
    """
    data = get_freqs_data(output_folder)
    zero_rank_data = {key: df[df['rank'] == 0] for key, df in data.items()}
    plot_coverage(zero_rank_data['py'], zero_rank_data['pe'], output_folder)
    mismatching_bases = _get_mismatching_bases(zero_rank_data)
    mismatching_bases.to_csv('mismatching_bases.csv')
    print(f"There are {len(mismatching_bases)} mismatching bases!")


def main(args):
    input_data_folder = args.input_data_folder
    output_folder = args.output_folder
    reference_file = args.reference_file
    alias = 'ComparePipelines'
    log = pipeline_logger(alias, output_folder)
    log.info(f"Comparing pipelines on data from {input_data_folder} and outputting to {output_folder}")
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
    args = parser.parse_args()
    main(args)
