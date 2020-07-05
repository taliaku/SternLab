import datetime
import getpass
import unittest
import filecmp
import os
import sys
import glob
import subprocess
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utils.logger import pipeline_logger
from utils.pbs_jobs import USER_FOLDER_DICT


def _assign_output_dir():
    username = getpass.getuser()
    user_folder = USER_FOLDER_DICT[username]
    testing_folder = os.path.join(user_folder, 'testing')
    if not os.path.exists(testing_folder):
        os.mkdir(testing_folder)
    return f"{user_folder}/testing/TestProjectRunner-{datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}"


def _omit_file(filename):
    if filename.find("tau.ac.il") == -1:
        return False
    else:
        return True


class TestProjectRunner(unittest.TestCase):

    def setUp(self):
        self.output_dir = _assign_output_dir()
        log = pipeline_logger('TestProjectRunner', self.output_dir)
        log.info('Starting TestProjectRunner!')
        self.input_dir = '/sternadi/home/volume3/ita/pipelineTester/small_data_samples/'
        reference = '/sternadi/home/volume3/ita/pipelineTester/test_data_reference.fasta'
        project_runner_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                           'Project_Runner.py')
        bash_command = f"python {project_runner_path} -o {self.output_dir} -i {self.input_dir} -r {reference} -c 0"
        log.info(f"Running bash command: {bash_command}")
        log.info(f"This should take 3-10 minutes...")
        subprocess.run(bash_command.split(), stdout=subprocess.PIPE)

    def test_files_in_dir(self):
        output_files = [file for file in glob.glob(self.output_dir) if not _omit_file(file)]
        example_files = [file for file in glob.glob(self.input_dir) if not _omit_file(file)]
        missing_files = [file for file in example_files if file not in output_files]
        self.assertTrue(len(missing_files) != 0, f"Whoops! we are missing these files: {missing_files}")

    def test_aggregated_summary(self):
        AggregatedSummaryExample = '/sternadi/home/volume3/ita/pipelineTester/small_sample_results/AggregatedSummary.csv'
        AggregatedSummaryTest = f'{self.output_dir}/AggregatedSummary.csv'
        self.assertTrue(os.path.isfile(AggregatedSummaryTest), "AggregatedSummary.csv does not exist!")
        self.assertTrue(filecmp.cmp(AggregatedSummaryTest, AggregatedSummaryExample),
                        'AggregatedSummary.csv does not match example..!')


if __name__ == '__main__':
    unittest.main()
