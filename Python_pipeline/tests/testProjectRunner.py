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
    testing_path = f"/sternadi/nobackup/volume1/tests/{username}"
    if not os.path.exists(testing_path):
        os.mkdir(testing_path)
    output_dir_name = f"TestProjectRunner-{datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}"
    output_dir_path = os.path.join(testing_path, output_dir_name)
    os.mkdir(output_dir_path)
    return output_dir_path

def _omit_file(filename):
    if filename.find("tau.ac.il") == -1:
        return False
    else:
        return True

def _get_relevant_file_names(path):
    ret = {}
    for dirpath, dirnames, filenames in os.walk(path):
        relevant_file_names = [file for file in filenames if not _omit_file(file)]
        relative_path = os.path.relpath(dirpath, path)
        ret[relative_path] = relevant_file_names
    return ret


class TestProjectRunner(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.output_dir = _assign_output_dir()
        log = pipeline_logger('TestProjectRunner', self.output_dir)
        log.info('Starting TestProjectRunner!')
        self.input_dir = '/sternadi/home/volume3/ita/pipelineTester/small_data_samples/'
        self.example_output = '/sternadi/home/volume3/ita/pipelineTester/small_sample_results/'
        reference = '/sternadi/home/volume3/ita/pipelineTester/test_data_reference.fasta'
        project_runner_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                           'Project_Runner.py')
        bash_command = f"python {project_runner_path} -o {self.output_dir} -i {self.input_dir} -r {reference} -c 0"
        log.info(f"Running bash command: {bash_command}")
        log.info(f"This should take 3-10 minutes...")
        print("----------------------------------------------------------------------")
        subprocess.run(bash_command.split(), stdout=subprocess.PIPE)

    def test_files_in_dir(self):
        output_files = _get_relevant_file_names(self.output_dir)
        example_files = _get_relevant_file_names(self.example_output)
        #print(output_files)
        #print("___!!!!____!!!____!!!_____________")
        #print(example_files)
        missing_files = []
        missing_dirs = []
        for dir, files in example_files.items():
            #TODO: test this part.
            if dir not in output_files.keys():
                missing_dirs.append(dir)
                continue
            for file in files:
                if file not in output_files[dir]:
                    missing_files.append(os.path.join(dir, file))
        self.assertTrue(len(missing_dirs) == 0, f"Uh oh! These directories were not created: {missing_dirs}")
        self.assertTrue(len(missing_files) == 0, f"Whah! Looks like we are missing these files: {missing_files}")

    def test_aggregated_summary(self):
        AggregatedSummaryExample = os.path.join(self.example_output, 'AggregatedSummary.csv')
        AggregatedSummaryTest = os.path.join(self.output_dir, 'AggregatedSummary.csv')
        self.assertTrue(os.path.isfile(AggregatedSummaryTest), "AggregatedSummary.csv does not exist!")
        self.assertTrue(filecmp.cmp(AggregatedSummaryTest, AggregatedSummaryExample),
                        'AggregatedSummary.csv does not match example..!')


if __name__ == '__main__':
    unittest.main()
