import datetime
import getpass
import unittest
import filecmp
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Project_Runner import main, create_parser
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utils.logger import pipeline_logger
import subprocess
from utils.pbs_jobs import USER_FOLDER_DICT


class TestProjectRunner(unittest.TestCase):

    def __init__(self):
        super().__init__()
        username = getpass.getuser()
        user_folder = USER_FOLDER_DICT[username]
        testing_folder = os.path.join(user_folder, 'testing')
        if not os.path.exists(testing_folder):
            os.mkdir(testing_folder)
        self.output_dir = f"{user_folder}/testing/TestProjectRunner-{datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S')}"
        log = pipeline_logger('TestProjectRunner', self.output_dir)
        log.info('Starting TestProjectRunner!')
        input_dir = '/sternadi/home/volume3/ita/pipelineTester/small_data_samples/'
        reference = '/sternadi/home/volume3/ita/pipelineTester/test_data_reference.fasta'
        parameters_dict = {
            'input_dir': '/sternadi/home/volume3/ita/pipelineTester/small_data_samples/',
            'output_dir': self.output_dir,
            'ref': '/sternadi/home/volume3/ita/pipelineTester/test_data_reference.fasta'
        }
        args = create_parser()
        """args.input_dir = parameters_dict['input_dir']
        args.output_dir = parameters_dict['output_dir']
        args.ref = parameters_dict['ref']"""
        for parameter, value in parameters_dict.items():
            setattr(args, parameter, value)

        project_runner_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                           'Project_Runner.py')
        bash_command = f"python {project_runner_path} -o {self.output_dir} -i {input_dir} -r {reference} -w Y"
        log.info(f"Running bash command: {bash_command}")
        log.info(f"This should take 3-10 minutes...")
        process = subprocess.run(bash_command.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        #output, error = process.communicate()
        #main(args)

    def assertAggregatedSummaryIsEqual(self):
        AggregatedSummaryExample = '/sternadi/home/volume3/ita/pipelineTester/small_sample_results/AggregatedSummary.csv'
        AggregatedSummaryTest = f'{self.output_dir}/AggregatedSummary.csv'
        self.assertTrue(os.path.isfile(AggregatedSummaryTest), "AggregatedSummary.csv does not exist!")
        self.assertTrue(filecmp.cmp(AggregatedSummaryTest, AggregatedSummaryExample),
                        'AggregatedSummary.csv does not match example..!')


if __name__ == '__main__':
    test = TestProjectRunner()
    test.assertAggregatedSummaryIsEqual()
