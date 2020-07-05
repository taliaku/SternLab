import time
import unittest
import filecmp
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Project_Runner import main, create_parser
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utils.logger import pipeline_logger
import subprocess



class TestProjectRunner(unittest.TestCase):

    def __init__(self):
        super().__init__()
        self.output_dir = '/tmp/TestProjectRunner' #TODO: should it really be this folder..?
        log = pipeline_logger('TestProjectRunner', self.output_dir)
        log.info('Starting TestProjectRunner!')
        log.info(f"{time.time()}")
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
        subprocess.run(bash_command.split(), stdout=subprocess.PIPE)
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
