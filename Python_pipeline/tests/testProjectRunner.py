import argparse
import unittest
import filecmp
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from Project_Runner import main, create_parser


class TestProjectRunner(unittest.TestCase):

    def __init__(self):
        super().__init__()
        print('Im init!!!!')
        self.output_dir = '/tmp/TestProjectRunner/'
        parameters_dict = {
            'input_dir': '/sternadi/home/volume3/ita/pipelineTester/small_data_samples/',
            'output_dir': self.output_dir,
            'ref': '/sternadi/home/volume3/ita/pipelineTester/test_data_reference.fasta'
        }
        args = create_parser()
        for parameter, value in parameters_dict.items():
            setattr(args, parameter, value)
        main(args)

    def assertAggregatedSummary(self):
        AggregatedSummaryExample = '/sternadi/home/volume3/ita/pipelineTester/small_sample_results/AggregatedSummary.csv'
        self.assertTrue(filecmp.cmp(f'{self.output_dir}/AggregatedSummary.csv', AggregatedSummaryExample),
                        'AggregatedSummary does not match example..!')


if __name__ == '__main__':
    test = TestProjectRunner()
    test.assertAggregatedSummary()
