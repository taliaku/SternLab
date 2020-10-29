#! /usr/local/python_anaconda/bin/python3.4

import os,sys,inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0,parentdir)
from optparse import OptionParser
from file_utilities import check_filename
from seqFileTools import convert_fasta_to_phylip

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--input", dest="input", help="fasta sequence file")
    parser.add_option("-o", "--output", dest="output", help="phylip output file", default=None)

    (options, args) = parser.parse_args()
    input = options.input
    output = options.output
    input = check_filename(input)
    if output != None:
        ouput = check_filename(output, Truefile=False)

   
    output = convert_fasta_to_phylip(input, output)
    print("converted %s fasta file into %s" % (input, output))

if __name__ == "__main__":
    main()
