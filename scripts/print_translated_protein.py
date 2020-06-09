#! /usr/local/python_anaconda/bin/python3.4

from optparse import OptionParser
from file_utilities import check_filename, check_dirname
from seqFileTools import translate_file_print_seq

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-f", "--file", dest="file", help="file to translate")

    (options, args) = parser.parse_args()


    file = options.file

    translate_file_print_seq(file)

if __name__ == "__main__":
    main()

