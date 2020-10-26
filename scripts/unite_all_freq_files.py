import sys
sys.path.append("..")
from freqs_utilities import unite_all_freq_files
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--d", dest="directory", help="directory with freqs files")
    parser.add_option('-o', '--o', dest='output csv file')

    (options, args) = parser.parse_args()
    directory = options.directory
    output = options.output

    unite_all_freq_files(directoyr, output)
