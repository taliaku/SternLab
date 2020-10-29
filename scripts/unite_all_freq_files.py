import os,sys,inspect
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
sys.path.insert(0,parentdir)
from freqs_utilities import unite_all_freq_files
from optparse import OptionParser

if __name__ == '__main__':
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-d", "--d", dest="directory", help="directory with freqs files")
    parser.add_option('-o', '--o', dest='output', help='output csv file')

    (options, args) = parser.parse_args()
    directory = options.directory
    output = options.output

    unite_all_freq_files(directory, output)
