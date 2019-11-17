# collect all summaries - bi
import re
import sys, argparse
from os import listdir
from os.path import isfile, join

COL_FROM = 1
COL_TO = 2
COL_MEDIAN = 3
COL_MAD = 4
COL_MIN = 5
COL_MAX = 6
COL_LEVENES_P = 7




def fits_mutation_united(input_dir, output_file):
    mypath = input_dir
    #mypath = "/Users/talzinger/Nobackup/FITS_Final_Figures/FigureAccuracyFitness_20170718_summary_bi/"
    file_list = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    summary_regex = re.compile(r'summary')
    with open(output_file, "w") as out_handler:
        out_handler.write("pos\tinferred_mu\trev_inferred_mu\tlevenes_p\tfilename\n")

    for filename in file_list:
        file_match = summary_regex.search(filename)
        if file_match == None:
            continue

        table_row_regex = re.compile(r'(\S)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)')

        pos_regex = re.compile(r'mutation_syn_([\d]+)')
        pos_match = pos_regex.search(filename)
        pos = pos_match.group(1)

        tmpfile = open(mypath+"/"+filename, "r")
        tmpcontent = tmpfile.readlines()
        tmpfile.close()

        # print(tmpcontent[len(tmpcontent)-2])
        table_match = table_row_regex.search(tmpcontent[len(tmpcontent)-2])
        inferred_w_value = table_match.group(COL_MEDIAN)
        rev_table_match = table_row_regex.search(tmpcontent[len(tmpcontent)-1])
        rev_inferred_w_value = rev_table_match.group(COL_MEDIAN)

        #table_match = table_row_regex.search(tmpcontent[len(tmpcontent)-2])
        #max_ldist = table_match.group(COL_MAXDIST)
        # category = table_match.group(COL_CATEGORY)
        levenes = table_match.group(COL_LEVENES_P)

        with open(output_file, "a") as out_handler:
            out_handler.write("%s\t%s\t%s\t%s\t%s\n" % (pos, inferred_w_value, rev_inferred_w_value, levenes, filename))

def main(args):

    input_dir = args.input_dir # directory contains q23_data_mutation.csv
    output_file = args.output_file

    fits_mutation_united(input_dir, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=str, help="the path to the mutation output directory")
    parser.add_argument("output_file", type=str, help="the path of the conjugated file")
    args = parser.parse_args(sys.argv[1:])
    main(args)