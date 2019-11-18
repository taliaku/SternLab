# collect all summaries - bi
import re
import sys, argparse
from os import listdir
from os.path import isfile, join

COL_ALLELE = 1
COL_MEDIAN = 2
COL_MEAN = 3
COL_LOW = 4
COL_HIGH = 5
COL_DEL = 6
COL_NEU = 7
COL_ADV = 8
COL_CATEGORY = 9
COL_LEVENES_P = 10


def fits_fitness_united(input_dir, output_file):
    mypath = input_dir
    #mypath = "/Users/talzinger/Nobackup/FITS_Final_Figures/FigureAccuracyFitness_20170718_summary_bi/"
    file_list = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    summary_regex = re.compile(r'summary')
    with open(output_file, "w") as out_handler:
        out_handler.write("pos\tinferred_w\tcategory\tlevenes_p\tfilename\n")

    for filename in file_list:
        file_match = summary_regex.search( filename )
        if file_match == None:
            continue

        table_row_regex = re.compile(r'(\S)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)')

        pos_regex = re.compile(r'fitness_all_([\d]+)')
        pos_match = pos_regex.search(filename)
        pos = pos_match.group(1)

        tmpfile = open(mypath+"/"+filename, "r")
        tmpcontent = tmpfile.readlines()
        tmpfile.close()

        # print(tmpcontent[len(tmpcontent)-1])
        table_match = table_row_regex.search(tmpcontent[len(tmpcontent)-1])
        inferred_w_value = table_match.group(COL_MEDIAN)

        #table_match = table_row_regex.search(tmpcontent[len(tmpcontent)-2])
        #max_ldist = table_match.group(COL_MAXDIST)
        category = table_match.group(COL_CATEGORY)
        levenes = table_match.group(COL_LEVENES_P)

        with open(output_file, "a") as out_handler:
            out_handler.write("%s\t%s\t%s\t%s\t%s\n" % (pos, inferred_w_value, category, levenes, filename))

def main(args):

    input_dir = args.input_dir # directory contains q23_data_mutation.csv
    output_file = args.output_file

    fits_fitness_united(input_dir, output_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=str, help="the path to the fitness output directory")
    parser.add_argument("output_file", type=str, help="the path of the conjugated file")
    args = parser.parse_args(sys.argv[1:])
    main(args)