#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file
import pbs_runners

def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--index", dest="index", help="pbs index")
    parser.add_option("-n", "--num_run", dest="num_run")
    parser.add_option("-a", "--allele", dest="allele")


    (options, args) = parser.parse_args()
    index = options.index
    num_run =int(options.num_run)
    allele = options.allele
    index_dics = {}
    all_inputs = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/*translated_*")
    print(len(all_inputs))
    if num_run == 1:
        inputs = all_inputs[:10000]
    elif num_run == 2:
        inputs = all_inputs[10001:20000]
    elif num_run == 3:
        inputs = all_inputs[20001:30000]
    elif num_run == 4:
        inputs = all_inputs[30001:40000]
    elif num_run == 5:
        inputs = all_inputs[40001:50000]
    elif num_run == 6:
        inputs = all_inputs[50001:60000]
    print(len(inputs))
    output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/%s/" % allele
    if not os.path.isdir(output_dir):
        os.system("mkdir %s" % output_dir)
    for i in range(len(inputs)):
        number = inputs[i].split("translated_")[1].split(".fasta")[0]
        output =   output_dir + inputs[i].split("-translated")[0].split("/")[-1] + ".%s.%s.netMHCpan_results" % (allele, number)
        #output = inputs[i].split("-translated")[0] + ".%s.netMHCpan_results" % number
        index_dics[str(i)] = (inputs[i], output)
    print("/sternadi/home/volume1/taliakustin/software/netMHCpan-4.0/netMHCpan %s -t 0.2 -l 9 -a '%s'> %s"
              % (index_dics[index][0], allele, index_dics[index][1] ))
    if os.path.isfile(index_dics[index][1]):
        print("exists")
        return
    os.system("/sternadi/home/volume1/taliakustin/software/netMHCpan-4.0/netMHCpan %s -t 0.2 -l 9 -a '%s'> %s"
              % (index_dics[index][0], allele, index_dics[index][1] ))

if __name__ == "__main__":
    main()

