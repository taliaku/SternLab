#! /powerapps/share/python-anaconda-3.6/bin/python

import sys
sys.path.insert(0,'/sternadi/home/volume1/taliakustin/SternLab')
from optparse import OptionParser
import os
import glob
from file_utilities import check_filename, check_dirname
from dirSel_utilities import write_params_file
import pbs_runners
from netMHCpan_utilities import netMHCpan_to_csv
from pandas_utilities import merge_dfs
import pandas as pd
def main():
    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-a", "--allele", dest="allele")


    (options, args) = parser.parse_args()
    allele = options.allele
    inputs = glob.glob("/sternadi/home/volume3/taliakustin/phyVirus_analysis/codon_aln_back_to_fasta/*translated_*")
    print(len(inputs))

    if glob.glob(f"/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/{allele}_results.csv") != []:
        print(f"allele {allele} exists")
        return

    output_dir = "/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/%s/" % allele
    out_df = pd.DataFrame()
    if not os.path.isdir(output_dir):
        print(os.system("mkdir %s" % output_dir))
        print(output_dir)

    for i in range(len(inputs)):
        number = inputs[i].split("translated_")[1].split(".fasta")[0]
        output =   output_dir + inputs[i].split("-translated")[0].split("/")[-1] + ".%s.%s.netMHCpan_results" % (allele, number)
        #output = inputs[i].split("-translated")[0] + ".%s.netMHCpan_results" % number
        if os.path.isfile(output):
            print(output)
            res = netMHCpan_to_csv(output, to_return=True)
            if res.empty:
                print("no allele!!")
                continue
            out_df = out_df.append(res, ignore_index=True)
            continue
        os.system("/sternadi/home/volume1/taliakustin/software/netMHCpan-4.0/netMHCpan %s -t 0.2 -l 9 -a '%s'> %s"
                  % (inputs[i], allele, output))
        print(output)
        res = netMHCpan_to_csv(output, to_return=True)
        if res.empty:
            print("no allele!!")
            continue
        out_df = out_df.append(res, ignore_index=True)
    out_df.to_csv("/sternadi/home/volume3/taliakustin/phyVirus_analysis/netMHCpan4/{}_results.csv".format(allele))
    print("rm -r {}".format(output_dir))
    os.system("rm -r {}".format(output_dir))






if __name__ == "__main__":
    main()

