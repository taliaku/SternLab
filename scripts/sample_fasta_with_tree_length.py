#! /powerapps/share/python-anaconda-3.6/bin/python


from optparse import OptionParser
from file_utilities import  check_filename
from seqFileTools import unalign, remove_description
import phylogenetic_utilities
from utils.pbs_jobs import check_pbs
import pbs_runners
import os


def main():

    parser = OptionParser("usage: %prog [options]")
    parser.add_option("-i", "--input", dest="input_aln", help="input alignment file")
    parser.add_option("-n", "--number", dest="number", help="number of sequences to sample randomly")
    parser.add_option("--max_length", dest="max_length", help="total max branch length")
    parser.add_option("--min_length", dest="min_length", default = 0, help="total min branch length")


    (options, args) = parser.parse_args()
    input_aln = options.input_aln
    number = int(options.number)
    max_length = float(options.max_length)
    min_length = float(options.min_length)


    #check filenames
    if input_aln == None:
        parser.error("you must specify an alignemnt")
    input_aln = check_filename(input_aln)
    if number > 50:
        print("Warning - the number of sequences to sample is %i - this is large and may take time"  % number)


    count = 1
    sampled_file = input_aln.split(".")[0] + "_sampled_%i.aln" % number
    sampled_unaligned_file = input_aln.split(".")[0] + "_sampled_%i.fasta" % number
    #files for prank-codon
    codon_aln_file =  input_aln.split(".")[0] + "_sampled_%i.codon_aln.best.fas" % number
    codon_tree_file = input_aln.split(".")[0] + "_sampled_%i.codon_aln.best.phy_phyml_tree.txt" % number
    #files for mafft
    aln_file =  input_aln.split(".")[0] + "_sampled_%i.mafft_aln" % number
    tree_file = input_aln.split(".")[0] + "_sampled_%i.mafft_aln-NO_DESCRIPTION.phy_phyml_tree.txt" % number

    print("Iteration %i" % count)
    total_branch_length = iteration_mafft(input_aln, sampled_file, sampled_unaligned_file, aln_file, tree_file, number)

    success = True
    count += 1
    while (total_branch_length > max_length or  total_branch_length < min_length):
        print("Iteration %i" % count)
        total_branch_length = iteration_mafft(input_aln, sampled_file, sampled_unaligned_file, aln_file, tree_file, number)
        count += 1
        """
        if count > 101:
            print("Over 100 iterations!!!!")
            success = False
            break
        """
    if success:
        print("DONE")
    else:
        print("Did not succeed")


def iteration_mafft(input_aln, sampled_file, sampled_unaligned_file, aln_file, tree_file, number):
    job_id = pbs_runners.sampling_runner(input_aln, number, sampled_file, random=True)
    check_pbs(job_id)
    unalign(sampled_file, outfile=sampled_unaligned_file)
    while (os.path.isfile(sampled_unaligned_file) == False):
        continue
    job_id = pbs_runners.mafft_runner(sampled_unaligned_file, alignment=aln_file)
    while (os.path.isfile(aln_file) == False):
        continue
    no_description = remove_description(aln_file)
    #print(no_description, aln_file)
    check_pbs(job_id)
    job_id = pbs_runners.phyml_runner(no_description, phylip=False)
    check_pbs(job_id)
    while (os.path.isfile(tree_file) == False):
        continue
    total_branch_length = phylogenetic_utilities.total_branch_length(tree_file)
    print(total_branch_length)
    return(total_branch_length)

def iteration_prank(input_aln, sampled_file, sampled_unaligned_file, codon_aln_file, tree_file, number):
    job_id = pbs_runners.sampling_runner(input_aln, number, sampled_file, random=True)
    check_pbs(job_id)
    unalign(sampled_file, outfile=sampled_unaligned_file)
    job_id = pbs_runners.prank_codon_runner(sampled_unaligned_file)
    check_pbs(job_id)
    job_id = pbs_runners.phyml_runner(codon_aln_file, phylip=False)
    check_pbs(job_id)
    total_branch_length = phylogenetic_utilities.total_branch_length(tree_file)
    print(total_branch_length)
    return(total_branch_length)








if __name__ == "__main__":
    main()

