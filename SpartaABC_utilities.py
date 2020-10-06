#! /powerapps/share/python-anaconda-3.6/bin/python

import os
import glob
from file_utilities import check_dirname, check_filename
import subprocess

def make_SpartaABC_command_file(file_path, indelible_file_path, input_aln, tree_file, output_file, numberOfIterations=100000,
                                aln_mode = 0, similarity_mode=0):
    RL = subprocess.run(f"perl /sternadi/home/volume1/taliakustin/software/SpartaABC_with_INDELible_20170320_bundle/helper_scripts/get_min_max_boundaries_from_msa.pl {input_aln}".split(),
                   stdout=subprocess.PIPE)
    max_bound = RL.stdout.decode('utf-8').split("\n")[0].split(" ")[1]
    min_bound = RL.stdout.decode('utf-8').split("\n")[1].split(" ")[1]

    make_INDELible_control_file(indelible_file_path, tree_file)
    with open(file_path, "w") as handle:
        handle.write(f"_indelibleTemplateControlFile {indelible_file_path}\n")
        handle.write(f"_dawgSimulator 0\n")
        handle.write(f"_inputRealMSAFile {input_aln}\n")
        handle.write(f"_outputGoodParamsFile {output_file}\n")
        handle.write(f"_numberOfSamplesToKeep {numberOfIterations}\n")
        handle.write(f"_alignmentMode {aln_mode}\n")
        handle.write(f"_similarity_mode {similarity_mode}\n")
        handle.write(f"_minRLVal {min_bound}\n")
        handle.write(f"_maxRLVal {max_bound}\n")
        handle.write(f"_minAVal 1.001\n")
        handle.write(f"_maxAVal 2.0\n")
        handle.write(f"_minIRVal 0.0\n")
        handle.write(f"_maxIRVal 0.15\n")
        handle.write(f"_wAvgUniqueGapSize 0.0611508168659\n")
        handle.write(f"_wMSAMin 0.00428943267438\n")
        handle.write(f"_wNumGapsLenTwo 0.00183773167418\n")
        handle.write(f"_wAvgGapSize 0.158506435717\n")
        handle.write(f"_wTotNumGaps 0.000302688690478\n")
        handle.write(f"_wNumGapsLenAtLeastFour 0.000586312813355\n")
        handle.write(f"_wNumGapsLenOne 0.000943599261764\n")
        handle.write(f"_wMSAMax 0.00332419698552\n")
        handle.write(f"_wMSALen 0.0005140817814\n")
        handle.write(f"_wTotNumUniqueGaps 0.00248790123796\n")
        handle.write(f"_wNumGapsLenThree 0.00300439144851\n")


def make_INDELible_control_file(control_file_path, tree_file, type="NUCLEOTIDE", method=1, output="FASTA", fileperrep="TRUE", submodel="JC"):
    """
    for more information see: http://abacus.gene.ucl.ac.uk/software/indelible/manual/
    :param control_file_path: outpur path of indelible control file
    :param tree_file: tree file
    :param type: NUCLEOTIDE or AMINOACID or CODON (default: NUCLEOTIDE)
    :param method: 1 or 2 (default 1) (Both methods give identical results but in some situations one will be faster than the other.)
    :param output: default FASTA
    :param fileperrep:
    :param submodel:
    :return:
    """
    with open(tree_file, "r") as tree:
        tree_text = tree.read()
    with open(control_file_path, "w") as handle:
        handle.write(f"//  INDELible control file\n")
        handle.write(f"[TYPE] {type} {method}\n")
        handle.write(f"[SETTINGS]\n")
        handle.write(f"[output] {output}\n")
        handle.write(f"[fileperrep] {fileperrep}\n")
        handle.write(f"[MODEL] modelname\n")
        handle.write(f"[submodel] {submodel}\n")
        handle.write(f"[indelmodel]  POW  ? 50\n")
        handle.write(f"[indelrate]   ?\n")
        handle.write(f"[TREE] treename {tree_text}\n")
        handle.write(f"[PARTITIONS] partitionname\n")
        handle.write(f"[treename modelname  ?]\n")
        handle.write(f"[EVOLVE]     partitionname  1  ?\n")
