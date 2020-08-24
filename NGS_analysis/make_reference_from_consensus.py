#! /usr/local/python_anaconda/bin/python3.4
from Bio import SeqIO,Seq
from Bio.Seq import MutableSeq
from optparse import OptionParser
import numpy as np
import pandas as pd

def replace_base(row, mut_reference):
    mut_reference[int(row["Pos"])-1]=row["Base"]

def main():
    parser = OptionParser("usage: %prog [options]\nTry running %prog --help for more information")
    parser.add_option("-f", "--fasta", dest="fasta", help="reference fasta file")
    parser.add_option("-p", "--freqs", dest="freqs", help="frequency file")
    parser.add_option("-o", "--output_file", dest="output_file", help="output fasta file")
    parser.add_option("-c", "--min_coverage", dest="coverage", help="minimal coverage of position to be updated (default: 1000)", type="int")
    # TODO- allow configurable output sequence name (currently: Reference name)

    (options, args) = parser.parse_args()
    fasta = options.fasta
    freqs = options.freqs
    output_file = options.output_file
    min_coverage = 100
    if options.coverage:
        min_coverage= options.coverage

    if options.fasta is None or options.fasta is None or options.output_file is None:
        parser.error("Missing file input")

    make_reference_from_consensus(fasta, freqs, min_coverage, output_file)


def make_reference_from_consensus(fasta, freqs, min_coverage, output_file):
    reference = None
    input_record = None
    mutable_reference = None
    with open(fasta, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            reference = record.seq
            input_record = record
    mutable_reference = MutableSeq(str(reference))
    data_freqs = pd.read_table(freqs)
    data_freqs = data_freqs[data_freqs["Read_count"] >= min_coverage]
    data_freqs = data_freqs[data_freqs["Rank"] == 0]
    # data_freqs = data_freqs[data_freqs["Base"] != data_freqs["Ref"]]
    data_freqs = data_freqs[data_freqs["Pos"] == np.round(data_freqs['Pos'])]  # remove insertion
    data_freqs = data_freqs[data_freqs["Base"] != "-"]  # remove deletion
    data_freqs.apply(replace_base, args=(mutable_reference,),
                     axis=1)  # updates mutable reference to hold correct consensus
    new_sequence = Seq.Seq(str(mutable_reference))
    input_record.seq = new_sequence

    # TODO- change default output sequence name, from ref name to ref+sample
    # sample_name = os.path.splitext(os.path.basename(freqs))[0]
    # input_record.description = input_record.description + sample_name

    with open(output_file, "w") as output_handle:
        SeqIO.write(input_record, output_handle, "fasta")


if __name__ == "__main__":
    main()