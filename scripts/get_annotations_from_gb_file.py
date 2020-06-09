#! /powerapps/share/python-anaconda-3.6/bin/python

from Bio import SeqIO, Seq, SeqRecord
import tqdm
import argparse


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', metavar='i', type=str, dest="input_file",
						help='input gene bank file')
	parser.add_argument('-f', dest='features',  nargs='+',
						help='features')

	args = parser.parse_args()
	file = args.input_file
	feature_list = args.features
	print(feature_list)
	gb = list(SeqIO.parse(file,format='genbank'))
	if feature_list == None:
		all_features = True
		features = {}
	else:
		features = {x:[] for x in feature_list}
	for i in tqdm.tqdm(gb):
		for feature in i.features:
			if feature.type in features.keys() or all_features:
				seq = SeqRecord.SeqRecord(feature.location.extract(i.seq))
				seq.id = i.id
				seq.name = i.name
				seq.description = i.description
				if feature.type not in features.keys():
					features[feature.type] = []
				features[feature.type].append(seq)
	for f in features:
		if features[f] == []:
			continue
		output = file.split(".")[0] + "_{}.fasta".format(f.replace("'", ""))
		SeqIO.write(features[f], output, "fasta")


if __name__ == "__main__":
    main()

