'''Need to specify: path to file with ids, path to folder with profiles and dssps, where to pickle the dataset
I want the command line to look like:
python3 create-dataset.py -i <path> -p <path> -s <path> -o <path>'''

import ml.dataset_preprocess as dsp
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('ids', help = 'path to file containing the ids to process')
parser.add_argument('pssm', help = 'path to folder containing the pssm files')
parser.add_argument('sec_struc', help = 'path to folder containing the dssp strings')
parser.add_argument('--outfile', help = 'path to file where the dataset object will be saved')
parser.add_argument('--test', help = 'include if the dataset is built for the purpose of testing; the value of the optional argument is the path to the fasta folder for one hot encoding')
args = parser.parse_args()

idss = open(args.ids)
datab = dsp.Database(idss, args.pssm, args.sec_struc)
print(len(datab.db))
if args.test:
    datab.onehot(fastapath=args.test)
if args.outfile:
    with open(args.outfile, 'wb') as outf:
        pickle.dump(datab, outf)
