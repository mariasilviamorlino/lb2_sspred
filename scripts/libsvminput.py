'''Need to specify: path to pickled database, path to file where I want to save the input
I want the command line to look like:
python3 libsvminput.py -i <path> -o <path>'''

import ml.dataset_preprocess as dsp
import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('dataset', help = 'path to pickled dataset')
parser.add_argument('--outfile', help = 'path to file where the correctly formatted libsvm input will be saved')
args = parser.parse_args()

with open(args.dataset, 'rb') as dst:
    datab = pickle.load(dst)
if args.outfile:
    datab.makelibsvminput(filepath=args.outfile)
