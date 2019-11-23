'''Opens pandas series containing performance results and builds a dataframe including the stacked pandas series
+ a row with mean and standard error for each measure
Usage:
python3 performance_aggregatecv.py <idfile> <path_to_folder> --outfile <report_path> --gor/--svm'''

import pickle
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('idfile', help = 'path to list of dataset results to aggregate')
parser.add_argument('result_folder', help = 'path to folder containing the pickled results')
parser.add_argument('--outfile', '-o', help = 'path to the output file if you wish to save the dataframe', default = False)
parser.add_mutually_exclusive_group()
parser.add_argument('--gor', help = 'specify if the results are referred to gor', action='store_true')
parser.add_argument('--svm', help = 'specify if the results are referred to svm', action='store_true')
args = parser.parse_args()

with open(args.idfile) as ids:
    ids_list = ids.readlines()

list_of_series, list_of_labels = list(), list()
# means, stderrs = [0 for i in range(14)], [0 for i in range(14)]

idx = 1
if args.gor: extension = 'gor'
if args.svm: extension = 'svm'

for i in ids_list:
    i = i.strip()
    if i != '':
        respath = args.result_folder + '/' + i + extension + '.pickle'
        with open(respath, 'rb') as reshandle:
            s = pickle.load(reshandle)
        list_of_series.append(s)
        list_of_labels.append('cv '+str(idx))
        idx += 1

# list_of_labels.append('mean')
# list_of_labels.append('stderr')
# list_of_series.append(means)
# list_of_series.append(stderrs)

aggregate = pd.DataFrame(list_of_series, index=list_of_labels)
means = aggregate.mean()
stderr = aggregate.sem()
aggregate = aggregate.append(means, ignore_index = True)
aggregate = aggregate.append(stderr, ignore_index = True)

if args.outfile:
    with open(args.outfile, 'wb') as outhandle:
        pickle.dump(aggregate, outhandle)
else:
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(aggregate)
