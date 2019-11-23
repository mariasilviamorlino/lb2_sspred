'''what the command line should look like:
python3 performance_assessment.py <database object> --outfile <outfile_path> --gor/--svm
assesses the performance of one data set for one of the two methods (gor or svm)
'''

import ml.dataset_preprocess as dbm
import ml.performance as perf
import ml.sov as sov
import argparse
import pickle
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('dbpath', help = 'path to pickled database containing the predictions')
parser.add_argument('--outfile', '-o', help = 'path to the output file if you wish to save the results', default = False)
parser.add_mutually_exclusive_group()
parser.add_argument('--gor', help = 'specify if you want to test the performance for gor', action='store_true')
parser.add_argument('--svm', help = 'specify if you want to test the performance for svm', action='store_true')
args = parser.parse_args()

with open(args.dbpath, 'rb') as dbhandle:
    datab = pickle.load(dbhandle)

sov_measure = {'H': 0, 'E': 0, '-': 0, 'tot': 0}
sov_ctr = {'H': 0, 'E': 0, '-': 0, 'tot': 0}
gor_confmat = perf.confusion_matrix()

if args.gor: pred = 3
if args.svm: pred = 4

for entry in datab:
    gor_confmat.populate(entry[pred], entry[1])
    for ss in 'HE-':
        if ss in entry[1]:
            sov_measure[ss] += sov.sov(entry[pred], entry[1], ss)
            sov_ctr[ss] += 1
        # else:
        #     sov_measure[ss] = float('nan')
    sov_measure['tot'] += sov.multiclass_sov(entry[pred], entry[1])
    sov_ctr['tot'] += 1
for i in sov_measure.keys():
    sov_measure[i] = sov_measure[i]/sov_ctr[i]

measures = {'tpr_H': gor_confmat.sen('H'), 'tpr_E': gor_confmat.sen('E'), 'tpr_C': gor_confmat.sen('-'),
            'ppv_H': gor_confmat.ppv('H'), 'ppv_E': gor_confmat.ppv('E'), 'ppv_C': gor_confmat.ppv('-'),
            'mcc_H': gor_confmat.mcc('H'), 'mcc_E': gor_confmat.mcc('E'), 'mcc_C': gor_confmat.mcc('-'), 'q3': gor_confmat.q3(),
            'sov_H': sov_measure['H'], 'sov_E': sov_measure['E'], 'sov_C': sov_measure['-'], 'sov_tot': sov_measure['tot']}
measures = pd.Series(measures)


if args.outfile:
    with open(args.outfile, 'wb') as outhandle:
        pickle.dump(measures, outhandle)
else:
    print("Measures of performance:\n{0}".format(measures))

# TODO:
# produce a script to compute aggregated cv performances
