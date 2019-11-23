'''Launch upon completing a prediction with libsvm to incorporate the results in the dataset. Command:
python3 add-svm-pred.py <dbpath> <svmfile>'''

import pickle
import ml.dataset_preprocess as dbm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('dbpath', help = 'path to pickled database to endow the predictions')
parser.add_argument('svmpath', help = 'path to svm file containing the predictions')
args = parser.parse_args()

with open(args.dbpath, 'rb') as dbhandle:
    datab = pickle.load(dbhandle)

datab.svmpred(args.svmpath)

with open(args.dbpath, 'wb') as dbhandle:
    pickle.dump(datab, dbhandle)
