'''What I want the command line to look like:
python3 gor-predict database_object path_to_gor_model'''

import argparse
import ml.cls_gor as gor_utils
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('dbpath', help = 'path to pickled database to update with predictions')
parser.add_argument('gorpath', help = 'path to the pickled GOR model')
args = parser.parse_args()

with open(args.gorpath, 'rb') as gorhandle:
    model = pickle.load(gorhandle)
with open(args.dbpath, 'rb') as dbhandle:
    datab = pickle.load(dbhandle)

updated_datab = model.predict(datab)
dbhandle = open(args.dbpath, 'wb')
pickle.dump(updated_datab, dbhandle)
dbhandle.close()
