'''What I want the command line to look like:
python3 gor-train database_object -w windowsize --outfile path_to_gor_model'''

import argparse
import ml.cls_gor as gor_utils
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('dbpath', help = 'path to pickled database')
parser.add_argument('--windowsize', '-w', help = 'integer; size of the window for the GOR training', type = int, default=17)
parser.add_argument('--outfile', '-o', help = 'path to pickle the GOR model')
args = parser.parse_args()

model = gor_utils.Gor(w=args.windowsize)
dbfile = open(args.dbpath, 'rb')
datab = pickle.load(dbfile)
dbfile.close()
model.train(datab)
# print(model)
if args.outfile:
    with open(args.outfile, 'wb') as out:
        pickle.dump(model, out)
