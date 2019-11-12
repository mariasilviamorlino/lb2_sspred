"""usage:
    1st argument = csv file
    2nd argument = blastclust file
    3rd argument = output file"""

from statstools import dictpprint
import sys
from math import inf


def minimize(contestants, ft_mapping, obj_ft):
    """
    1st argument = iterable of items
    2nd argument = mapping (dict) of each item to a dict of features
    3rd argument = objective feature to minimize
    traverse the iterables in search of the one with the feature with minimum value"""
    best_ft = inf
    best = None
    for i in contestants:
        if float(ft_mapping[i][obj_ft]) < best_ft:
            best_ft = float(ft_mapping[i][obj_ft])
            best = i
    return (best, best_ft)


# parsing of .csv file (must contain fields: PDB ID,Chain ID,Resolution,Sequence,Chain Length)

features = dict()
with open(sys.argv[1]) as blind_set:
    # dorky trick to remove the headers which we don't care about
    blind_set.readline()
    for line in blind_set:
        ft_list = line.strip().split('"')
        if len(ft_list) > 1:  # another dorky trick to remove the last line which is just \n
            key = ft_list[1] + '_' + ft_list[3]
            features[key] = {
                    'resol': ft_list[5],
                    'fasta': ft_list[7],
                    'chainlen': ft_list[9]
                    }

# once generated, the features variable can also be imported elsewhere

# mapping resolution and choosing the best (and writing the list in a file)

with open(sys.argv[2]) as clusters:
    representatives = list()
    for line in clusters:
        pdb_ids = line.strip().split()
        bst = minimize(pdb_ids, features, 'resol')
        representatives.append(bst)
#print(representatives)
#print(len(representatives))  # OK, a representative is retrieved from each cluster

with open(sys.argv[3], 'w+') as outfile:
    for i in representatives:
        outfile.write(i[0])
        outfile.write('\n')

