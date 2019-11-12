import sys


def dfparse(dsspfile, chain):
    """
    1st argument = dssp file handler
    2nd argument = chain ID
    Traverse the dssp file and write a concatenated string with the secondary structure (16th char)"""
    mapper = {'H': 'H', 'B': 'E', 'E': 'E', 'G': 'H', 'I': 'H', 'T': '-', 'S': '-', ' ': '-'}
    flag = False
    ss = ''
    for line in dsspfile:
        if flag and line[11] == chain:
            ss = ss + line[16]
        if '#' in line:
            flag = True
    return ss


def fastaheader(pdb_id):
    return '>' + pdb_id


if __name__ == '__main__':
    '''testcode usage:
        1st argument = file with pdbid_chain list
        2nd argument = path to folder with dssp files to parse'''
    pdb_id = sys.argv[1].strip().split('_')
    dssp_file = sys.argv[2] + pdb_id[0].lower() + '.pdb.dssp'
    with open(dssp_file) as df:
        outfile = sys.argv[2] + 'blind_dssp_formatted/' + sys.argv[1] + '.dssp'
        of = open(outfile, 'w')
        of.write(fastaheader(sys.argv[1].strip()))
        of.write('\n')
        of.write(dfparse(df, pdb_id[1]))
        of.write('\n')
        of.close()

