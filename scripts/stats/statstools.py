import sys


amino_acids = 'ACDEFGHIKLMNPQRSTVWXY'


def compo(seq, counters):
    """1st argument: fasta or dssp file
    2nd argument: dictionary of counters
    return: updated dictionary
    Traverse the string and count the characters
    In our framework, can be used both for residue composition and ss abundancy,
    as well as SCOP and taxonomy stats"""
    for i in seq:
        counters[i] = counters.get(i, 0) + 1
    return counters


def ssperresidue(fasta, dssp, stats):
    """1st argument: fasta file
    2nd argument: dssp file
    3rd argument: dictionary
    return: updated dictionary
    Traverse both the fasta and the dssp and increment the corresponding values"""
    chars = list(zip(dssp, fasta))
    for i in chars:
        stats[i[0]] = stats.get(i[0], dict())
        stats[i[0]][i[1]] = stats[i[0]].get(i[1], 0) + 1
    return stats


def compo_stats(stats):
    """argument: dictionary of counters
    return: dictionary of percentages"""
    tot = sum(list(stats.values()))
    percent = dict()
    for i in stats.keys():
        percent[i] = stats[i]/tot * 100
    return percent


def ssres_stats(stats):
    """argument: dictionary of counters
    return: dictionary of dictionaries of percentages
    the dictionary format is {ss:{aa:counter}}
    Compute the percentage of each amino acid with respect to all aa in a certain ss"""
    percent = dict()
    for i in stats.keys():
        percent[i] = compo_stats(stats[i])
    return percent


def dictpprint(dic):
    """pretty print"""
    for k in sorted(dic.keys()):
        print(k, dic[k])


def inheat(rows, cols):
    """1st argument: iterable w/keys of rows
    2nd argument: iterable w/ keys of columns
    Create a dictionary of lists of counters initalized to zero"""
    heatmap = dict()
    for i in rows:
        heatmap[i] = [0 for j in cols]
    return heatmap

def slidingwindow(aa, ss, heatmap, desired_ss):
    """1st argument: amino acid seq
    2nd argument: dssp seq
    3rd argument: dictionary of lists of counters
    4th argument specifying whether we are doing the heatmap for helices or strands
    return: updated dictionary"""
    # this function shows a somewhat daunting bug which I will not fix now since it's 2am
    # FIXED: the bug was due to the implementation of inheat
    for i in range(8, len(aa)-8):
        if ss[i] == desired_ss:
            w_ptr = 0
            seq_ptr = i - 8
            while w_ptr < 17:
                heatmap[aa[seq_ptr]][w_ptr] = heatmap.get(aa[seq_ptr])[w_ptr]+ 1
                w_ptr += 1
                seq_ptr += 1
    return heatmap


if __name__ == "__main__":
    # test code usage:
    # first argument = file with ids
    # second argument = folder w/ fasta
    # third argument = folder w/ dssp

    ids = open(sys.argv[1])
    overall = dict()
    ss_compo = dict()
    aa_compo = dict()
    heatmapH = inheat(amino_acids, range(1,18))
    heatmapE = inheat(amino_acids, range(1,18))

    for id in ids:
        # we know that each file only has two lines, the first one is always the header
        # the second one is always the sequence we are interested in
        id = id.rstrip()
        fasta_path = open(sys.argv[2]+id+'.fasta')
        fasta = fasta_path.readlines()[1].rstrip()
        dssp_path = open(sys.argv[3]+id+'.dssp')
        dssp = dssp_path.readlines()[1].rstrip()
        ss_compo = compo(dssp, ss_compo)
        aa_compo = compo(fasta, aa_compo)
        overall = ssperresidue(fasta, dssp, overall)
        heatmapH = slidingwindow(fasta, dssp, heatmapH, 'H')
        heatmapE = slidingwindow(fasta, dssp, heatmapE, 'E')
        fasta_path.close()
        dssp_path.close()
#    dictpprint(ss_compo)
#    dictpprint(aa_compo)
#    dictpprint(overall)
#    print(compo_stats(aa_compo))
#    dictpprint(ssres_stats(overall))
#    dictpprint(heatmapH)
#    dictpprint(heatmapE)
    ids.close()

# all checks ok
