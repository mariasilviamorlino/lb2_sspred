import sys


def taxonomy_parse(line):
    """argument: line of csv file
    return: PDB id and """
    l = line.rstrip().split(',')
    return (l[0].strip('"'), l[2].strip('"'))


def taxonomy_stats(id_taxo, stats, redundancy):
    if redundancy.get(id_taxo[0]) == None:
        stats[id_taxo[1]] = stats.get(id_taxo[1], 0) + 1
        redundancy[id_taxo[0]] = 'foo'
    return stats, redundancy


if __name__ == '__main__':
    # test code usage:
    # 1st argument = path to PDB--taxonomy mapping
    f = open(sys.argv[1])
    taxostats = dict()
    red = dict()
    for line in f:
        if line.startswith('"'):
            taxostats, red = taxonomy_stats(taxonomy_parse(line), taxostats, red)
    print(taxostats)
