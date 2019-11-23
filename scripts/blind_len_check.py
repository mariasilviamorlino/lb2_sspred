'''Aim:
for the whole set of blind candidates:
take the fasta file and the dssp file and check if the length of the sequences is the same
generate a report including all sequences passing the test
Usage: python3 blind_len_check.py <ids_to_check> <folder_to_fasta> <folder_to_dssp> <report_output>'''

import sys

# sys.argv[1]: file with PDB_chain ids   /home/mary/bioinfo/LB2/project/new_blind_testset/blind_candidates.txt
# sys.argv[2]: folder to fasta files     /home/mary/bioinfo/LB2/project/blind_testset/blind_fasta
# sys.argv[3]: folder to dssp files      /home/mary/bioinfo/LB2/project/blind_testset/blind_dssp_formatted
# sys.argv[4]: output file               /home/mary/bioinfo/LB2/project/blind_testset/blind_good_candidates.txt

with open(sys.argv[1]) as id_file:
    id_list = id_file.readlines()

fastapath = sys.argv[2]
dssppath = sys.argv[3]

with open(sys.argv[4], 'w') as outfile:
    for i in id_list:
        i = i.strip()
        fastafile = fastapath + '/' + i + '.fasta'
        dsspfile = dssppath + '/' + i + '.dssp'
        with open(fastafile) as fastahandle:
            fasta = fastahandle.readlines()[1].strip()
        with open(dsspfile) as dssphandle:
            dssp = dssphandle.readlines()[1].strip()
        
        if len(fasta) == len(dssp):
            outfile.write(i+'\n')
