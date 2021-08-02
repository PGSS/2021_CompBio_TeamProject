from types import List


from Bio.Seq import Seq
my_seq = Seq('AAAAAAAAAAAAAA')

Seq('AAAAAAAAAAAAAAGAATTCAAAAAAAAAAAAAA')

# creating RestrictionBatch
import Bio
from Bio import Restriction
from Bio.Restriction import *
rb = RestrictionBatch(['AluI', 'HaeIII', 'MboI'])

# analyze sequence w/ RestrictionBatch
sequence = Seq('''TCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATG
CCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTA
AAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCC
AGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCA
CGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTG
ATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGC
GTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAAACTGACGCTCA
GGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGA
TGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGAC
CGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAG
CGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCC
ACAGAACCTTGTAGAGATACGAGGGTGCCTTCGGGAACTGTGAGACAGGTGCTGCATGGC
TGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATC
CTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGA
AGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAA
TGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTA
GTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATC
AGAATGTCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAG
TGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTC 
ATGACTGGGGTGAAGTCGTAACAAGGTAACC''')
sequence = sequence.replace(" ","")
sequence = sequence.replace("\n", "")
seq_dict = rb.search(sequence)
#print(seq_dict)

def find_Lengths(sequence, seq_dict):
    #create new dictionary with same keys as seq_dict
    length_dict = {} #empty dictionary
    for key in seq_dict:
        lengths = []
        for i in range(len(seq_dict[key])):
            if i==0:
                lengths.append(seq_dict[key][0])
            else:
                lengths.append(seq_dict[key][i]-seq_dict[key][i-1])
                if i==len(seq_dict[key])-1:
                    lengths.append(len(sequence) - seq_dict[key][i])
        length_dict[key] = lengths
    return length_dict

#print(find_Lengths(sequence,seq_dict))

import Bio
from Bio import SeqIO
from os.path import expanduser

home = expanduser("~")
bacteria_database = SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta")

#make a class containing bacteria name, id, and length dict
class BacteriaInfo:

    #initializer
    def __init__(self, name, id, length_dict):
        self.name = name
        self.id = id
        self.lengths = length_dict

    #change representation to string
    def __str__(self):
        return f"{{name: {self.name}, id: {self.id}, lengths: {self.lengths}}}\n"
    def __repr__(self):
        return self.__str__()

seq_name = []


import pandas as pd
import csv
home = expanduser("~")

#col_list = ['Bacteria Name', 'Bacteria Sequence File Name', 'Length of PCR Product', 'PCR Product', 'AluI', 'HaeIII', 'MboI']
df = pd.read_csv(home + "/test_data_bacteria.csv", header=0)

def create_dict_from_row(row):
    test = {'AluI': [], 'HaeIII': [], 'MboI': []}
    test['AluI'] = [int(x) for x in row['AluI'].split(',')]
    test['HaeIII'] = [int(x) for x in row['HaeIII'].split(',')]
    test['MboI'] = [int(x) for x in row['MboI'].split(',')]

    return test

test_data = []

for idx, row in df.iterrows():
    length_dict = create_dict_from_row(row)
    name = row['Bacteria Name']
    id = name
    test_data.append(BacteriaInfo(name, id, length_dict))

print(test_data)

for seq_record in bacteria_database:
    seq_name.append(seq_record)
    sequence = seq_record.seq
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    seq_dict = rb.search(sequence)
    length_dict = find_Lengths(sequence, seq_dict)
    bac = BacteriaInfo(seq_record.description, seq_record.id, length_dict)

counter = 0

# Fragment Comparer
def fragment_comparer(candidate_list : List[int], lab_list : List[int]) :

    #exclude smaller than Alejandro's number: 100, 50
    #exclude if it [] and other has

    if
    #within 10%

    return False

# input: takes in list of bacteriaInfo objects (database)
#       takes in a single bacteriaInfo object (lab data)

#returns: list of bacteriaInfo objects from database that match lab data (empty-#)

def pos_reducer(database_list : List[BacteriaInfo], lab_data : BacteriaInfo):

    possible_matches = []

    for bacteria in database_list:
        haeIII = fragment_comparer(bacteria.lengths["HaeIII"], lab_data.lengths["HaeIII"])
        mboI = fragment_comparer(bacteria.lengths["MboI"], lab_data.lengths["MboI"])
        aluI = fragment_comparer(bacteria.lengths["AluI"], lab_data.lengths["AluI"])

        if haeIII and mboI and aluI:
            possible_matches.append(bacteria)

    return possible_matches

#pos_reducer()
#print(counter)






