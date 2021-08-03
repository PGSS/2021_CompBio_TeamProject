import Bio
from Bio import SeqIO
from os.path import expanduser

from Bio import Restriction
from Bio.Restriction import *

import pandas as pd

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

def create_dict_from_row(row):
    test = {'AluI': [], 'HaeIII': [], 'MboI': []}
    test['AluI'] = [int(x) for x in row['AluI'].split(',')]
    test['HaeIII'] = [int(x) for x in row['HaeIII'].split(',')]
    test['MboI'] = [int(x) for x in row['MboI'].split(',')]
    return test

MINIMUM_BASE_PAIR = 50
###################################################
# Fragment Comparer
def fragment_comparer(candidate_list, lab_list): #: List[int], lab_list : List[int]):

    #exclude smaller than Alejandro's number: 100, 50
    #exclude if it [] and other has

    if (len(candidate_list) == 0 and len(lab_list) != 0):
        print("one dooesn't digest", candidate_list, lab_list)
        return False

    filtered_candidates = [x for x in candidate_list if x > MINIMUM_BASE_PAIR]
    filtered_candidates.sort()

    filtered_lab = [x for x in lab_list if x > MINIMUM_BASE_PAIR]
    filtered_lab.sort()

    if len(filtered_candidates) != len(filtered_lab):
        print("different products after filter", filtered_candidates, filtered_lab)
        return False

    #1: comparing each number to see if it fit within range
    for i in range(len(filtered_candidates)):
        if (abs(filtered_candidates[i]-filtered_lab[i]) > 0.1 * filtered_lab[i]):
            print("invalid comparison", filtered_candidates[i], filtered_lab[i], filtered_candidates, filtered_lab)
            return False

    return True

##############################################################################################
# input: takes in lis
# t of bacteriaInfo objects (database)
#       takes in a single bacteriaInfo object (lab data)
#returns: list of bacteriaInfo objects from database that match lab data (empty-#)


def pos_reducer(database_list, lab_data): #: List[BacteriaInfo], lab_data : BacteriaInfo):

    possible_matches = []
    processed = 0
    for bacteria in database_list:
        for result in lab_data:
            haeIII = fragment_comparer(bacteria.lengths[HaeIII], result.lengths['HaeIII'])
            mboI = fragment_comparer(bacteria.lengths[MboI], result.lengths['MboI'])
            aluI = fragment_comparer(bacteria.lengths[AluI], result.lengths['AluI'])
            if haeIII and mboI and aluI:
                print(bacteria)
                possible_matches.append(bacteria)
            processed += 1
            if processed % 1000:
                print(f"proceseed {processed} records, found: {len(possible_matches)}, matches")
    return possible_matches



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
home = expanduser("~")
#bacteria_database = SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta")
count = 0
#record in bacteria_database:
#    count = count + 1
#    if count == 1:
#        print(seq_record)
#        break
#print(count)
#count = 0
#print(seq_record.seq)

df = pd.read_csv(home + "/test_data_bacteria.csv", header=0)
test_data = []

for idx, row in df.iterrows():
    length_dict = create_dict_from_row(row)
    name = row['Bacteria Name']
    id = name
    test_data.append(BacteriaInfo(name, id, length_dict))

print(test_data[0])

rb = RestrictionBatch(['AluI', 'HaeIII', 'MboI'])
'''
print("starting database analysis: ")

ecoli_seq = []

for seq_record in bacteria_database:
    if 'Escherichia coli' in seq_record.description:
        ecoli_seq.append(seq_record)
        if 'TCCTACGGGAGGCAGCAGT' in seq_record.seq:
            print("contains forward primer", seq_record.id)

        print(seq_record.description)
        sequence = seq_record.seq
        sequence = sequence.replace(" ", "")
        sequence = sequence.replace("\n", "")
        seq_dict = rb.search(sequence)
        length_dict = find_Lengths(sequence, seq_dict)
        print(length_dict)

print("ecoli sequences", len(ecoli_seq))
SeqIO.write(ecoli_seq, home+'/ecoli.fasta', 'fasta')
exit(0)
'''

#manually making test records
test_record = {}


bac_list = []
possible_matches = []
seq_count = 0
match_count = 0

test_matches = {}
for bac in test_data:
    test_matches[bac.name] = []

ecoli_database = list(SeqIO.parse(home+'/ecoli.fasta','fasta'))
print(len(ecoli_database))
for seq_record in ecoli_database:
    if (seq_count==0):
        print(seq_record)

    sequence = seq_record.seq
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    seq_dict = rb.search(sequence)
    #print(seq_dict)
    length_dict = find_Lengths(sequence, seq_dict)
    #print(length_dict)
    bac = BacteriaInfo(seq_record.description, seq_record.id, length_dict)
    #bac_list.append(bac)
    #break

    for result in [test_data[0]]:

        haeIII = fragment_comparer(bac.lengths[HaeIII], result.lengths['HaeIII'])
        mboI = fragment_comparer(bac.lengths[MboI], result.lengths['MboI'])
        aluI = fragment_comparer(bac.lengths[AluI], result.lengths['AluI'])
        if haeIII and mboI and aluI:
            match_count = match_count + 1
            test_matches[result.name].append(bac)
    seq_count = seq_count + 1
    if seq_count % 10000 == 0:
        print('in ',seq_count,' sequences, found ',match_count, ' matches')

print(test_matches)
#print('In ', seq_count, 'sequences ', match_count, 'matched')
#print(possible_matches)


rb = RestrictionBatch(['AluI', 'HaeIII', 'MboI'])

'''
bac_list = []
possible_matches = []
seq_count = 0
match_count = 0
for seq_record in bacteria_database:
    sequence = seq_record.seq
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    seq_dict = rb.search(sequence)
    #print(seq_dict)
    length_dict = find_Lengths(sequence, seq_dict)
    #print(length_dict)
    bac = BacteriaInfo(seq_record.description, seq_record.id, length_dict)
    #bac_list.append(bac)
    #break
    for result in test_data:
        haeIII = fragment_comparer(bac.lengths[HaeIII], result.lengths['HaeIII'])
        mboI = fragment_comparer(bac.lengths[MboI], result.lengths['MboI'])
        aluI = fragment_comparer(bac.lengths[AluI], result.lengths['AluI'])
        if haeIII and mboI and aluI:
            match_count = match_count + 1
            possible_matches.append(bac)
    seq_count = seq_count + 1
    if seq_count % 10000 == 0:
        print('in ',seq_count,' sequences, found ',match_count, ' matches')

    if seq_count == 20000:
        for members in possible_matches:
            print(members)

print('In ', seq_count, 'sequences ', match_count, 'matched')
print(possible_matches)
'''
#how to not enter the douplicate sequences





