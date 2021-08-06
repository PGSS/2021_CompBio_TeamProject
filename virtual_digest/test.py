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
                lengths.append(seq_dict[key][0]-1)
            else:
                lengths.append(seq_dict[key][i]-seq_dict[key][i-1])
                if i==len(seq_dict[key])-1:
                    lengths.append(len(sequence) - seq_dict[key][i] +1)
        length_dict[key] = lengths
    return length_dict

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

def create_dict_from_row(row):
    test = {'AluI': [], 'HaeIII': [], 'MboI': []}
    test['AluI'] = [int(x) for x in row['AluI'].split(',')]
    test['HaeIII'] = [int(x) for x in row['HaeIII'].split(',')]
    test['MboI'] = [int(x) for x in row['MboI'].split(',')]
    return test

test_data = []
home = expanduser("~")
bacteria_database = SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta")
#df = pd.read_csv(home + "/test_bacteria_data.csv", header=0)
df = pd.read_csv(home + "/real_lab_data.csv", header=0)
#df_new = df[df['AluI'].notnull()]

for idx, row in df.iterrows():
    length_dict = create_dict_from_row(row)
    name = row['Bacteria Name']
    id = name
    test_data.append(BacteriaInfo(name, id, length_dict))

MINIMUM_BASE_PAIR = 100
###################################################
# Fragment Comparer
def fragment_comparer(candidate_list, lab_list): #: List[int], lab_list : List[int]):

    #exclude smaller than Alejandro's number: 100, 50
    #exclude if it [] and other has

    if (len(candidate_list) == 0 and len(lab_list) != 0):
        #print("one dooesn't digest", candidate_list, lab_list)
        return False

    filtered_candidates = [x for x in candidate_list if x > MINIMUM_BASE_PAIR]
    filtered_candidates.sort()

    filtered_lab = [x for x in lab_list if x > MINIMUM_BASE_PAIR]
    filtered_lab.sort()

    if len(filtered_candidates) != len(filtered_lab):
        #print("different products after filter", filtered_candidates, filtered_lab)
        return False

    #1: comparing each number to see if it fit within range
    for i in range(len(filtered_candidates)):
        if (abs(filtered_candidates[i]-filtered_lab[i]) > 0.1 * filtered_lab[i]):
            #print("invalid comparison", filtered_candidates[i], filtered_lab[i], filtered_candidates, filtered_lab)
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





home = expanduser("~")
bacteria_database = SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta")
count = 0
#record in bacteria_database:
#    count = count + 1
#    if count == 1:
#        print(seq_record)
#        break
#print(count)
#count = 0
#print(seq_record.seq)

#df = pd.read_csv(home + "/test_data_bacteria.csv", header=0)

bacteria_database = list(SeqIO.parse(home+'/pcr_product.fa','fasta'))

test_data = []

#df = pd.read_csv(home + "/real_lab_data.csv", header=0)
df = pd.read_csv(home + "/test_data_bacteria.csv", header=0)
#df_new = df[df['AluI'].notnull()]

for idx, row in df.iterrows():
    length_dict = create_dict_from_row(row)
    name = row['Bacteria Name']
    id = name
    test_data.append(BacteriaInfo(name, id, length_dict))

#print(test_data[0])

rb = RestrictionBatch(['AluI', 'HaeIII', 'MboI'])

#print("starting database analysis: ")

ecoli_seq = []

#manually making test records
test_record = {}


print("starting database analysis: ")
pcr_seq = []
bac_list = []
possible_matches = []
seq_count = 0
match_count = 0
test_matches = {}
pcr_count = 0
pcr_database = {}
test_bac_found = False
test_bac_list = []
test_match_count = {}

for result in test_data:
    test_match_count[result.name] = 0

# '''
# check if these are in the database

count = 0

pcr_database = list(SeqIO.parse(home + '/pcr_product.fa', 'fasta'))

print('Total PCR sequences: ', len(pcr_database))
for seq_record in pcr_database:
    pcr_count = pcr_count + 1
    sequence = seq_record.seq
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    seq_dict = rb.search(sequence)
    length_dict = find_Lengths(sequence, seq_dict)
    pcr = BacteriaInfo(seq_record.description, seq_record.id, length_dict)

    test_count = 0

    for result in test_data:
        test_count = test_count + 1
        #        print ('test ', test_count, ':', end='')
        #        print('HaeII: ', end='')
        haeIII = fragment_comparer(pcr.lengths[HaeIII], result.lengths['HaeIII'])
        #        print('         MboI: ', end='')
        mboI = fragment_comparer(pcr.lengths[MboI], result.lengths['MboI'])
        #        print(' AluI: ', end='')
        aluI = fragment_comparer(pcr.lengths[AluI], result.lengths['AluI'])

        test_bac_found = False
        if haeIII and mboI and aluI:
            match_count = match_count + 1
            possible_matches.append(seq_record)

            if match_count == 1:
                test_bac_found = False
            else:
                for element in test_bac_list:
                    if element.name == result.name:
                        test_bac_found = True
                        test_match_count[result.name] = test_match_count[result.name] + 1

            if test_bac_found == False:
                test_bac_list.append(result)
                test_match_count[result.name] = test_match_count[result.name] + 1

    if pcr_count % 100000 == 0:
        print('in ', pcr_count, ' PCR sequences, found ', match_count, ' matches')

# print(possible_matches)
print(test_bac_list)

total = 0
for match in test_data:
    total = total + test_match_count[match.name]
# print('total matches: ', total)
print(test_match_count)


rb = RestrictionBatch(['AluI', 'HaeIII', 'MboI'])
