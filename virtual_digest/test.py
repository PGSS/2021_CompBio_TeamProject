##########################################################################################################
#PROJECT TITLE: Microbial Forensics: Identifying Bacteria and Yeast Using Ribosomal DNA Fingerprints
#Author: Hannah Chang and Claire Zheng
#Instrucor: Mr. Andrew McGuier
#PGSS COMPUTATIONAL BIOLOGLY

###################################################
#IMPORTING PACKAGES

import Bio
from Bio import SeqIO
from os.path import expanduser

from Bio import Restriction
from Bio.Restriction import *

import pandas as pd


###################################################
#Function find lengths calculates fragment lengths for bacteria DNA ribosmal sequences given cut sights
#inputs: takes in a sequences
#returns: a dictionary of fragment lengths

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

###################################################
#bacteria class contains bacteria name, id, and length dict
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

###################################################
# function creates a dictionary with enzymes set as keys
# inputs: takes in a row and adds creates a list of fragment legnths from each enzyme
# returns: the dictionary of fragment lengths
def create_dict_from_row(row):
    test = {'AluI': [], 'HaeIII': [], 'MboI': []}
    test['AluI'] = [int(x) for x in row['AluI'].split(',')]
    test['HaeIII'] = [int(x) for x in row['HaeIII'].split(',')]
    test['MboI'] = [int(x) for x in row['MboI'].split(',')]
    return test

MINIMUM_BASE_PAIR = 100

###################################################
# function compares two lists of fragment lengths
# returns: true (lists of fragment lengths match) or false (lists of fragment lengths do not match)
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
        return False

    #1: comparing each number to see if it fit within range
    for i in range(len(filtered_candidates)):
        if (abs(filtered_candidates[i]-filtered_lab[i]) > 0.1 * filtered_lab[i]):
            return False

    return True

###################################################
# LOADING IN TEST DATA AND DATABASE
home = expanduser("~")
#bacteria_database = SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta")
count = 0

#bacteria_database = list(SeqIO.parse(home+'/pcr_product.fa','fasta'))

#creating a bacteria objects for each bacteria in test data
test_data = []
df = pd.read_csv(home + "/test_data_bacteria.csv", header=0)
for idx, row in df.iterrows():
    length_dict = create_dict_from_row(row)
    name = row['Bacteria Name']
    id = name
    test_data.append(BacteriaInfo(name, id, length_dict))


#################################################   MAIN PROGRAM    ####################################################

rb = RestrictionBatch(['AluI', 'HaeIII', 'MboI'])

print("starting database analysis: ")

match_count = 0
possible_matches = []
pcr_count = 0
pcr_database = {}
test_bac_registered = False
test_bac_list = []
test_matches = {}
test_match_count = {}

for result in test_data:
    test_match_count[result.name] = 0

pcr_database = list(SeqIO.parse(home + '/pcr_product.fa', 'fasta'))

#iterate through the database and create bacteria objects
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

    #iterate through the test dataset and compares fragment lengths from each enzyme with the database
    for result in test_data:
        test_count = test_count + 1
        haeIII = fragment_comparer(pcr.lengths[HaeIII], result.lengths['HaeIII'])
        mboI = fragment_comparer(pcr.lengths[MboI], result.lengths['MboI'])
        aluI = fragment_comparer(pcr.lengths[AluI], result.lengths['AluI'])

        test_bac_registered = False

        #sequences are a match of fragment lengths in all three enzymes match
        if haeIII and mboI and aluI:
            match_count = match_count + 1
            possible_matches.append(seq_record)

            if match_count == 1:
                test_bac_registered = False
            else:
                for element in test_bac_list:
                    if element.name == result.name:
                        test_bac_registered = True
                        test_match_count[result.name] = test_match_count[result.name] + 1

            if test_bac_registered == False:
                test_bac_list.append(result)
                test_match_count[result.name] = test_match_count[result.name] + 1

    if pcr_count % 100000 == 0:
        print('in ', pcr_count, ' PCR sequences, found ', match_count, ' matches')

# PRINTING RESULTS
print(test_bac_list)
print(test_match_count)


