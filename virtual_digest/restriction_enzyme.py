from Bio.Seq import Seq

# creating RestrictionBatch

import Bio
from Bio import SeqIO
from os.path import expanduser

from Bio import Restriction
from Bio.Restriction import *
rb = RestrictionBatch(['AluI', 'HaeIII', 'MboI'])

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


#-----------------------------------------------
#test data
import pandas as pd
home = expanduser("~")
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

#print(test_data)

#------------------------------------------------------

#bacteria in the database
bac_list = []
for seq_record in bacteria_database:
    #seq_name.append(seq_record)
    sequence = seq_record.seq
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    seq_dict = rb.search(sequence)
    length_dict = find_Lengths(sequence, seq_dict)
    bac = BacteriaInfo(seq_record.description, seq_record.id, length_dict)
    bac_list.append(bac)
    break
print(bac_list)

###################################################
# Fragment Comparer
def fragment_comparer(candidate_list, lab_list): #: List[int], lab_list : List[int]):

    #exclude smaller than Alejandro's number: 100, 50
    #exclude if it [] and other has

    num_true = 0
    #1: comparing each number to see if it fit within range
    for i in range(len(candidate_list)):
        for j in range(len(lab_list)):
            if (abs(candidate_list[i]-lab_list[j]) <= 15):
                num_true += 1
                break
        #stats
        if (num_true/len(candidate_list) >= 0.9):
            return True
    return False

##############################################################################################
# input: takes in list of bacteriaInfo objects (database)
#       takes in a single bacteriaInfo object (lab data)
#returns: list of bacteriaInfo objects from database that match lab data (empty-#)


def pos_reducer(database_list, lab_data): #: List[BacteriaInfo], lab_data : BacteriaInfo):

    possible_matches = []

    for bacteria in database_list:
        #haeIII = fragment_comparer(bacteria.lengths["HaeIII"], lab_data.lengths["HaeIII"])
        mboI = fragment_comparer(bacteria.lengths["MboI"], lab_data.lengths["MboI"])
        aluI = fragment_comparer(bacteria.lengths["AluI"], lab_data.lengths["AluI"])

        #print(haeIII)
        #print(bacteria.lengths)

        if haeIII and mboI and aluI:
            possible_matches.append(bacteria)

    return possible_matches

print(pos_reducer(bac_list,test_data))




