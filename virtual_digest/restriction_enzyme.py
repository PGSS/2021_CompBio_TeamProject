#dir()
#['Restriction', '__builtins__', '__doc__', '__name__', '__package__']
#Restriction.EcoRI
#EcoRI
#Restriction.EcoRI.site
#'GAATTC'


from Bio.Seq import Seq
my_seq = Seq('AAAAAAAAAAAAAA')
#my_seq
#Seq('AAAAAAAAAAAAAA')

#EcoRI.search(my_seq)


# Retrieving the sequences produced by a digestion

#ecoseq = my_seq + Seq(EcoRI.site) + my_seq
#ecoseq
Seq('AAAAAAAAAAAAAAGAATTCAAAAAAAAAAAAAA')
#EcoRI.search(ecoseq)
[16]

'''print(ecoseq[:15], ecoseq[15:])
print(EcoRI.catalyze(ecoseq))
print(EcoRI.search(ecoseq, linear=False))
#[16]
EcoRI.catalyse(ecoseq, linear=False)
#(Seq('AATTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAG'))
ecoseq  # for memory
#Seq('AAAAAAAAAAAAAAGAATTCAAAAAAAAAAAAAA')'''


'''EcoRI.search(ecoseq, linear=False)
#[16]
EcoRI.catalyse(ecoseq, linear=False)
#(Seq('AATTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAG'),)
ecoseq  # for memory'''

''''# circular sequences
new_seq = Seq('TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA')
EcoRI.search(new_seq)
EcoRI.search(new_seq, linear=False)'''

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
print(seq_dict)

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

    #use for loop for each key in seq_dict
        #make new list for fragment length values
        #another for loop to subtract between elements in the list of cut sites

print(find_Lengths(sequence, seq_dict))
print(len(sequence))

#import database

import Bio
from Bio import SeqIO
from os.path import expanduser

home = expanduser("~")

bacteria_database = SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta")

#make a class containing bacteria name, id, and length dict
class BacteriaInfo:
    def __init__(self, name, id, length_dict):
        self.name = name
        self.id = id
        self.lengths = length_dict

counter = 1

seq_name = []

for seq_record in bacteria_database:

    seq_name.append(seq_record)

    sequence = seq_record.seq
    sequence = sequence.replace(" ", "")
    sequence = sequence.replace("\n", "")
    seq_dict = rb.search(sequence)
    length_dict = find_Lengths(sequence, seq_dict)
    bac = BacteriaInfo(seq_record.description, seq_record.id, length_dict)
        #print(bac.name)
        #print(bac.id)
        #print(bac.lengths)
        #print()

    counter += 1
    if (counter > 10):
        break

for i in range(0,len(seq_name)):
    print(seq_name[i])

'''
def pos_reducer(seq_name, seq):
    # iterate through each sequence in database

    for i in range(0, len(seq)-1):

        if(seq[i].lengths['fragment_sizes'] == seq_name[i].lengths['fragment_sizes']):

            for key in seq[i].lengths:
                if (seq[i].lengths[key] == seq_name[i].lengths[key]) :
                    seq[i].name = seq_name[i].name

    return seq

pos_reducer(seq_name, )
'''