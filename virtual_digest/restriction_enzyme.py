#get a set of primers
    #for each bacteria record
    #get PCR out of the record
    #take exyme in, use chunk and run the virtual enqyme digest

    #METHOD
        #Func: takes PCR product and enzyme and gets chunks
        #func: takes a full record and gives back the PCR product

        #two index of forward and reverse
        #take everything between on the forward to reverse
        #slice and take that region out

        #finding the sight of where the forward primer is bound
        #gives a start index = this site
        #reverse = STOP index
        #PCR pro

# database group, run through data
    # get a record
    # list of bacteria
    # use as input, take forward and reverse seq
    # (dont need0 do PCR on each
    # output: seq, or nothing
    #do virtual digest on PCR


#loading in the enzyme seq


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

'''
for seq_record in SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
'''

counter = 1
for seq_record in bacteria_database:
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))
    print()
    counter += 1
    if(counter > 10):
        break
#iterate through each sequence in database
#make output object w/ bacteria name, number, & dictionary of fragment lengths


