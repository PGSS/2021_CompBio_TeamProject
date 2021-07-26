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