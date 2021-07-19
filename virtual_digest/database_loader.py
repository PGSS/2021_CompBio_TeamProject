import Bio
from Bio import SeqIO
from os.path import expanduser

home = expanduser("~")

for record in SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta"):
    print(record.id)