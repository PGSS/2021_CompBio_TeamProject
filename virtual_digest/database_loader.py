import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from os.path import expanduser

home = expanduser("~")

bacteria_database = SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta")

forward_primer = Seq('TCCTACGGGAGGCAGCAGT'.lower())
reverse_primer = Seq('GGTTACCTTGTTACGACTT'.lower()).reverse_complement()


processed = 0
filtered_bacteria = []
for seq in bacteria_database:
    processed += 1
    if processed % 1000 == 0:
        print(f"Processed {processed} records from the database"
              f"")
    forward_location = seq.seq.find(forward_primer)
    if forward_location < 0:
        continue
    reverse_location = seq.seq.find(reverse_primer)
    filtered_seq = None
    if reverse_location >= 0:
        filtered_seq = seq.seq[forward_location:reverse_location]
        filtered_seq = filtered_seq[len(forward_primer) :] if filtered_seq.startswith(forward_primer) else filtered_seq
    else:
        filtered_seq = seq.seq[forward_location:]
        filtered_seq = filtered_seq[len(forward_primer):] if filtered_seq.startswith(forward_primer) else filtered_seq

    seq.seq = filtered_seq
    filtered_bacteria.append(seq)
SeqIO.write(filtered_bacteria, home+'/pcr_product.fa', 'fasta')
print("original bacteria database length", processed)
print("Filtered length:",len(filtered_bacteria))




'''
for seq_record in SeqIO.parse(home + "/current_Bacteria_unaligned.fa", "fasta"):
    print(seq_record.id)
    print(repr(seq_record.seq))
    print(len(seq_record))

'''