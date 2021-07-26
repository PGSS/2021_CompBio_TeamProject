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

        #fidning the sight of where the forward primer is bound
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

from Bio import Restriction

from Bio.Restriction import EcoRI

dir()
['Restriction', '__builtins__', '__doc__', '__name__', '__package__']
Restriction.EcoRI
EcoRI
Restriction.EcoRI.site
'GAATTC'


from Bio.Seq import Seq
my_seq = Seq('AAAAAAAAAAAAAA')
my_seq
Seq('AAAAAAAAAAAAAA')

EcoRI.search(my_seq)


# Retrieving the sequences produced by a digestion

ecoseq = my_seq + Seq(EcoRI.site) + my_seq
ecoseq
Seq('AAAAAAAAAAAAAAGAATTCAAAAAAAAAAAAAA')
EcoRI.search(ecoseq)
[16]

print(ecoseq[:15], ecoseq[15:])
print(EcoRI.catalyze(ecoseq))
print(EcoRI.search(ecoseq, linear=False))
#[16]
EcoRI.catalyse(ecoseq, linear=False)
#(Seq('AATTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAG'))
ecoseq  # for memory
#Seq('AAAAAAAAAAAAAAGAATTCAAAAAAAAAAAAAA')


EcoRI.search(ecoseq, linear=False)
#[16]
EcoRI.catalyse(ecoseq, linear=False)
#(Seq('AATTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAG'),)
ecoseq  # for memory

# circular sequences
new_seq = Seq('TTCAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAA')
EcoRI.search(new_seq)
EcoRI.search(new_seq, linear=False)

#gies enqyme seq








