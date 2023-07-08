import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Nucleotides= ["A","C","G","T"]

randomStr=''.join([random.choice(Nucleotides)
                   for nuc in range(21)])




answer = input('which sequece you will type?\n ' + "(dna or 1)  ")

if answer.upper() in ["DNA","1"]:
    parse_or_seq = input("file or seq: ")
    if parse_or_seq.lower() == "file":
        file_path= input("file path: ")
        file_type= input("file type: ")
        for record in SeqIO.parse(file_path, file_type):
            seq = record.seq
            print(f"id = {record.description} \n your sequence is: {seq}")
    elif parse_or_seq.lower() == "seq":
        sequence = randomStr
        seq = Seq(sequence)
        print(f"sequence is: {seq}")
    while True:
        next_step = input("""what's next: \n 1) transcribe  \n 2) translate \n 3)complement seq(comp) \n 4)reverse_complement(rev-comp)\n 5)reverse_complement_rna(5) \n __stop or x:  """)
        if next_step.lower() in  ['transcribe' ,'1']:
            print(f"m-RNA is: {seq.transcribe()}")
            seq=seq.transcribe()
            records= SeqRecord(seq,id='',name='transcribe',description=f'{record.description}')

            with open("example.fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        elif next_step.lower() in ['translate','2']:
                print(f"Amino acid seq: {seq.translate()}") 
                seq=seq.translate()
                records= SeqRecord(seq,id='',name='translate',description=f'{record.description}')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:   
                   SeqIO.write(records, handle, "fasta")
        elif   next_step.lower() in ['comp','3'] :
                print(f"complement seq: {seq.complement()}")
                seq=seq.complement()
                records= SeqRecord(seq,id='',name='complementary sequence',description='complementary sequence')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")   
        elif next_step.lower() in ['rev-comp','4']:
            print(f"{seq.reverse_complement()}")
            seq=seq.reverse_complement()
            records= SeqRecord(seq,id='',name='reverse complementary sequence',description='reverse complementary sequence')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta") 
        elif next_step.lower() in ['rev_comp_rna','5']:
            seq=seq.complement_rna()
            print(f"complment rna:{seq}")   
            records= SeqRecord(seq,id='',name='reverse complementary_rna sequence',description='reverse complementary_rna sequence')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")         
        elif next_step.lower() in  ['stop' ,'x']:
            break
else:
        pass        



