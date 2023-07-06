import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Nucleotides= ["A","C","G","T"]

randomStr=''.join([random.choice(Nucleotides)
                   for nuc in range(20)])




answer = input('which sequece you will type?\n ' + "(dna or 1)  ")

if answer.upper() in ["DNA","1"]:
    parse_or_seq = input("file or seq: ")
    if parse_or_seq.lower() == "file":
        file_path= input("file path: ")
        file_type= input("file type: ")
        for record in SeqIO.parse(file_path, file_type):
            seq = record.seq
            print(f"id = {record.id} \n your sequence is: {seq}")
    elif parse_or_seq.lower() == "seq":
        sequence = randomStr
        seq = Seq(sequence)
        print(f"sequence is: {seq}")
    """next_step = input("""#what's next: \n 1) transcribe  \n 2) translate
                            ##""")"""
    while True:
        next_step = input("""what's next: \n 1) transcribe  \n 2) translate \n )stop or x:  """)
        if next_step.lower() in  ['transcribe' ,'1']:
            print(f"m-RNA is: {seq.transcribe()}")
            seq=seq.transcribe()
            records= SeqRecord(seq,id='1.0',name='transcribe',description='haemoglobin Rna')

            with open("example.fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        elif next_step.lower() in ['translate','2']:
                print(f"Amino acid seq: {seq.translate()}") 
                seq=seq.translate()
                records= SeqRecord(seq,id='',name='translate',description='haemoglobin AA')

                with open("example2.fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")
        elif next_step.lower() in  ['stop' ,'x']:
            break
else:
        pass        


