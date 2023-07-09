import random
from Bio import SeqIO 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import *
Nucleotides= ["A","C","G","T"]

randomStr=''.join([random.choice(Nucleotides)
                   for nuc in range(21)])




answer = input('which sequece you will type?\n ' + "(dna or 1)  ")

if answer.upper() in ["DNA","1"]:
    parse_or_seq = input("file or seq: ")
    if parse_or_seq.lower() in ["file","1"]:
        file_path= input("file path: ")
        file_type= input("file type: ")
        for record in SeqIO.parse(file_path, file_type):
            seq = record.seq
            print(f"id = {record.description} \n your sequence is: {seq}")
    elif parse_or_seq.lower() in ["seq","2"]:
        sequence = randomStr
        seq = Seq(sequence)
        print(f"sequence is: {seq}")
    while True:
        next_step = input(
    """what's next: 
    1) transcribe   
    2) translate  
    3)complement seq(comp) 
    4)reverse_complement(rev-comp) 
    5)reverse_complement_rna(5) 
    6)GC content(%) 
    7)Molecular weight(mw)
    8)six frame translation
    __stop or x:  """)
        if next_step.lower() in  ['transcribe' ,'1']:
            print(f"m-RNA is: {seq.transcribe()}")
            records= SeqRecord(seq.transcribe(),id='',name='transcribe',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        elif next_step.lower() in ['translate','2']:
                print(f"Amino acid seq: {seq.translate()}") 
                records= SeqRecord(seq.translate(),id='',name='translate',description=f'{record.description}')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")
        elif   next_step.lower() in ['comp','3'] :
                print(f"complement seq: {seq.complement()}")
                records= SeqRecord(seq.complement(),id='',name='complementary sequence',description=f'{record.description}')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")   
        elif next_step.lower() in ['rev-comp','4']:
            print(f"{seq.reverse_complement()}")
            records= SeqRecord(seq.reverse_complement(),id='',name='reverse complementary sequence',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta") 
        elif next_step.lower() in ['comp_rna','5']:
            print(f"complment rna:{seq.complement_rna()}")   
            records= SeqRecord(seq.complement_rna(),id='',name='reverse complementary_rna sequence',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        elif next_step.lower() in  ['gc','6']:
            print(f"GC %: {round(gc_fraction(seq) * 100) }")
        elif next_step.lower() in ['mw', '7']:
            print(f'Moleceular weight: {molecular_weight(seq)}')
        elif next_step.lower() in ['sixframe', '8']:
            print(f'{six_frame_translations(seq) }')     
            name = input("name of file:")
            test = open(name+".txt",'w')
            test.write(str(six_frame_translations(seq) ))
            test.close()
        elif next_step.lower() in  ['stop' ,'x']:
            break
else:
        pass        



