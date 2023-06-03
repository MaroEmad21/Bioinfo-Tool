import random
from Bio import SeqIO
from Bio.Seq import Seq
Nucleotides= ["A","C","G","T"]

randomStr=''.join([random.choice(Nucleotides)
                   for nuc in range(20)])




answer = input('which sequece you will type?\n ' + "(dna , rna, aa)  ")
def main():
    if answer.upper() == "DNA":
        parse_or_seq = input("file or seq: ")
        if parse_or_seq.lower() == "file":
            file_path= input("file path: ")
            file_type= input("file type: ")
            for record in SeqIO.parse(file_path, file_type):
                seq = record.seq
                print(f"your sequence is: {seq}")
        elif parse_or_seq.lower() == "seq":
            sequence = randomStr
            seq = Seq(sequence)
            print(f"sequence is: {seq}")
        next_step = input("what's next: ")
        if next_step.lower() == 'transcribe':
                print(f"m-RNA is: {seq.transcribe()}")
        else: 
                print(" no choice")    
    else:
        pass    

main()            