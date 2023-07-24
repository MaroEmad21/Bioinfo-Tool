from Bio import AlignIO , SeqIO  ,Align
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import *
from Bio.Align import PairwiseAlignment, MultipleSeqAlignment , PairwiseAligner
import os

def make_alignment(sequence1):
    Q1= int(input("how many sequences? "))
    if Q1 == 2:
        aligner = PairwiseAligner()
        Q2 = input("mode (Global or local)")
        if Q2.lower() == 'global':
            aligner.mode = "global"
        else: 
            aligner.mode = "local"
        Q3 = input("with other  or 2 different  ")
        if Q3.lower() == "with":
            file_path= input("file path: ")
            file_type= input("file type: ")
            for record in SeqIO.parse(file_path, file_type):
                sequence2 = record.seq   
            alignment=aligner.align(sequence1,sequence2)
            print(alignment[0])
            print(f"score is : {alignment.score}")
            
            test=open("test2.aln", "Clustal")
            test.write(str(alignment[0]))
            test.close()             
        elif Q3.lower() == "diff":
            file_path= input("file path: ")
            file_type= input("file type: ")
            for record in SeqIO.parse(file_path, file_type):
                sequenceA = record.seq
            file_path= input("file path: ")
            file_type= input("file type: ")
            for record in SeqIO.parse(file_path, file_type):
                sequenceB = record.seq
            alignment=aligner.align(sequenceA,sequenceB)
            print(alignment[0])            
    elif Q1 >= 3:
            import glob
            path=input("path of files: ")
            files = []
            seqs=[]
            for file in glob.glob(f"{path}"+"\*.fasta"):
                files.append(file)
            for i in files:
                for records in  SeqIO.parse(f"{i}","fasta"):
                    seqs.append(SeqRecord(records.seq,records.id))
            align= MultipleSeqAlignment(seqs)
            name = input("name of file:")
            with open(name +".aln", "w") as handle:
                SeqIO.write(align, handle, "clustal")