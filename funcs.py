"""
this file will contain all functions needed to be made for the tool that 
doesn't exist in libraries
"""

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import   SeqIO  ,Phylo
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import *
from Bio.Align import  MultipleSeqAlignment , PairwiseAligner
import os
import glob

# function that makes alignment and saves it in a clustal file
# NOTE!!! you should put all files to be aligned in one folder
def make_alignment(seq):
    Q1= int(input("how many sequences? "))
    if Q1 == 2:
        # pairwise alignment
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
            alignment=aligner.align(seq,sequence2)
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
    # Multiple sequence alignment (MSA)     
    # remember to put all needed files in one folder               
    elif Q1 >= 3:
            # asks if it's a file or a folder             
            Q5 = input("File or Folder (1 or 2) ")
            if Q5 == "1":
                seqs=[]
                file_path= input("file name: ")
                file_type= input("file type: ")
                for records in SeqIO.parse(f"{file_path}", f"{file_type}"):
                    seqs.append(SeqRecord(records.seq,records.id))
                align= MultipleSeqAlignment(seqs)
                name = input("name of file:")
                with open(name +".aln", "w") as handle:
                    SeqIO.write(align, handle, "clustal")
            elif Q5 == "2":
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
            else:    
                print("WRONG!!!!")




# makes phylo genetic tree
def make_phylo(seq,ids):
    info = input("do you want to include your sequence? (y or n) ")
    if info == "y":
        info2 = input("File or Folder? (1 or 2) ")
        if info2 == "1":
            seqs=[]
            file_path= input("file name: ")
            file_type= input("file type: ")
            for records in SeqIO.parse(f"{file_path}", f"{file_type}"):
                seqs.append(SeqRecord(records.seq,records.id))
            seqs.append(SeqRecord(seq,ids))
            align= MultipleSeqAlignment(seqs)
            # Calculate the distance matrix
            calculator = DistanceCalculator('identity')
            distMatrix = calculator.get_distance(align)
            # Create a DistanceTreeConstructor object
            constructor = DistanceTreeConstructor()
            # Construct the phlyogenetic tree using UPGMA algorithm
            UPGMATree = constructor.upgma(distMatrix)
            # Construct the phlyogenetic tree using NJ algorithm
            NJTree = constructor.nj(distMatrix)
            # Draw the phlyogenetic tree
            Phylo.draw(UPGMATree)
            # Draw the phlyogenetic tree using terminal
            Phylo.draw_ascii(NJTree)
        elif info2 == "2":
            """
                it takes files in folder and makes a phylo tree 
            """
            path=input("path of files: ")
            files = []
            seqs=[]
            for file in glob.glob(f"{path}"+"\*.fasta"):
                files.append(file)
            for i in files:
                for records in  SeqIO.parse(f"{i}","fasta"):
                    seqs.append(SeqRecord(records.seq,records.id))
            seqs.append(SeqRecord(seq,ids))
            align= MultipleSeqAlignment(seqs)
            # Calculate the distance matrix
            calculator = DistanceCalculator('identity')
            distMatrix = calculator.get_distance(align)
            # Create a DistanceTreeConstructor object
            constructor = DistanceTreeConstructor()
            # Construct the phlyogenetic tree using UPGMA algorithm
            UPGMATree = constructor.upgma(distMatrix)
            # Construct the phlyogenetic tree using NJ algorithm
            NJTree = constructor.nj(distMatrix)
            # Draw the phlyogenetic tree
            Phylo.draw(UPGMATree)
            # Draw the phlyogenetic tree using terminal
            Phylo.draw_ascii(NJTree)
            raise ValueError("try again")
    elif info == "n":
        info2 = input("File or Folder? (1 or 2) ")
        if info2 == "1":
            seqs=[]
            file_path= input("file name: ")
            file_type= input("file type: ")
            for records in SeqIO.parse(f"{file_path}", f"{file_type}"):
                seqs.append(SeqRecord(records.seq,records.id))
            align= MultipleSeqAlignment(seqs)
            # Calculate the distance matrix
            calculator = DistanceCalculator('identity')
            distMatrix = calculator.get_distance(align)
            # Create a DistanceTreeConstructor object
            constructor = DistanceTreeConstructor()
            # Construct the phlyogenetic tree using UPGMA algorithm
            UPGMATree = constructor.upgma(distMatrix)
            # Construct the phlyogenetic tree using NJ algorithm
            NJTree = constructor.nj(distMatrix)
            # Draw the phlyogenetic tree
            Phylo.draw(UPGMATree)
            # Draw the phlyogenetic tree using terminal
            Phylo.draw_ascii(NJTree)
        elif info2 == "2":
            """
                it takes files in folder and makes a phylo tree 
            """
            path=input("path of files: ")
            files = []
            seqs=[]
            for file in glob.glob(f"{path}"+"\*.fasta"):
                files.append(file)
            for i in files:
                for records in  SeqIO.parse(f"{i}","fasta"):
                    seqs.append(SeqRecord(records.seq,records.id))
            align= MultipleSeqAlignment(seqs)
            # Calculate the distance matrix
            calculator = DistanceCalculator('identity')
            distMatrix = calculator.get_distance(align)
            # Create a DistanceTreeConstructor object
            constructor = DistanceTreeConstructor()
            # Construct the phlyogenetic tree using UPGMA algorithm
            UPGMATree = constructor.upgma(distMatrix)
            # Construct the phlyogenetic tree using NJ algorithm
            NJTree = constructor.nj(distMatrix)
            # Draw the phlyogenetic tree
            Phylo.draw(UPGMATree)
            # Draw the phlyogenetic tree using terminal
            Phylo.draw_ascii(NJTree)
            raise ValueError("try again!!")
        else:
            print("wrong choice")
    else:
        print("choose well!!!!!")



