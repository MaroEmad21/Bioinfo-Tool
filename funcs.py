"""
this file will contain all functions needed to be made for the tool that 
doesn't exist in libraries
"""

from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import   SeqIO  ,Phylo
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import *
from Bio.Align import  MultipleSeqAlignment , PairwiseAligner
import os
import glob
import random 
from primer3 import calc_tm ,calc_hairpin
from Bio.Restriction import AllEnzymes
from Bio.Restriction.Restriction import RestrictionBatch
from pydna.gel import gel 
from pydna.ladders import GeneRuler_1kb
from pydna.dseqrecord import Dseqrecord
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
"""
it takes the sequence inputed in case user wanted to compare it with other files 
        after making all needed analysis or expression
"""
def make_phylo(seq,ids):
    info = input("do you want to include your sequence? (y or n) ")
    if info.lower() == "y":
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
    elif info.lower() == "n":
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

""" 
    Primer desgining that depends on the Restriction Enzyme 
    But till now it will generate random sequence with certain conditions
"""
#primer generator
"""
    Generates as random sequence that could be added to the sequence
    as an extended part to clone th whole sequence
    with properties of a good primer and not forming a hairpin structure
""" 



def get_primer(sequence,id):
    Q1 = input("Primer Within sequence or out of sequences(1 or 2)? ")
    if Q1 == "1":
        print(Q1)
    elif Q1 == "2":
        size= int(input("size [18 : 24]"))
        if size > 24 or size < 18:
           print("this will be problem")
        else:
            Nucleotides= ["A","C","G","T"]
            Forward=''.join([random.choice(Nucleotides)
                                for nuc in range(size)])
            reverse=''.join([random.choice(Nucleotides)
                                for nuc in range(size)]) 
        while True:
            if round(gc_fraction(Forward)*100) in range(40,60) and round(calc_tm(Forward)) in range(50,60) and round(gc_fraction(reverse)*100) in range(40,60) and round(calc_tm(reverse)) in range(50,60) and str(calc_hairpin(Forward)) == "ThermoResult(structure_found=False, tm=0.00, dg=0.00, dh=0.00, ds=0.00)" and  str(calc_hairpin(reverse)) == "ThermoResult(structure_found=False, tm=0.00, dg=0.00, dh=0.00, ds=0.00)": 
                print(Forward) 
                print(round(gc_fraction(Forward)*100))
                print(round(calc_tm(Forward)))
                print(calc_hairpin(Forward))
                print(reverse) 
                print(round(gc_fraction(reverse)*100))
                print(round(calc_tm(reverse)))
                print(calc_hairpin(reverse))
                align= PairwiseAligner()
                align.mode == "local"
                print(align.align(Forward,reverse)[0])
                #mainQ=input("do you want")
                add_primer(sequence,Forward,reverse,id)
                break
            else:
                Forward=''.join([random.choice(Nucleotides)
                                for nuc in range(size)])
                reverse=''.join([random.choice(Nucleotides)
                                for nuc in range(size)])
                continue
    else :
        raise ValueError("WRONG VALUE")        
    




#this function adds primer to the sequnece 
def add_primer(seq,forward,reverse,id):
    print(f"{colored(forward)}{seq}{colored(reverse)}")
    newsequence= Seq(f"{forward}{seq}{reverse}")
    newsequence= SeqRecord(newsequence,id=id,description="sequence with primer")
    name = input("name of file:")
    with open(name +".fasta", "w") as handle:
        SeqIO.write(newsequence, handle, "fasta")




# to colorize sequences
def colored(seq):
    bcolors = {
        'A': '\033[92m',
        'C': '\033[92m',
        'G': '\033[92m',
        'T': '\033[92m',
        'U': '\033[92m',
        'reset': '\033[0;0m'
    }

    tmpStr = ""

    for nuc in seq:
        if nuc in bcolors:
            tmpStr += bcolors[nuc] + nuc
        else:
            tmpStr += bcolors['reset'] + nuc

    return tmpStr + '\033[0;0m'
            
            



Nucleotides= ["A","C","G","T"]

Forward=''.join([random.choice(Nucleotides)
                                for nuc in range(20)])

