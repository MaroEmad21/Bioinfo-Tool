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
from Bio.Data.CodonTable import TranslationError
import glob
import random
from primer3 import calc_tm ,calc_hairpin
from Bio.Restriction import *
from Bio.Restriction.Restriction import RestrictionBatch
from pydna.gel import gel 
from pydna.ladders import GeneRuler_1kb
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from pydna.amplify import pcr,Anneal
from pydna.primer import Primer
from reportlab.lib import colors
from reportlab.lib.units import cm
from pydna.utils import rc
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
def get_primer_seq(sequence,id):
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






# a true primer designer 
def primer_designer(seq):
    for nuc in range(0,len(seq),1):
        if seq[nuc].upper() == "U":
            seq = seq.back_transcribe()
            amplicon =  primer_design(Dseqrecord(seq))
            return amplicon
        else:
            continue    
    amplicon =  primer_design(Dseqrecord(seq))
    test = open("primer.txt",'w')
    test.write(str(amplicon.figure()))
    test.write("\n \n")
    test.write(f"forward: {amplicon.forward_primer.seq}")
    test.write("\n")
    test.write(f"reverse: {amplicon.reverse_primer.seq}")
    test.write("\n \n")
    test.write(f"Program of PCR using Taq polymerase: {amplicon.program()}")
    test.close()
    return amplicon

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






# A new translate Function that produces 
# a protein with the index of the start codon
def new_translate(seq):
    n=0
    while True:
        if len(seq) // 3 == 0:
            pass
        else:
            seq = seq+"N"
            break
    for n in range(0,len(seq),1):
        if seq[n:n+3].translate()==TranslationError:
            n+=1
            continue
        elif seq[n:n+3].translate() == "M":
            protein= seq[n:].translate(table=2,to_stop=True)
            print(f"index of the start codon: {n , n+3}")
            return protein  
        else:
            continue


# enzyme map with mostly all commersial enzymes
def enzyme_map(sequence):
    # here it maps all the enzymes that cuts the sequence
    C_contain=RestrictionBatch([],['X','B','I'])
    Analog= Analysis(C_contain,sequence,False)
    # different types of cuts  
    # loop to allow different times 
    while True:
            ask =  input("""which type of analysis? 
full(1)
blunt(2)
overhang5(3)
overhang3(4)
defined(5)
with_sites(6)
without_site(7)
with_N_sites(8)
with_site_size(9)
stop(x)
""")
            if ask == "1": 
                answer = Analog.full()
            elif ask == "2":
                answer = Analog.blunt()
            elif ask == "3":
                answer = Analog.overhang5() 
            elif ask == "4":
                answer = Analog.overhang3() 
            elif ask == "5":
                answer = Analog.defined() 
            elif ask == "6":
                answer = Analog.with_sites() 
            elif ask == "7":
                answer = Analog.without_site() 
            elif ask == "8":
                N =int(input("number of sites: ")) 
                answer = Analog.with_N_sites(N) 
            elif ask == "9":
                N =int(input("site size: ")) 
                answer = Analog.with_site_size(N)
            elif ask.lower() == "x":
                break
            Analog.print_that(answer)
            Analog.print_as("map")
            # saves it in a file of text format
            test = open("map.txt",'w')
            test.write(str(Analog.format_output(answer)))
            test.close()
    #this will be made later 
    """ # testing making a draw
    gd_diagram = GenomeDiagram.Diagram(id)
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()
    for enzyme in C_contain:
        for site, name, color in [
        (enzyme.site,str(enzyme),colors.lightsteelblue)
        ]:
            index = 0
            while True:
                index = sequence.find(site, start=index)
                if index == -1:
                    break
                feature = SeqFeature(SimpleLocation(index, index + len(site)))
                gd_feature_set.add_feature(
                feature,
                color=color,
                name=name,
                label=True,
                label_size=10,
                

)
                index += len(site)
                gd_diagram.draw(
                            format="circular",
                            circular=True,
                            pagesize=(40 * cm, 40 * cm),
                            start=0,
                            end=len(sequence),
                            circle_core=0.5,

)               
                
                #feature = SeqFeature(SimpleLocation(index_of_M(sequence),index+len(sequence)),type="ORF").translate(sequence,start_offset=1,to_stop=True,)
                #gd_feature_set.add_feature(feature, sigil="ARROW", color="blue", arrowhead_length=0.25)
                #gd_diagram.write("test.pdf", "Pdf",dpi=300)
"""
                



# PCR maker
def make_pcr(sequence):
    # first check if primer anneal and how 
    amplicon=primer_designer(sequence)
    fp = Primer(amplicon.forward_primer.seq)
    rp = Primer(amplicon.reverse_primer.seq)
    anneal = Anneal([fp,rp],Dseqrecord(sequence))
    print(f"report: {anneal.report()}\n")
    #now we add sites of the restriction enzyme
    main_q= input("want to add enzyme?(y or n) ")
    if main_q.lower() == "y":
        sec_q = int(input("how many enzymes?(1 or 2) "))
        if sec_q == 1:
            name= input("write enzyme name ")
            enz = RestrictionBatch([name])
            for i in enz:
                site = i.site
            site = Seq(site)
            fp = site  + fp
            rp = site.reverse_complement() + rp  
            product= pcr(fp,rp,sequence)
            print(product.figure())
            return product
        elif sec_q == 2:
            name1=   input("write first enzyme name(forward enzyme) ")                  
            enz1= RestrictionBatch([name1])
            for a in enz1:
                site1 = Seq(a.site)
            name2=  input("write second enzyme name ")
            enz2= RestrictionBatch([name2])
            for b in  enz2:
                site2 = Seq(b.site)
            fp = site1 + fp
            rp = site2.reverse_complement() + rp
            product = pcr(fp,rp,sequence)
            print(product.figure())
            return product
        else:
            raise ValueError("Wrong value (1 or 2) only ")
    elif main_q.lower() == "n":
        product = pcr(fp,rp,sequence)
        print(product.figure())
        return product
    else:
        raise ValueError("please on only choose y or n")    