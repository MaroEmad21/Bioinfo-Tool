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
from Bio.Restriction import CommOnly,AllEnzymes
import glob
import random
from pydna.readers import read
from primer3 import calc_tm ,calc_hairpin
from Bio.Restriction import Analysis
from Bio.Restriction.Restriction import RestrictionBatch
from pydna.gel import gel
from pydna.ladders import GeneRuler_1kb
from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from pydna.amplify import pcr,Anneal
from pydna.common_sub_strings import terminal_overlap
from pydna.primer import Primer
from pydna.assembly import Assembly
from pydna.design import assembly_fragments
from pydna.amplicon import Amplicon
from reportlab.lib import colors
from reportlab.lib.units import cm
from pydna.utils import rc
from pydna.parsers import parse

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
            protein= seq[n:].translate(table=2)
            print(f"index of the start codon: {n , n+3}")
            return protein
        else:
            continue


# enzyme map with mostly all commersial enzymes
def enzyme_map(sequence):
    # here it maps all the enzymes that cuts the sequence
    C_contain=RestrictionBatch([],['X','B','I'])
    Analog= Analysis(AllEnzymes,sequence,False)
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
            #type= input("fromat type: ")
            Analog.print_as("list")
            Analog.print_that(answer)
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

def Orf(sequence):
    proteins = []
    i=1
    for orf in Dseqrecord(sequence).orfs():
        proteins.append(f"protein no. {i}: {orf.translate().seq}")
        i+=1
    return proteins
"""  And now it's to for making the cloning Function.
Since Most plasmids are
"""
def cloning(sequence):
    answer = int(input("Sub cloning or homolgy(1 or 2): "))
    if answer == 1:
        """Restriction cloning by restriction digestion and ligation"""
        #vector_path = input("file name: ")
        vector =  read("puc.gb")
        cutters = []
        # once cutters in vector
        for cutter in vector.n_cutters(n=1):
            cutters.append(cutter)
        print(f" one cutters with sticky ends: {cutters}")
        #user must choose an enzyme to linearize the vector
        main_Q = input("write name of the enzyme to linearize the vector: ")
        enzyme = RestrictionBatch([main_Q])
        for a in enzyme:
            linear_vector  = vector.linearize(a)
        # makes a map of enzymes
        enzyme_map(sequence)
        # making pcr
        prod=make_pcr(sequence)
        sec_q= input("name of the enzyme that cuts the insert: ")
        enzyme2 = RestrictionBatch([sec_q])
        for a in enzyme2:
            stuffer1, insert, stuffer2 = prod.cut(a)
        # then we assemble all together
        # recombined product is saved in txt file
        try:
            Recomb = (linear_vector+insert).looped()
            Recomb.write("test3.txt")
        # to retry if errors occured  if ends are not compatible
        except TypeError or ValueError:
            print("try again STICKY ENDS are not compaible!!!!!")
            cloning(sequence)
    # will be updated later
    elif answer == 2:
        """Cloning by homologous recombination"""
        #vector_path = input("file name: ")
        vector =  read("sequence.gb")
        # enzyme map for the vector
        enzyme_map(vector.seq)
        #user must choose an enzyme to cut the vector
        main_Q = input("write name of the enzyme to remove  the part from vector: ")
        enzyme= RestrictionBatch([main_Q])
        for a in enzyme:
            # remove the neede part to separate
            rOF_vector, homologous   = vector.cut(a)
        #search for the zero cutter in sequence
        enzyme_map(sequence)
        # add the site of
        sec_q= input("name of the enzyme that will be added: ")
        enzyme2 = RestrictionBatch([sec_q])
        for a in enzyme2:
            site = Dseqrecord(a.site)
        """ site is added to ensure assembly and to be used in gel"""
        seq_amp = primer_design(Dseqrecord(sequence))
        #assemble every thing together
        fragment_list=assembly_fragments((rOF_vector,site,seq_amp,rOF_vector))
        fragment_list=fragment_list[:-1]
        try:
            asm= Assembly(fragment_list)
            Recomb =asm.assemble_circular()[0]
            print(Recomb)
            Recomb.write("test4.txt")
        # to retry if errors occured  if ends are not compatible
        except TypeError or ValueError:
            print("try again !!!!!")
            cloning(sequence)
        # NOTE there will be further updates

# assembly maker
def make_assembly(sequence):
    """ Gibson assembly"""
    #vector_path = input("file name of vector: ")
    vector =  read("vector.gb")
    if (vector.linear) == True:
        vector = vector.looped()
    else:
        pass
    cutters = []
    # once cutters in vector
    for cutter in vector.once_cutters():
        cutters.append(cutter)
    print(f" once cutters: {cutters}")
    #user must choose an enzyme to linearize the vector
    main_Q = input("write name of the enzyme to linearize the vector: ")
    enzyme = RestrictionBatch([main_Q])
    for a in enzyme:
        linear_vector  = vector.linearize(a)
    #other_sec_path = input("file name of other sequence: ")
    other_seq= read("GFP.gb")
    print(other_seq.list_features())
    # n= int(input("no. of the feature to be extracted and replaced"))
    feature = other_seq.extract_feature(5)
    # to debug if the right feature
    print(feature.isorf())
    #make primer for the sequence
    seq_amplicon = primer_design(sequence)
    #show the figure
    print(seq_amplicon.figure())
    feature_amplicon = primer_design(feature)
    print(feature_amplicon.figure())

    #add all sequences together
    fragments = (linear_vector,seq_amplicon,feature_amplicon,linear_vector)
    # check for zero cutters in all fagments
    no_cutters=[]
    C_contain=RestrictionBatch([],['X','B','I'])
    for enzyme in C_contain:
        if not any( x.cut(enzyme) for x in fragments ):
            no_cutters.append(enzyme)
        else:
            continue
    print(no_cutters)
    #add this site to be used later in gel
    sec_q= input("name of the enzyme that will be added: ")
    enzyme2 = RestrictionBatch([sec_q])
    for a in enzyme2:
        site = Dseqrecord(a.site)

    fragments = (linear_vector,site,seq_amplicon,feature_amplicon,linear_vector)
    # assemble them all
    fragment_list= assembly_fragments(fragments)
    fragment_list=fragment_list[:-1]
    # assembly
    asm=Assembly(fragment_list)
    assembeled= asm.assemble_circular()
    print(assembeled[0].figure())
    assembeled[0].write("test5.txt")


Nucleotides= ["A","C","G","T"]

linkers=''.join([random.choice(Nucleotides)
                                for nuc in range(4)])







# this class wil contain two algorithms for making Golden Gate assembly
# The class will be updated soon if any issues occur 
# according on feedback we will test new algorithms
class GoldenGateAssembly:
    def __init__(self, vector_file, insert_files):
        self.vector_file = vector_file
        self.insert_files = insert_files
    def read_files(self):
        # Load the vector sequence from file
        vector_seq = read(self.vector_file)

        # Load the insert sequences from files
        insert_seqs = []
        for insert_file in self.insert_files:
            insert_seq = read(insert_file)
            insert_seqs.append(insert_seq)
        return insert_seqs,vector_seq

    def splitIntoFrags(self,enzyme):
    # Calculate and output the length of the total input sequence in both bp and kb
        self.insert_seqs,self.vector_seq  = self.read_files()

        seq_length = len(self.vector_seq)
        print('The length of the input sequence is %(length)i bp, or %(length_kb).3f kb' % {'length': seq_length, 'length_kb': seq_length/1000})
        # Let the user choose the number of fragments and set it to target_num
        # Output the sequences of split fragments according to user input
        enzyme = RestrictionBatch([enzyme])
        frag_list = []
        for i in (Dseqrecord(self.vector_seq).cut(enzyme)):
            save = i
            frag_list.append(Dseqrecord(save.seq))
        return  enzyme,frag_list, self.insert_seqs

    def makeReverseComplement(self,sequence):
        comp_Base_Pairs = {
        "A" : "T",
        "C" : "G",
        "T" : "A",
        "G" : "C"
        }
        complement = ''
        length = len(sequence)

        for i in range(length):
            complement = complement + comp_Base_Pairs[sequence[length - i - 1]]
        return complement

    #the function ADDOVERHANGS takes in the parameter frag_list, which is obtained by parsing fasta
    #and the parameter cut_site, which is returned from the ENZYMECUTSITESELECTION function
    #and outputs the frag_list containing the fragments with the correct overhangs
    def addOverhangs(self,frag_list,cut_site):
        result_list = []
        forward_cute_site = cut_site
        reverse_cut_site = self.makeReverseComplement(cut_site)
        random_codon = random.choice(['A','C','T','G'])
        for i in range(len(frag_list)):
            if i == len(frag_list)-1:
                start = frag_list[0][0:4]
            else:
                start = frag_list[i+1][0:4]
            new_frag = forward_cute_site + random_codon + frag_list[i] + start + random_codon + reverse_cut_site
            result_list.append(new_frag)
        return result_list


    def main(self):
        test,vector_seq= self.read_files()
        # first we make enzyme map for the seq
        enzyme_map(vector_seq.seq)
        #enzyme =  input("write the name of the enzyme: to cut the vector: ")
        enzyme,frag_list, insert_list = self.splitIntoFrags("BamHI")
        # to check for number of cuts
        for i in frag_list:
            print(i)
        # to save enzyme name
        for i in enzyme:
            enzyme = i
        # if more than 2 cuts remove the middle part    
        if len(frag_list) == 2:
            pass
        else:
            frag_list.pop(1)
        # makes random seq
        random_codon = ''.join([random.choice(Nucleotides)
                                for nuc in range(4)])

        #the loop for the adding inserts seq after cuts
        n = 1
        for insert in insert_list:
            insert = primer_design(Dseqrecord(insert))
            insert= enzyme.site + random_codon+ insert.forward_primer + insert.template + insert.reverse_primer + random_codon + enzyme.site 
            none,insert,none = Dseqrecord(insert).cut(enzyme)
            frag_list.insert(n,insert)
            n += 1
        # assembly for all of them    
        asm = Assembly(frag_list)
        Recomb= asm.assemble_circular()
        # to show the meaning of the empty brackets
        if Recomb == []:
            print("No sign for any assembly")
        else:
            print(Recomb)
            pass    
        # the other algorthim
        def golden_gate(self,sequence):

            """golden gate assembly"""
            #vector_path = input("file name of vector: ")
            vector =  read("vector.gb")
            if (vector.linear) == True:
                vector = vector.looped()
            else:
                pass
            enzyme_map(vector.seq)
            #user must choose an enzyme to linearize the vector
            main_Q = input("write name of the enzyme to linearize the vector: ")
            enzyme = RestrictionBatch([main_Q])
            for a in enzyme:
                enz = a
                # remove the neede part to separate
                vector1, vector2   = vector.cut(a)
            if vector1 >vector2:
                main_vector = vector1
            else:
                main_vector = vector2
            print(main_vector)
            test =  main_vector[:4]
            print(test)
            #other_sec_path = input("file name of other sequence: ")
            other_seq= read("sequence.fasta")
            #make primer for the sequence
            sequence = Dseqrecord(sequence)
            seq_amplicon = primer_design(sequence)
            #show the figure
            print(seq_amplicon.figure())
            new_amplicon= primer_design(other_seq)
            for i in (seq_amplicon,new_amplicon):
                i.forward_primer = enz.site + "A" + complement("GCCA") + i.forward_primer
                i.reverse_primer = enz.site + "A"+ self.makeReverseComplement("CTGT")  + i.reverse_primer
                i=primer_design(i,forward_primer=i.forward_primer,reverse_primer=i.reverse_primer)


            seq_amplicon = pcr(seq_amplicon.forward_primer,seq_amplicon.reverse_primer,seq_amplicon.template)
            new_amplicon = pcr(new_amplicon.forward_primer,new_amplicon.reverse_primer,new_amplicon.template)

            # cut both with the enzyme
            none,seq_amplicon,seq_amplicon,none=seq_amplicon.cut(enz)
            none,new_amplicon,none=new_amplicon.cut(enz)
            print(seq_amplicon)
            print(new_amplicon)
            try:
                #add them all together
                fragment_list = (main_vector,seq_amplicon)
                asm= Assembly(fragment_list,algorithm=terminal_overlap,limit=14)
                Recomb =asm.assemble_circular()
                print(Recomb)
                #Recomb.write("test4.txt")
            # to retry if errors occured  if ends are not compatible
            except TypeError or ValueError:
                print("try again !!!!!")
                golden_gate(sequence)
            # NOTE there will be further updates
