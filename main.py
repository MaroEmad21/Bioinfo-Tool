"""
    importing  libraries needed to do the tool
    this is a simple software that could help you to get better info  and data 
    that could be used In Bioinformatics
    this is Version 1.0  
"""
from Bio import SeqIO 
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import *
# to use functions made that not exist in Biopython
from funcs import *   



# First Question
answer = input( 
               """which sequece you will type?
(Dna or 1)
(Rna or 2)  
(enzyme map or 3) 
(Pcr or 4 )  
(Cloning  or 5) 
(assembly or 6) 
(Gel electrophoresis or 7) """)
print("KINDLY KNOW THAT THIS IS A BETA VERSION AND MORE UPDATES WILL COME")
if answer.upper() in ["DNA","1"]:
    file_path= input("file name: ")
    file_type= input("file type: ")
    for record in SeqIO.parse(f"{file_path}", f"{file_type}"):
        seq = record.seq
        ids = record.id
        print(f"id = {record.description} \n your sequence is: {seq}")

    # this loops to get more info from the same file
    while True:
        next_step = input(
    """what's next: 
    1)transcribe   
    2)translate  
    3)complement seq(comp) 
    4)reverse_complement(rev-comp) 
    5)reverse_complement_rna(5) 
    6)GC content(%) 
    7)Molecular weight(mw)
    8) Open Reading frames
    9)Alignment
    10)phylognetic tree
    11) Primer_designer
    __stop or x:  """)
        # Back transcription
        if next_step.lower() in  ['transcribe' ,'1']:
            print(f"m-RNA is: {seq.transcribe()}")
            records= SeqRecord(seq.transcribe(),id='',name='back_transcribe',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        # translate        
        elif next_step.lower() in ['translate','2']:
                print(f"Amino acid seq: {new_translate(seq)}") 
                records= SeqRecord(new_translate(seq),id='',name='translate',description=f'{record.description}')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")
        # complementary sequence           
        elif   next_step.lower() in ['comp','3'] :
                print(f"complement seq: {seq.complement()}")
                records= SeqRecord(seq.complement(),id='',name='complementary sequence',description=f'{record.description}')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")  
        #reverse complementary            
        elif next_step.lower() in ['rev-comp','4']:
            print(f"{seq.reverse_complement()}")
            records= SeqRecord(seq.reverse_complement(),id='',name='reverse complementary sequence',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta") 
        # Rna complementary sequence
        elif next_step.lower() in ['comp_rna','5']:
            print(f"complment rna:{seq.complement_rna()}")   
            records= SeqRecord(seq.complement_rna(),id='',name='reverse complementary_rna sequence',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        #GC content        
        elif next_step.lower() in  ['gc','6']:
            print(f"GC %: {round(gc_fraction(seq) * 100) }%")
        # molecular weight    
        elif next_step.lower() in ['mw', '7']:
            print(f'Moleceular weight: {molecular_weight(seq)}')
        # ORFs    
        elif next_step.lower() in ['sixframe', '8']:
            print(f'{Orf(seq) }')     
            name = input("name of file:")
            test = open(name+".txt",'w')
            for protein in Orf(seq):
                test.write(f"{protein}\n")
            test.close()
        # making alignment    
        elif next_step.lower() in ['al', '9']:
            make_alignment(seq)
        # makes and reads phylogenetic tree    
        elif next_step.lower() in ['phylo', '10']:
            make_phylo(seq,ids)
        elif next_step.lower() in ['pr', '11']:
            primer_designer(seq)     
        # to break the loop   
        elif next_step.lower() in  ['stop' ,'x']:
            break

    """
        this makes nearly same things as in Dna 
        but for Rna Except back transcribe 
    """

elif answer.upper() in ["rna","2"]:
    file_path= input("file name: ")
    file_type= input("file type: ")
    for record in SeqIO.parse(f"{file_path}", f"{file_type}"):
        seq = record.seq
        ids = record.id
        print(f"id = {record.description} \n your sequence is: {seq}")
    # this loops to get more info from the same file
    while True:
        next_step = input(
    """what's next: 
    1)Back transcribe   
    2)translate  
    3)complement seq(comp) 
    4)reverse_complement(rev-comp) 
    5)reverse_complement_rna(5) 
    6)GC content(%) 
    7)Molecular weight(mw)
    8) Open Reading frames
    9)Alignment
    10)phylognetic tree
    11) Primer_designer
    __stop or x:  """)
        # Back transcription
        if next_step.lower() in  ['transcribe' ,'1']:
            print(f"m-RNA is: {seq.back_transcribe()}")
            records= SeqRecord(seq.back_transcribe(),id='',name='back_transcribe',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        # translate        
        elif next_step.lower() in ['translate','2']:
                print(f"Amino acid seq: {new_translate(seq)}") 
                records= SeqRecord(new_translate(seq),id='',name='translate',description=f'{record.description}')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")
        # complementary sequence           
        elif   next_step.lower() in ['comp','3'] :
                print(f"complement seq: {seq.complement()}")
                records= SeqRecord(seq.complement(),id='',name='complementary sequence',description=f'{record.description}')
                name = input("name of file:")
                with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta")  
        #reverse complementary            
        elif next_step.lower() in ['rev-comp','4']:
            print(f"{seq.reverse_complement()}")
            records= SeqRecord(seq.reverse_complement(),id='',name='reverse complementary sequence',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                   SeqIO.write(records, handle, "fasta") 
        # Rna complementary sequence
        elif next_step.lower() in ['comp_rna','5']:
            print(f"complment rna:{seq.complement_rna()}")   
            records= SeqRecord(seq.complement_rna(),id='',name='reverse complementary_rna sequence',description=f'{record.description}')
            name = input("name of file:")
            with open(name +".fasta", "w") as handle:
                SeqIO.write(records, handle, "fasta")
        #GC content        
        elif next_step.lower() in  ['gc','6']:
            print(f"GC %: {round(gc_fraction(seq) * 100) }")
        # molecular weight    
        elif next_step.lower() in ['mw', '7']:
            print(f'Moleceular weight: {molecular_weight(seq)}')
        # ORFs    
        elif next_step.lower() in ['sixframe', '8']:
            print(f'{Orf(seq) }')     
            name = input("name of file:")
            test = open(name+".txt",'w')
            test.write(str(Orf(seq) ))
            test.close()
        # making alignment    
        elif next_step.lower() in ['al', '9']:
            make_alignment(seq)
        # makes and reads phylogenetic tree    
        elif next_step.lower() in ['phylo', '10']:
            make_phylo(seq,ids)
        # produce Primer    
        elif next_step.lower() in ['pr', '11']:
            primer_designer(seq)
        # to break the loop   
        elif next_step.lower() in  ['stop' ,'x']:
            break



# this to ask if making enzyme map        
elif answer.upper() in ["enz","3"]:
    file_path= input("file name: ")
    file_type= input("file type: ")
    for record in SeqIO.parse(f"{file_path}", f"{file_type}"):
        seq = record.seq
        ids = record.id
        print(f"your sequence is: {seq}")
    enzyme_map(seq)


 
# Pcr Product
elif answer.upper() in ["pcr","4"]:
    file_path= input("file name: ")
    file_type= input("file type: ")
    for record in SeqIO.parse(f"{file_path}", f"{file_type}"):
        seq = record.seq
        ids = record.id
        print(f"your sequence is: {seq}")
    make_pcr(seq)


#Cloning maker
elif answer.upper() in ["cloning","5"]: 
    file_path= input("file name: ")
    seq = read(f"{file_path}")
    cloning(seq)
    answer = input("Do you want to make Gel elctrophoresis(y or n):")
    if answer.lower() == "y":
        names = int(input("Number of enzymes used for gel"))
        enymes= []
        i = 0
        for i in range(i,names):
            enzyme = input("ENTER ENZYME NAME: ")
            enymes.append(enzyme)
        test = Gel(seq,enzymes=enymes)
        test.make_gel()
    elif answer.lower() == "n":
        pass    

# Here 2 types or 2 algorithms of assembly
elif answer.upper() in ["assembly","6"]: 
    file_path= input("file name: ")
    seq = Dseqrecord(read(f"{file_path}").seq)
    answer  = int(input("Gibson assembly or Golden Gate(1 or 2): "))
    if answer == 1:    
        make_assembly(seq)
        answer = input("Do you want to make Gel elctrophoresis(y or n):")
    if answer.lower() == "y":
        names = int(input("Number of enzymes used for gel"))
        enymes= []
        i = 0
        for i in range(i,names):
            enzyme = input("ENTER ENZYME NAME: ")
            enymes.append(enzyme)
        test = Gel(seq,enzymes=enymes)
        test.make_gel()
    elif answer.lower() == "n":
        pass    
    elif answer == 2:    
        print("""please note that there 2 algorithms used so there will be 2 products 
              let me know which one works better!! (Beta program) """)
        GoldenGateAssembly.golden_gate(seq)
        # Example usage
        vector_path = input("ENTER FILE PATH:")
        vector_file = f"{vector_path}"
        n = int(input("ENTER NUMBER OF FILES WILL BE USED"))
        insert_files = []
        for i in n:
            file = input("ENTER FILE NAME: ")
            insert_files.append(file)
        assembly = GoldenGateAssembly(vector_file, insert_files)
        recombined_vector = assembly.main()
        print(recombined_vector)
        answer = input("Do you want to make Gel elctrophoresis(y or n):")
    if answer.lower() == "y":
        names = int(input("Number of enzymes used for gel"))
        enymes= []
        i = 0
        for i in range(i,names):
            enzyme = input("ENTER ENZYME NAME: ")
            enymes.append(enzyme)
        test = Gel(seq,enzymes=enymes)
        test.make_gel()
    elif answer.lower() == "n":
        pass    
    else:
        raise ValueError("write the number 1 or 2!!!")
elif answer.upper() in ["Gel","7"]: 
    file_path= input("file name: ")
    seq = read(f"{file_path}")
    enzyme_map(seq.seq)
    names = int(input("Number of enzymes used for gel"))
    enymes= []
    i = 0
    for i in range(i,names):
        enzyme = input("ENTER ENZYME NAME: ")
        enymes.append(enzyme)
    test = Gel(seq,enzymes=enymes)
    test.make_gel()
else:    
    raise ValueError("Please Pick from the Choices!!")        



