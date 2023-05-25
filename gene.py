# making a program that makes gene 
import collections

# first of all we make a class that conatin the meaning of gene 
# the sequence object that contain all nucleotides


Nucleotides= ["A","C","G","T"]


class Dna:
    def __init__(self, sequence):
        self.sequence = sequence
    
    # a function to add sequence
    def add_sequence(sequence):
            try:
                seq =Dna(sequence)
                Dna.get_sequence(seq)  
            except ValueError:
                print("YOU DIDN'T TYPE ANYTHING!!!!")
    # function that returns the sequence for testing
    def get_sequence(self):
        seq = self.sequence
        #sequnce must be in upper case
        if seq:
            seq=seq.upper()
            # to check sequence if correct will continue else program will break
            Dna.check_sequence(seq)
        else:
            print("error")
        
    # a funcrion that check the sequence to see if it's correct or not (not random typing)
    def check_sequence(self):
        for nucleotide in self:
            if nucleotide in Nucleotides:
                pass
            else:
                print(f"that's not correct the nucleotide {nucleotide} is wrong")
                quit()
        print(f"sequence is [{self}] and the length is {len(self)} bp")


    # a function to count the frequency of each nucleotide
    def NucFreq(sequence):
        return dict(collections.Counter(sequence))
    


    #transcription function
    def transcription(seq):

         return seq.replace("T","U")     




#gene is part of dna so I will make a subclass that gets a nucleotide sequence and then adds promoter and gives transcripts
'''class Gene(Dna):
    def check_promoter():
        pass
    def add_promoter():
        pass
    pass 
'''
