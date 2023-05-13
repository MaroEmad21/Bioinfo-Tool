# making a program that makes gene 
# first of all we make a class that conatin the meaning of gene 

# the sequence object that contain all nucleotide
class Dna:
    def __init__(self, sequence):
        self.sequence = sequence
    
    def add_sequence():
        try:
            addSeq = input("what is your sequence: ")
            x =Dna(addSeq)
            print(Dna.get_sequence(x))
        except ValueError:
            print("YOU DIDN'T TYPE ANYTHING!!!!")


    def get_sequence(self):
        return self.sequence
    def check_sequence(self):
        pass





class Gene(Dna):
    def recog_promoter():
        pass
    def add_promoter():
        pass
    pass 
