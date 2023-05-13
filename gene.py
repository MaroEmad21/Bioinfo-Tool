# making a program that makes gene 
# first of all we make a class that conatin the meaning of gene 

# the sequence object that contain all nucleotides
class Dna:
    def __init__(self, sequence):
        self.sequence = sequence
    
    # a function to add sequence
    def add_sequence():
        try:
            addSeq = input("what is your sequence: ")
            seq =Dna(addSeq)
            return Dna.get_sequence(seq)
        except None:
            print("YOU DIDN'T TYPE ANYTHING!!!!")

    # function that returns the sequence for testing
    def get_sequence(self):
        seq = self.sequence
        #sequnce must be in upper case
        if seq:
            return seq.upper()
        else:
            return seq
    # a funcrion that check the sequence to see if it's correct or not (not random typing)
    def check_sequence(self):
        pass





class Gene(Dna):
    def recog_promoter():
        pass
    def add_promoter():
        pass
    pass 
