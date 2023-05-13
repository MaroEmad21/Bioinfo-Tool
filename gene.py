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
            Dna.get_sequence(seq)
            
        except ValueError:
            print("YOU DIDN'T TYPE ANYTHING!!!!")

    # function that returns the sequence for testing
    def get_sequence(self):
        seq = self.sequence
        #sequnce must be in upper case
        if seq:
            print (seq.upper())
        else:
            print("error")
        
        print( Dna.check_sequence(seq))
        
    # a funcrion that check the sequence to see if it's correct or not (not random typing)
    def check_sequence(self):
        seq= self.sequence
        seq_list=[]
        seq_list.append(seq)
        for nucleotide in seq_list:
            if nucleotide == "A" or nucleotide == "C" or nucleotide == "G" or nucleotide == "T":
                print("that's correct")
            else:
                print("something is wrong")
        





class Gene(Dna):
    def check_promoter():
        pass
    def add_promoter():
        pass
    pass 
