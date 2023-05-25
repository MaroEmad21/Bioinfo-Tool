# here is th file where the app works

# first we import all classes
# this file is the main file whic will run all other files till now
#it will be for testing at first then after finishing wi
from gene import *
import random



randomStr=''.join([random.choice(Nucleotides)
                   for nuc in range(20)])
# AS I WILL USE THIS MAIN FUCTION TO BE THE MODERATOR OF ALL FUNCTIONS ANDD TO RUN ALL OPERTAIONS NEEDED
def main():
    print(Dna.add_sequence(randomStr))
    print(Dna.NucFreq(randomStr))
    print(Dna.transcription(randomStr))    
            



main()            