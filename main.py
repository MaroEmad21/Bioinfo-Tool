# here is th file where the app works
# first we import all classes
# this file is the main file whic will run all other files till now
from gene import *
# AS I WILL USE THIS MAIN FUCTION TO BE THE MODERATOR OF ALL FUNCTIONS ANDD TO RUN ALL OPERTAIONS NEEDED
def main():
    answer = input("Which sequence you want to add? ")
    if answer.lower() == 'dna':
        Dna.add_sequence()
    else:
        print("that's not true")
        
            



main()            