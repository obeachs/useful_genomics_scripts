from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
from Bio.Data import CodonTable
import argparse
import sys
import linecache

def run(args):
    l = Seq(args.input)
    print("Your sequence = " + "5- " + l + " -3")
    print("Reverse = " + "3- " +l[::-1] + " -5")
    print("Complement = " + "3- " + l.complement() + " -5")
    print("Transcribed = " + "5- " + l.transcribe() + " -3")
    if len(l)%3==0:
        print("Translated = " + l.translate())
    #This is a cool function that lets us look at the codon table
    table = CodonTable.unambiguous_dna_by_name["Standard"]
    #We can also only show the stop codons.
    #Seq objects are not mutable, we cannot change any nucleotides in them
    #unless we call them with the .tomutable() - which involves importing
    #yet another module MutableSeq. The main concern with using MutableSeq is that
    #they cannot be used as dictionaries - if they're changed,they're changed
    #This is why it may be wise to call a Seq variable initially and then mutate it
    #with MutableSeq. A MutableSeq object can be reverted to a Seq object fairly
    #easily with object.toseq()
    print(table.stop_codons)
    print(table)



def main():
	parser=argparse.ArgumentParser(description="Trim a reads to a specific length ")
	#Adding arguments is much simpler than expected, the 'dest' component of it
	#is what is called in the functions above, allows for as many arguments input
	#as you  want.
	parser.add_argument("-in",help="Sequence 5prime to 3prime " ,dest="input", type=str, required=True)
	parser.set_defaults(func=run)
	args=parser.parse_args()
	args.func(args)

if __name__=="__main__":
	main()
