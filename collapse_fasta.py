#!/usr/bin/python
'''
Created October 2020
@author: chasew
Collapses a fasta file to be two-line: <info and sequence

Usage: python collapse_fasta.py --fasta FASTA.fasta

'''

import os,sys,argparse

def CollapseFasta(fasta):

	with open(fasta,'r') as infile:
		lines = infile.readlines()

	with open(os.path.splitext(fasta)[0]+".collapsed.fa",'w') as outfile:
		faString = ""
		
		outfile.write(lines[0]) # write first line

		for line in lines[1:]:
			if line[0] == ">":
				outfile.write(faString+"\n") # add sequence
				faString = ""
				outfile.write(line) # add next fasta name
			else:
				faString += line.rstrip() # concatenate sequence

		outfile.write(faString+"\n") # write last line


if __name__ == '__main__': # allows another python script to import the functions

	parser = argparse.ArgumentParser(description="Corrects a fasta file to be two-line: <info and sequence", usage='%(prog)s --fasta FASTA',add_help=False) # a description of the function

	required = parser.add_argument_group('Required Input', 'Specify the line-ended fasta to correct')

	required.add_argument("--fasta",help="Input fasta file to correct", required=True)

	help_opt = parser.add_argument_group('Help')
	help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")

	args = parser.parse_args()

	CollapseFasta(args.fasta)



