#!/usr/bin/python
'''
Created November 2021
@author: chasew
takes in a specifically formatted fasta file and text file with kmers to find chromosomal and sequence locations of all kmer instances in fasta.

Usage: python Locate_kmers.py --fasta FASTA.collapsed.reformat.fa --Search kmers.txt


'''

import os,sys,argparse

def parse_kmer_txt(txt):

	with open(txt,'r') as infile:
		lines = infile.readlines()

	to_search = []
	for line in lines:
		to_search += [line.rstrip()]

	return to_search

def getLoc(info,index,k): #["FBgn00300343-1","loc=X:19961689..19961845","length=157","parent=FBgn0031081","FBtr=FBtr0070000,FBtr0307554"]
	transcript = ":".join(info[4].split('=')[1].split(','))
	parent = info[3].split("=")[1]

	chrom,locinfo = info[1].split(":")
	chrom = chrom.split("=")[1]
	
	if "complement" in locinfo:
		if "))," in locinfo:
			strand = "+-"
			locinfo = locinfo[0:locinfo.find(")),")] + locinfo[locinfo.find(")),")+2:]
		else:	
			strand = "-"
	else:
		strand = "+"

	locinfo = locinfo.strip("complementjoin( )").split(',')
	
	region = []

	for chunk in locinfo:
		start,end = chunk.split('..')
		region += [i for i in range(int(start),int(end)+1)]

	region = ','.join(str(region[index:index+k]).strip('[]').split(', '))

	return [transcript,parent,chrom,strand,region]

def checkSearched(kmer,searched):
	check = False
	for seq in searched:
		if seq in kmer:
			check = True

	return check

def labelSite(k):
	nts = []
	for i in range(k):
		nts += ["nt{}".format(i+1)]
	nts += ["RNAlength","RNApos1"]
	return ",".join(nts)

def Read_kmers(fasta,searched,outname):
	
	k = len(searched[0])
	for term in searched[1:]:
		assert len(term) == k, 'All searched kmers must be equal lengths'

	with open(fasta,'r') as infile:
		lines = infile.readlines()

	with open("{}_KMERS.csv".format(outname),'w') as outfile:
		outfile.write("transcript,parent,kmer,chrom,strand,{},sequence\n".format(labelSite(k)))	

	for i in range(0,len(lines),2):
		info = lines[i].rstrip().split(";")
		info = info[:5] + [info[-2]]
		info = info[0].split(" ") + info[1:]
		info[0] = info[0].lstrip('>')
		del info[1]
		del info[2]
		for j in range(len(info)):
			info[j] = info[j].strip() # this has created a list of info for the sequence ["FBgn00300343-1","loc=X:19961689..19961845","length=157","parent=FBgn0031081","FBtr=FBtr0070000,FBtr0307554"]

		seq = lines[i+1].rstrip()
		length = int(info[2].split("=")[1])

		for pos in range(0,len(seq)-k+1):

			kmer = seq[pos:pos+k]

			if checkSearched(kmer,searched):

				transcript,parent,chrom,strand,region = getLoc(info,pos,k)
				outstring = ",".join([transcript,parent,kmer,chrom,strand]) + ",{},{},{},{}\n".format(region,length,pos+1,seq)

				with open("{}_KMERS.csv".format(outname),'a') as outfile:
					outfile.write(outstring)

if __name__ == '__main__': # allows another python script to import the functions

	parser = argparse.ArgumentParser(description="Locates all instances of specified kmers in a formatted fasta", usage='%(prog)s --fasta --Search',add_help=False) # a description of the function

	required = parser.add_argument_group('Required Input', 'These specifications are necessary to run.')

	required.add_argument("--fasta",help="Input fasta file including sequences to search", required=True)
	required.add_argument("--Search",help='Text file with kmers of equal length on each line: for example AGTTGA,AGTTAA,AGTGGA.', required=True)

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the script is generally run.')
	data_opt.add_argument('-O',"--OutPrefix",type=str,default="SEARCHED",help="Specify the prefix of the output _KMERS.csv. Default is SEARCHED")

	help_opt = parser.add_argument_group('Help')
	help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")

	args = parser.parse_args()

	Read_kmers(args.fasta,parse_kmer_txt(args.Search),args.OutPrefix)