#!/usr/bin/python
'''
Created November 2020
@author: chasew
Corrects a fasta file to specific format

Usage: python reformatdmel_region_fasta.py --fasta FASTA.fasta

'''

import os,sys,argparse

def TupInfo(Type,region,MD5,length,FBgn,relspec,seqline):
	tupinfo = tuple(Type + region + MD5 + length + FBgn + relspec + [seqline]) 
	# Makes ('type=CDS','loc=X:19963955..19964071','MD5=5c3c0b466b23e32d99dfff926a5e8c6b', ' length=2361',' parent=FBgn0031081',' release=r6.36';' species=Dmel',seqline)

	return tupinfo

def ParseUTR(info,seqline):
	split1 = info[0].find("type")
	Type = [info[0][split1:]]
	FBtr = info[0][:split1].strip("> ")
	region = [info[1]]
	MD5 = [info[2].strip()]
	length = [info[3]]
	FBgn = [info[4]]
	relspec = info[-3:-1]

	newinfo = TupInfo(Type,region,MD5,length,FBgn,relspec,seqline)

	return [newinfo,FBtr]

def ParseCDS(info,seqline):

	if "(mdg4)" not in info[0]: # excludes certain unique trans-spliced transcripts
		Type = [info[0][info[0].find("type"):]]
		region = [info[1]]
		MD5 = [info[4].strip()]
		length = [info[5]]
		FBgn,FBtr = info[6].split(',')
		FBgn = [FBgn]
		relspec = info[-3:-1]

		newinfo = TupInfo(Type,region,MD5,length,FBgn,relspec,seqline)

	else:
		newinfo = None
		FBtr = None

	return [newinfo,FBtr]

def ParsemRNA(info,seqline):
	if "(mdg4)" not in info[0]:
		Type = [info[0][info[0].find("type"):]]
		region = [info[1]]
		MD5 = [info[5].strip()]
		length = [info[6]]
		FBgn = [info[7]]
		FBtr = info[2][info[2].find("=")+1:]
		relspec = info[-3:-1]

		newinfo = TupInfo(Type,region,MD5,length,FBgn,relspec,seqline)

	else:
		newinfo = None
		FBtr = None

	return [newinfo,FBtr]


def GetUnique(fasta,Type="CDS"):

	with open(fasta,'r') as infile:
		lines = infile.readlines()

	AllSeqs = {}

	trs = []
	gns = []
		
	for i in range(0,len(lines),2):
		info = lines[i].rstrip().split(';')
		seqline = lines[i+1]
		
		if Type == "CDS":
			newinfo,FBtr = ParseCDS(info,seqline)

		elif Type == "UTR":
			newinfo,FBtr = ParseUTR(info,seqline)

		elif Type == "mRNA":
			newinfo,FBtr = ParsemRNA(info,seqline)
		
		else:
			newinfo,FBtr = [None,None]

		if newinfo is not None:

			if FBtr not in trs:
				trs += [FBtr]
			if newinfo[4] not in gns:
				gns += [newinfo[4]]

			if newinfo not in AllSeqs:
				AllSeqs[newinfo] = [FBtr]
			else:
				AllSeqs[newinfo] += [FBtr]

	return [AllSeqs,len(trs),len(gns)]

def CorrectFasta(fasta,Type="CDS"):
	
	AllSeqs,num_txs,num_gns = GetUnique(fasta,Type=Type)

	FBgns = {}

	with open(os.path.splitext(fasta)[0]+".reformat.fa",'w') as outfile:
		for seq in AllSeqs:
			seqline = seq[-1]
			newinfo = seq[:-1]
			FBgn = newinfo[4].split('=')[1]
			
			if FBgn in FBgns:
				FBgns[FBgn] += 1
			else:
				FBgns[FBgn] = 1

			count = FBgns[FBgn]

			newinfo = tuple([">{}-{} ".format(FBgn,count)+newinfo[0]]) + newinfo[1:]
			newline = ";".join(newinfo) + ";FBtr={};\n".format(','.join(AllSeqs[seq]))

			outfile.write(newline)
			outfile.write(seqline)

	return [num_txs,num_gns]

if __name__ == '__main__': # allows another python script to import the functions

	parser = argparse.ArgumentParser(description="Reformats a flybase transcript, CDS, or UTR file to be used with Count_Sites", usage='%(prog)s --fasta FASTA',add_help=False) # a description of the function

	required = parser.add_argument_group('Required Input', 'Specify the line-ended fasta to reformat')

	required.add_argument("--fasta",help="Input fasta file to reformat", required=True)

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the program runs.')
	data_opt.add_argument('-T',"--Type",default="CDS",help="specify type of fasta as mRNA, UTR, or CDS")

	help_opt = parser.add_argument_group('Help')
	help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")

	args = parser.parse_args()

	num_txs,num_gns = CorrectFasta(args.fasta,Type=args.Type)

	print("{} transcripts were reformatted including {} genes.".format(num_txs,num_gns))



