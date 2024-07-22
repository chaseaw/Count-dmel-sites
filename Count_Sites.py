#!/usr/bin/python
'''
Created November 2021
@author: chasew
Takes in Kmer-location csvs and outputs per-transcript and per-gene counts.

Usage: python Count_Sites.py --csv KLOC.csv


'''

import os,sys,argparse

def create_dicts(csvs): # needs a list of csv filenames

	TransDicts = {}

	for csv in csvs:

		with open(csv,'r') as infile:
			lines = infile.readlines()

		for line in lines[1:]:
			info = line.strip().split(',') # 0,1,-3, and -2 are transcript, gene, length, and position info
			transcripts = info[0].split(':')

			for transcript in transcripts:
				KEY = (transcript,info[1],int(info[-3]))
				if KEY in TransDicts:
					TransDicts[KEY] += [int(info[-2])]
				else:
					TransDicts[KEY] = [int(info[-2])]

	# Now have dictionaries with keys = (transcript,gene,length) -> [site1 start, site2 start,...]

	for transcript in TransDicts:
		TransDicts[transcript].sort()

	return TransDicts

def GetStats(transdict): 

	StatDicts = {}

	
	for transcript in transdict:

		length = transcript[2]
		sites = transdict[transcript]
		num_sites = len(sites)
		
		if num_sites:
			StatDicts[transcript] = num_sites

	return StatDicts

def GeneAverage(statdicts):
	# (transcript,gene,length) -> num_sites
	
	avgdicts = {}
	running_dicts = {}
	transcript_counts = {}
	gene_max = {}
	length_max = {}

	for transcript in statdicts:

		gene = transcript[1]
		
		if gene not in running_dicts:
		
			running_dicts[gene] = [transcript[2]] + [statdicts[transcript]].copy() + [[transcript[0]]] # [length,num_sites,[transcript]]
			transcript_counts[gene] = 1

			gene_max[gene] = statdicts[transcript]
			length_max[gene] = transcript[2]

		else:
			
			gene_max[gene] = max(statdicts[transcript],gene_max[gene])
			length_max[gene] = max(transcript[2],length_max[gene])

			running_dicts[gene][1] += statdicts[transcript]
			running_dicts[gene][0] += transcript[2]
			running_dicts[gene][-1] += [transcript[0]] 
			transcript_counts[gene] += 1

	for gene in running_dicts:
		avgdicts[gene] = [":".join(running_dicts[gene][-1])] + [float(x)/transcript_counts[gene] for x in running_dicts[gene][:-1]] + [length_max[gene]] + [gene_max[gene]]
		
	return avgdicts # gene -> transcript_list,avg_length,avg_num_sites,max_length,max_num_sites

def syn_dicts(syn_csv):

	if len(syn_csv):
		with open(syn_csv,'r') as infile:
			lines = infile.readlines()[1:]

		syndict = {}

		for line in lines:
			info = line.strip().split(',')
			syndict[info[0]] = (info[1],info[2])

	else:
		syndict = None

	return syndict

def title_add(convert_dict):

	if convert_dict is not None:
		return ",symbol,fullname"
	else:
		return ""

def line_add(FBgn,convert_dict):

	if (convert_dict is not None) and (FBgn in convert_dict):
		return ",{},{}".format(convert_dict[FBgn][0],convert_dict[FBgn][1])
	else:
		return ""


def print_Gene_dict(genestats,outname,convert_dict):

	with open("{}.csv".format(outname),'w') as outfile:
		
		outfile.write("parent,transcripts,avg_length,avg_num_sites,max_length,max_num_sites{}\n".format(title_add(convert_dict)))

		for gene in genestats:
			
			outfile.write("{},".format(gene) + ",".join(str(e) for e in genestats[gene]) + "{}\n".format(line_add(gene,convert_dict)))

def print_transcript_dict(statdict,outname,convert_dict):	
	
	# (transcript,gene,length) -> num_sites

	with open("{}.csv".format(outname),'w') as outfile:

		outfile.write("transcript,parent,length,num_sites{}\n".format(title_add(convert_dict)))

		for transcript in statdict:

			outfile.write("{},{},{},".format(transcript[0],transcript[1],transcript[2]) + "{}".format(statdict[transcript]) + "{}\n".format(line_add(transcript[1],convert_dict)))

if __name__ == '__main__': # allows another python script to import the functions

	parser = argparse.ArgumentParser(description="Outputs stats for kmer prevalence in transcripts and gene averages.", usage='%(prog)s --csv',add_help=False) # a description of the function

	required = parser.add_argument_group('Required Input', 'These specifications are necessary to run.')
	required.add_argument("--csv",type=str,help="Input csvs to get stats on.", required=True)

	data_opt = parser.add_argument_group('Basic Arguments', 'These options can be used to change how the program runs.')
	data_opt.add_argument("-C","--convert",type=str,default="",help="Three column csv file with FBgn,symbol,gene name.")

	help_opt = parser.add_argument_group('Help')
	help_opt.add_argument('-h', '--help', action="help", help="show this help message and exit")

	args = parser.parse_args()

	statlist = args.csv.split(',')

	convert_dict = syn_dicts(args.convert)

	to_stat = create_dicts(statlist)

	statdict = GetStats(to_stat)

	print_transcript_dict(statdict,"TRANSCRIPT_K_STATS",convert_dict)

	genedict = GeneAverage(statdict)

	print_Gene_dict(genedict,"GENE_K_STATS",convert_dict)
