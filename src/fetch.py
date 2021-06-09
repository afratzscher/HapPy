'''
FILE: fetch.py
PURPOSE: connects to 1000GP ftp and fetches data based on range computed in fetch.py
INPUT: none
OUTPUT: error code (-1 if error)
CREATED BY: Anne-Sophie Fratzscher
'''

import config
import os
import pandas as pd
from pathlib import Path
import subprocess
import read_vcf
import fetch

def checkAutosomes(filepath):
	config.__FILENAME__ = "1000G_chr" + (str(int(config.__CHR__)) + "_" +
						str(config.__START__) + "-" + 
						str(config.__END__) + ".vcf")
	return

def createFolder(filepath):
	config.__FILEPATH__ = filepath + "/results/" + config.__FOLDERNAME__ + "/"
	geneFolder = filepath + "/results/" + config.__FOLDERNAME__ 
	#if no folder, create folder
	if not (os.path.isdir(geneFolder)):
		try:
		    Path(geneFolder).mkdir(parents=True, exist_ok=True)
		except OSError:
		    print ("Creation of the directory %s failed" % geneFolder)
	return

def run_commands(*commands):
	os.system(' ; '.join(commands))

def combine():
	raw = read_vcf.read_vcf(config.__FILEPATH__+'raw_chr' + str(config.__CHR__) + "_" + str(config.__START__) + "-" 
			+ str(config.__END__) + ".vcf")
	chrom = str(config.__CHR__)

	dbSNP = pd.read_csv(config.__FILEPATH__+'dbSNP_chr'+ chrom + "_" + str(config.__START__) + "-" 
		+ str(config.__END__) + ".vcf", sep='\t', header=None) # ISSUE HERE
	dbSNP.columns = ['POS', 'ID', 'REF', 'ALT']

	for i in range(0, len(raw.index)):
		# if NOT in dbSNP, dont replace (remove later)
		if not (dbSNP.loc[dbSNP['POS'] == raw.loc[i]['POS']]).empty: # check if pos match
			row = dbSNP.loc[dbSNP['POS'] == raw.loc[i]['POS']].values[0]
			if ((raw.iloc[i,3] == row[2]) and (raw.iloc[i,4] == row[3])): # only replace if REF/ALT match
				raw.iloc[i,2] = row[1] # ID
				raw.iloc[i,3] = row[2] # REF
				raw.iloc[i,4] = row[3] # ALT
	df = raw[raw.ID != '.'] # remove if not in dbSNP
	df.to_csv((config.__FILEPATH__ + config.__FILENAME__), sep="\t", mode='a', index=False)

	#remove files
	# cmds = []
	# baseName =  config.__FILENAME__[len('1000G_'):] 
	# cmds.append('rm ' + config.__FILEPATH__ + "raw_" + baseName)
	# cmds.append('rm ' + config.__FILEPATH__ + "dbSNP_" + baseName)
	# run_commands(*cmds)
	

# for P3 grch38 b154
def makeCommands(name, ftp, cmds, xy):
	if cmds == "":
		cmds = []
	#define ftp
	cmds.append(ftp+name)

	#define dbSNP ftp -> GRCH38 BUILD 154 
	baseUrl = 'dbSNPFTP=ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/'
	refVersion = 'GCF_000001405.38.gz'
	cmds.append(baseUrl+refVersion)

	#command to get data file
	baseName =  config.__FILENAME__[len('1000G_'):] 
	cmds.append('bcftools view "$ftp" -I -r chr' + str(config.__CHR__) + ":" 
			+ str(config.__START__) + "-" + str(config.__END__) 
			+ " > " + config.__FILEPATH__ + "raw_" + baseName)
	cmds.append('bcftools query "$dbSNPFTP" -f "%POS\t%ID\t%REF\t%ALT\n" -r ' + config.__CHRVERSION__
					+ ":" + str(config.__START__) + "-" + str(config.__END__) + "> " + config.__FILEPATH__ + "dbSNP_" + baseName)
	# clean up
	cmds.append('rm ' + name + '.tbi')
	cmds.append('rm ' + refVersion + '.tbi')
	return cmds

def fetchSeq(filepath):
	createFolder(filepath) # create folder if doesnt exist

	#if DONT have data, fetch
	if not Path(config.__FILEPATH__ + config.__FILENAME__).is_file():
		# if dont have raw data, fetch
		rawName = "raw_" + config.__FILENAME__[len('1000G_'):] 
		if not Path(config.__FILEPATH__ + rawName).is_file():
			# ftp = "ftp=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/"
			# vcfgzName = "ALL.chr" + str(config.__CHR__) + "_GRCh38.genotypes.20170504.vcf.gz"
			ftp = "ftp=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/"
			vcfgzName = "CCDG_14151_B01_GRM_WGS_2020-08-05_chr" + str(config.__CHR__) + ".filtered.shapeit2-duohmm-phased.vcf.gz"
			cmd = makeCommands(vcfgzName, ftp, "", "")
			run_commands(*cmd)
		
		#then combine raw and dbSNP
		combine()

# MAIN DECISION: picks which getData and makeCommands to use 
def getData(filepath):
	checkAutosomes(filepath)
	fetchSeq(filepath)

def selectGene(filepath):
	# fetch data from 1000GP
	getData(filepath)

def main():
	print("*****STARTING SELECTION*****")
	config.__FOLDERPATH__ = os.getcwd()[:-(len('src/'))]
	errCode = selectGene(config.__FOLDERPATH__)
	