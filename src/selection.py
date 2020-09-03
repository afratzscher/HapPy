'''
FILE: selection.py
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

def autosomes(filepath):
	config.__FILENAME__ = "1000G_chr" + (str(int(config.__CHR__)) + "_" +
						str(config.__START__) + "-" + 
						str(config.__END__) + ".vcf")
	return

def createFolder(filepath):
	config.__FILEPATH__ = filepath + "/results/" + config.__GENENAME__ + "/"
	geneFolder = filepath + "/results/" + config.__GENENAME__ 
	#if no folder, create folder
	if not (os.path.isdir(geneFolder)):
		try:
		    Path(geneFolder).mkdir(parents=True, exist_ok=True)
		except OSError:
		    print ("Creation of the directory %s failed" % geneFolder)
	return

def run_commands(*commands):
	os.system(' ; '.join(commands))
    # subprocess.run(' ; '.join(commands), shell=True)
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

# for P3 grch37
def makeCommandsP3(name, ftp, cmds, xy):
	if cmds == "":
		cmds = []
	#define ftp
	cmds.append(ftp+name)

	cmds.append('bcftools view "$ftp" -I -r ' + str(config.__CHR__) + ":" 
			+ str(config.__START__) + "-" + str(config.__END__) 
			+ " > " + config.__FILEPATH__ + config.__FILENAME__ )
	# clean up
	cmds.append('rm ' + name + '.tbi')

	return cmds
def getDataP3(filepath):
	createFolder(filepath) # create folder if doesnt exist

	#if DONT have data, fetch
	if not Path(config.__FILEPATH__ + config.__FILENAME__).is_file():
		ftp = "ftp=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
		vcfgzName = "ALL.chr" + str(config.__CHR__) + ".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
		cmd = makeCommandsP3(vcfgzName, ftp, "", "")
		run_commands(*cmd)

# for P3 grch38 b149
def makeCommandsP3_38(name, ftp, cmds, xy):
	if cmds == "":
		cmds = []
	#define ftp
	cmds.append(ftp+name)

	cmds.append('bcftools view "$ftp" -I -r ' + str(config.__CHR__) + ":" 
			+ str(config.__START__) + "-" + str(config.__END__) 
			+ " > " + config.__FILEPATH__ + config.__FILENAME__ )
	# clean up
	cmds.append('rm ' + name + '.tbi')

	return cmds
def getDataP3_38(filepath):
	createFolder(filepath) # create folder if doesnt exist

	#if DONT have data, fetch
	if not Path(config.__FILEPATH__ + config.__FILENAME__).is_file():
		# get file name
		ftp = "ftp=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/"
		vcfgzName = "ALL.chr" + str(config.__CHR__) + "_GRCh38.genotypes.20170504.vcf.gz"
		cmd = makeCommandsP3_38(vcfgzName, ftp, "", "")
		run_commands(*cmd)

# for P3 grch38 b154
def makeCommandsP3_38_154(name, ftp, cmds, xy):
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
	cmds.append('bcftools view "$ftp" -I -r ' + str(config.__CHR__) + ":" 
			+ str(config.__START__) + "-" + str(config.__END__) 
			+ " > " + config.__FILEPATH__ + "raw_" + baseName)
	cmds.append('bcftools query "$dbSNPFTP" -f "%POS\t%ID\t%REF\t%ALT\n" -r ' + config.__CHRVERSION__
					+ ":" + str(config.__START__) + "-" + str(config.__END__) + "> " + config.__FILEPATH__ + "dbSNP_" + baseName)
	
	# clean up
	cmds.append('rm ' + name + '.tbi')
	cmds.append('rm ' + refVersion + '.tbi')
	return cmds
def getDataP3_38_154(filepath):
	createFolder(filepath) # create folder if doesnt exist

	#if DONT have data, fetch
	if not Path(config.__FILEPATH__ + config.__FILENAME__).is_file():
		# if dont have raw data, fetch
		rawName = "raw_" + config.__FILENAME__[len('1000G_'):] 
		if not Path(config.__FILEPATH__ + rawName).is_file():
			ftp = "ftp=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/"
			vcfgzName = "ALL.chr" + str(config.__CHR__) + "_GRCh38.genotypes.20170504.vcf.gz"
			cmd = makeCommandsP3_38_154(vcfgzName, ftp, "", "")
			run_commands(*cmd)
		
		#then combine raw and dbSNP
		combine()

# MAIN DECISION: picks which getData and makeCommands to use 
def getData(filepath):
	autosomes(filepath)
	if config.__REFVER__ == '38':
		getDataP3_38_154(filepath)
	elif config.__REFVER__ == '37':
		getDataP3(filepath)
def selectGene(filepath):
	# decide if need to fetch info or already in file
	infoFile = 'gene_info.txt'
	data = pd.read_csv(infoFile, sep="\t")
	data = data.loc[data['gene'] == (config.__GENENAME__ + "_" + config.__REFVER__)]
	data = data.loc[data['vrs'] == int(config.__REFVER__)]
	data = data.sort_values('chr', ascending=True) # X first, Y second
	data = data.reset_index(drop=True)

	if (not data.empty):
		config.__START__ = int((data.loc[0]['start']).item())
		config.__END__ = int((data.loc[0]['end']).item())
		config.__GENESTART__ = int((data.loc[0]['gene_start']).item())
		config.__GENEEND__ = int((data.loc[0]['gene_end']).item())
		config.__CHR__ = int(data.loc[0]['chr'])
		config.__CHRVERSION__ = (data.loc[0]['chrvers'])

	else: # if not found yet, fetch info
		errCode = fetch.main()
		if errCode == -1: # gene not found
			return(-1)
		else:
			data = data.append({'vrs': config.__REFVER__, 'gene': config.__GENENAME__ + "_" + config.__REFVER__,
									'chr': config.__CHR__, 'chrvers': config.__CHRVERSION__,
									 'start': config.__START__, 'end': config.__END__, 
									 "gene_start": config.__GENESTART__, "gene_end": config.__GENEEND__}, 
									 ignore_index = True)
			data.to_csv(infoFile, sep='\t', mode='a', index = False, header = False)
	config.__GENENAME__ = config.__REFVER__ + "_" + config.__GENENAME__
	# fetch data from 1000GP
	getData(filepath)

def extract():
	config.__GENENAME__ = 'ACKR1'
	config.__MOSTFREQ__ = 1 # includes most frequent
	return

def main():
	print("*****STARTING SELECTION*****")
	config.__FOLDERPATH__ = os.getcwd()[:-(len('src/'))]
	extract()
	errCode = selectGene(config.__FOLDERPATH__)
	