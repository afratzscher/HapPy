'''
FILE: main.py
PURPOSE: main entry point to application
INPUT: none
OUTPUT: none
'''

import config
from optparse import OptionParser
import search
import pandas as pd
# import time
import fetch
import cleaner
import haplotypes
import distinct
import visualization
import popcounts
import userfile
import haplotypeGraph
import os

def getOptions():
	parser = OptionParser()
	parser.add_option("-g", "--gene", dest="gene", default=False, type='string',
					  help="select by gene", metavar="GENE")    
	parser.add_option("-r", "--region", dest="region", default=False, type='string',
					  help="select by region (GRCh38), input as chr#:#-#", metavar="chr#:#-#")
	parser.add_option("-c", "--custom", dest="foldername", default=False, type='string',
					  help="select a custom folder name (otherwise, is chr#:#-#)")
	parser.add_option("-f", "--filter", dest="filter", default=False, action="store_true",
					  help="filters data to only include original 2504 individuals")
	parser.add_option("-d", "--data", dest="data", default=False, type='string',
					  help="location of dbSNP and 1000 Genome project data ")
	(opts, args) = parser.parse_args()
	return opts

def selectionOptions(opts):
	option_dict = vars(opts)

	if option_dict.get('filter'):
		config.__FILTERED__ = True

	if option_dict.get('gene') and option_dict.get('region'):
		print("Select either gene or region (not both)")
		exit(0)

	if option_dict.get('foldername'):
		config.__FOLDERNAME__ = option_dict.get('foldername')

	# sets path where dbSNP and 1000G data have been downloaded (IF NOT SET, CALLS DIRECTLY ->
	# HAD ISSUES WITH DIRECT CALLS SINCE SEPTEMBER 2021)
	if option_dict.get('data'):
		config.__LOCAL__ = True
		config.__LOCALPATH__ = option_dict.get('data')
		
	# if given gene, search for region in gene database
	if option_dict.get('gene'):
		config.__GENENAME__ = option_dict.get('gene')
		if config.__FOLDERNAME__ == '':
			config.__FOLDERNAME__ = config.__GENENAME__
		errCode = search.main()
		if errCode == -1: # gene not found
			return(-1)
		return(0)

	elif option_dict.get('region'):
		config.__REGIONFLAG__ = True
		if not 'chr' in option_dict.get('region'):
			print("Incorrect format. Use '-h' to get help")
			exit(0)
		string = option_dict.get('region').split(':')
		config.__CHR__ = '11'#string[0][3]
		string = string[1].split("-")
		if(string[0]=='' or string[1]==''):
			print("Incorrect format. Use '-h' to get help")
			exit(-1)

		# CURRENTLY: looks at whole region as gene and as start/end
		config.__START__ = int(string[0])
		config.__END__ = int(string[1])
		config.__GENESTART__ = int(string[0])
		config.__GENEEND__ = int(string[1])

		if config.__FOLDERNAME__ == '':
			config.__FOLDERNAME__ = 'chr'+config.__CHR__+"_"+str(config.__START__)+"-"+str(config.__END__)

		df = pd.read_csv(os.getcwd() + '/data/GRCh38_chr_versions.txt', sep='\t')
		df = df[df['chr'] == config.__CHR__]
		config.__CHRVERSION__ = df['version'][0]
		return(0)
	return(-1)

def execute():
	print('Run info')
	print(config.__CHR__, 'chr')
	print(config.__START__, 'start')
	print(config.__GENESTART__, 'gene start')
	print(config.__GENEEND__, 'gene end')
	print(config.__END__, 'end')
	print(config.__GENENAME__)
	print(config.__CHRVERSION__)
	print(config.__FILENAME__)

	fetch.main() #KEEP
	cleaner.main() #KEEP
	haplotypes.main()
	distinct.main() #KEEP
	popcounts.main()
	visualization.main() #KEEP graph of populations for most frequent haplotype
	userfile.main() # need to get metadata for haplotypeGraph
	haplotypeGraph.main() #KEEP graph showing length of haplotyps
	 
def main():	
	# get options from user (gene name or region)
	opts = getOptions() 

	# parse options
	errcode = selectionOptions(opts) 

	# decide if run once or for multiple genes
	if (errcode == -1):
		print("Incorrect format. Use '-h' to get help")
	elif (errcode == 0): 
		execute()
	elif (errcode == 2):
		print('ERROR: unexpected error')
		exit()
	
if __name__ == '__main__':
	main()
