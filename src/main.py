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
import time
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
					  help="select by region (GRCh38)", metavar="chr#:#-#")
	parser.add_option("-a", "--all", dest="all", default=False, action="store_true",
					  help="run for all of chr 1")
	parser.add_option("-c", "--custom", dest="foldername", default=False, type='string',
					  help="select a custom folder name (otherwise, is chr#:#-#)")
	parser.add_option("-f", "--filter", dest="filter", default=False, action="store_true",
					  help="filters data to only include original 2504 individuals")
	
	(opts, args) = parser.parse_args()

	return opts

def selectionOptions(opts):
	option_dict = vars(opts)

	if option_dict.get('filter'):
		config.__FILTERED__ = True
	if option_dict.get('all'):
		print('run for all of chromosome 1')
		print("TO IMPLEMENT (settings.py)")
		return(2)

	if option_dict.get('gene') and option_dict.get('region'):
		print("Select either gene or region (not both)")
		exit(0)

	if option_dict.get('foldername'):
		config.__FOLDERNAME__ = option_dict.get('foldername')

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
		print("TO FIX -> add gene AND whole region as vars...")
		print("TO DO: also fix haplotype (before/after DNE b/c same as gene)")
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
	print('data notes')
	print(config.__CHR__, 'chr')
	print(config.__START__, 'start')
	print(config.__GENESTART__, 'gene start')
	print(config.__GENEEND__, 'gene end')
	print(config.__END__, 'end')
	print(config.__GENENAME__)
	print(config.__CHRVERSION__)
	print(config.__FILENAME__)

	start_time = time.time()
	fetch.main() #KEEP
	cleaner.main() #KEEP
	haplotypes.main()
	distinct.main() #KEEP
	popcounts.main()
	visualization.main()
	userfile.main()
	haplotypeGraph.main()
	print("--- total time %s seconds ---" % (time.time() - start_time))
			 
def main():	
	opts = getOptions()
	start_time = time.time()
	errcode = selectionOptions(opts) 
	print("--- settings %s seconds ---" % (time.time() - start_time))

	# TO REMOVE LATER 
	config.__LOCAL__ = True
	print(config.__LOCAL__, ' main run local or not')

	# decide if run once or for multiple genes
	if (errcode == -1):
		print("Incorrect format. Use '-h' to get help")
	elif (errcode == 0): # single gene
		print('here')
		execute()
	elif (errcode == 2): # all genes
		print('multi')
		execute()
	
if __name__ == '__main__':
	main()
