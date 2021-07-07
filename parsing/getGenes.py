# TO DO: read annotRelease from file and UPDATE once in a while
## TO DO: ask for email
#NOTE: changes start/end values so that 1st nt has value of 1 (vs. stored as 1st nt has index of 0)#
from Bio import Entrez
import xmltodict
import ujson
import os
from datetime import date, time 
import pandas as pd
import getRange

# def getInfo():


'''
FUNCTION: eSearch(ext, nt)
PURPOSE: makes call to database
INPUT: ext (term to search) and nt (database to search)
OUTPUT: results from search
'''
def eSearch():
	Entrez.email = "afratzscher@yahoo.com"
	paramEutils = {'usehistory': 'Y'} # use entrez search history to cache results
	
	db = 'gene'
	term =  '"Homo sapiens"[Organism] AND ("genetype protein coding"[Properties] AND alive[prop] \
		AND (NC_000001[nucl_accn] AND 000000000001[CHRPOS] : 000248956422[CHRPOS]))'

	eSearch = Entrez.esearch(db=db, term=term, sort='Location', **paramEutils)

	res = Entrez.read(eSearch)

	paramEutils['WebEnv'] = res['WebEnv']
	paramEutils['query_key'] = res['QueryKey']
	paramEutils['rettype'] = 'xml' #report as xml
	paramEutils['retstart'] = 0

	result = Entrez.esummary(db=db, **paramEutils)
	xml = result.read()

	# convert xml to python dict object for convenient parsing
	dsdocs = xmltodict.parse(xml)
	return dsdocs

def query(chrom):
	lst = []
	# getInfo()
	dsdocs = eSearch()
	length=0
	#get set of dbVar DocumentSummary (dsdocs) and print report for each (ds)
	for ds in dsdocs ['eSummaryResult']['DocumentSummarySet']['DocumentSummary']:
		length+=1
		geneName = ds['Name']
		# print(geneName)
		if 'pseudogene' in ds['Description']:
			pseudoflag = True
		else:
			pseudoflag = False
		if 'pseudogene' in ds['Description']: # excludes some pseudogenes that are not filtered out, including OR10J4
			pseudo = 'pseudogene'
		else:
			pseudo = 'gene'
		loc = ds['MapLocation']

		#NOTE: FAM151A and CYP4X1 are on short arm, but missing annotation in file for some reason...
		if geneName == 'FAM151A':
			loc = 'p'
		if geneName == 'CYP4X1':
			loc = 'p'

		loc = ''.join(filter(str.isalpha, loc))
		loc = ''.join(set(loc)) # removes duplicate letter (sometimes get qq or pp)
		if loc == 'p': # p is short arm
			loc = 'short'
		if loc == 'q': # q is long arm
			loc = 'long'

		for p in ds['LocationHist']['LocationHistType']:
			if p['AnnotationRelease'] == annotRelease:
				start = int(p['ChrStart']) + 1 # bc 1st nt has index of 0 (instead of 1)
				end = int(p['ChrStop']) + 1
				if start > end:
					strand = 'minus'
				else:
					strand = 'plus'
				# if loc == 'long':
				# 	lst.append({'gene': geneName, 'start': start, 'end': end})
				if not pseudoflag:
					lst.append({'gene': geneName, 'chr': chrom, 'start': start, 'end': end, 'arm': loc, 'type': pseudo, 'strand': strand})
					# lst.append({'gene': geneName, 'chr': chrom, 'start': start, 'end': end, 'arm': loc})
				break
			break
	lst = sorted(lst, key = lambda i: i['start'])
	head, tail = os.path.split(os.getcwd())
	fname = head + "/src/data.json"
	with open(fname, 'w', encoding='utf-8') as f:
		ujson.dump(lst, f, ensure_ascii=False, indent=4)
							
def jsonToTxt():
	today = date.today().strftime("%B %d, %Y")
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/src"
	file = direct +  "/data.json"
	txtfile = open(direct + "/data.txt", "w")
	txtfile.write("# Updated on: " + today + "\n")
	df = pd.read_json(file)
	df.to_csv(txtfile, index = False, mode = 'a')

def main():
	print("*****STARTING QUERY*****")
	global annotRelease 
	global chrom
	annotRelease = '109.20210514'

	# currently just for chrom 1
	chrom = 1
	query(chrom)
	jsonToTxt() 
	getRange.main()


if __name__ == '__main__':
	main()
