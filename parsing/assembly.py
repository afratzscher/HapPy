# TO DO: read annotRelease from file and UPDATE once in a while
## TO DO: ask for email
from Bio import Entrez
import xmltodict
import ujson


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
	
	db = 'nucleotide'
	term =  '"Homo sapiens chromosome 21, GRCh38.p13 Primary Assembly"'

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
	# print(dsdocs)
	return dsdocs

def query():
	lst = []
	# getInfo()
	dsdocs = eSearch()
	for ds in dsdocs['eSummaryResult']:
		print(ds)
		print('next')
	# print(dsdocs['eSummaryResult']['DocSum'][3])
def main():
	print("*****STARTING QUERY*****")
	query()

if __name__ == '__main__':
	main()
