'''
FILE: fetch.py
PURPOSE: finds location of gene of interest; then, finds closest protein-encoding genes (NOT pseudogenes) to find range to search
INPUT: none
OUTPUT: error code (-1 if error)
CREATED BY: Anne-Sophie Fratzscher
NOTES: Implementation based on https://www.ncbi.nlm.nih.gov/dbvar/content/tools/entrez/
'''
import config
import pandas as pd
import sys
from Bio import Entrez
import xmltodict

'''
FUNCTION: setConfigValues(p)
PURPOSE: sets config values based on info found either in db or in file
INPUT: p (row of info from file or db)
OUTPUT: none
'''
def setConfigValues(p):
	config.__CHRVERSION__ = p['ChrAccVer']
	chrom = ((config.__CHRVERSION__[4:].lstrip("0")).split(".", 1))[0]
	config.__CHR__ = int(chrom)
	start = int(p['ChrStart']) + 1
	end = int(p['ChrStop']) + 1
	if start > end: # need to do this check because sometimes the values are flipped?
		temp = start
		start = end
		end = temp
	config.__GENESTART__ = int(start)
	config.__GENEEND__ = int(end)

'''
FUNCTION: eSearch(ext, nt)
PURPOSE: makes call to database
INPUT: ext (term to search) and nt (database to search)
OUTPUT: results from search
'''
def eSearch(ext, nt):
	Entrez.email = config.__EMAIL__
	paramEutils = {'usehistory': 'Y'} # use entrez search history to cache results
	
	if nt == 'nucleotide':
		db = 'nucleotide'
		term = ext
	else:
		db = 'gene'
		term = '("homo sapiens"[Organism]) AND ' + ext

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

'''
FUNCTION: getGeneRange()
PURPOSE: finds start and end nucelotide for gene
INPUT: none
OUTPUT: error code (if cnnot be found)
'''
def getGeneRange():
	dsdocs = eSearch(config.__GENENAME__, "")

	fail = False
	errCode = 0

	try:
		dsdocs['eSummaryResult']['DocumentSummarySet']
	except KeyError: #search failed -> return that gene not found
		fail = True
		errCode = 1
		return(errCode)

	if not fail:
		#get set of dbVar DocumentSummary (dsdocs) and print report for each (ds)
		for ds in dsdocs ['eSummaryResult']['DocumentSummarySet']['DocumentSummary']:
			if ds['Name'] == config.__GENENAME__:
				for p in ds['LocationHist']['LocationHistType']:
					if p['ChrAccVer'][:7] == 'NC_0000':
						if (config.__REFVER__ == '37'): #want first 105 (since is sorted)
							config.__ANNOTATIONRELEASE__ = '105'
							if (p['AnnotationRelease'] == '105'):
								setConfigValues(p)
						elif (config.__REFVER__ == '38'):
							config.__ANNOTATIONRELEASE__ = p['AnnotationRelease']
							setConfigValues(p)
							
	return(errCode)

'''
FUNCTION: getLength()
PURPOSE: gets length of chromosome
INPUT: none
OUTPUT: length of chromosome
'''
def getLength():
	chromV = config.__CHRVERSION__
	chromosome = config.__CHR__
	
	if config.__REFVER__ == '38': # search nucleotide database
		dsdocs = eSearch(chromV, 'nucleotide')
		for ds in dsdocs['eSummaryResult']['DocSum']['Item']:
			if (ds['@Name'] == 'Length'):
				length = ds['#text']
		return(int(length))
	elif config.__REFVER__ == '37': # read from file
		versionFile = 'GRCh37_version_length.txt'
		data = pd.read_csv(versionFile, sep="\t")
		length = data.iloc[int(chromosome)-1]['length']
		return(int(length))

'''
FUNCTION: getClosest(upDown, length)
PURPOSE: gets length of chromosome
INPUT: upDown (look upstream of downstream of gene), length of chromosome
OUTPUT: none
'''
def getClosest(upDown, length): # looks for nearest protein-coding gene (NOT pseudogene)
	increment = 100000 # default is search 100,000 bases away

	check = 0
	while (check == 0): # if not found, search further
		if upDown == 'up':
			start = config.__GENESTART__ - increment
			end = config.__GENESTART__
		elif upDown == 'down':
			start = config.__GENEEND__
			end = config.__GENEEND__ + increment

		if (start < 1):
			config.__START__ = 1
			return
		elif (end > length):
			config.__END__ = int(length)
			return
		df = pd.DataFrame(columns=['name', 'start', 'end'])
		chromosome = config.__CHRVERSION__ 
		ext = (chromosome.split(".", 1)[0]
			+ "[nucl_accn] AND (" + str(start) + "[CHRPOS]"
			+ " : " + str(end) + "[CHRPOS]" + ') AND ("genetype protein coding"[Properties])')
		dsdocs = eSearch(ext, "")
		
		fail = False
		try:
			dsdocs['eSummaryResult']['DocumentSummarySet']
		except KeyError: #search failed (nothing found)
			fail = True
			increment += 100000

		successful = False
		if not fail:
			try: # only 1 gene found (edge case)
				name = dsdocs['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['Name']
				successful = True
			except TypeError: # multiple genes found
				for ds in dsdocs ['eSummaryResult']['DocumentSummarySet']['DocumentSummary']:
					name = ds['Name']
					for p in ds['LocationHist']['LocationHistType']:	
						if (p['ChrAccVer'][:2] == 'NC'):
							if (p['AnnotationRelease'] == config.__ANNOTATIONRELEASE__):
								if (p['ChrAccVer'] == config.__CHRVERSION__):
									fStart = int(p['ChrStart']) + 1
									fEnd = int(p['ChrStop']) + 1
									if fStart > fEnd:
										temp = fStart
										fStart = fEnd
										fEnd = temp
									if (fStart > start):
										df = df.append({'name': name, 'start': fStart, 'end' : fEnd}, ignore_index = True)
			if successful:
				for p in dsdocs ['eSummaryResult']['DocumentSummarySet']['DocumentSummary']['LocationHist']['LocationHistType']:	
					if (p['ChrAccVer'][:2] == 'NC'):
						if (p['AnnotationRelease'] == config.__ANNOTATIONRELEASE__):
							if (p['ChrAccVer'] == config.__CHRVERSION__):
								fStart = int(p['ChrStart']) + 1
								fEnd = int(p['ChrStop']) + 1
								if fStart > fEnd: # sometimes end and start swapped in db
										temp = fStart
										fStart = fEnd
										fEnd = temp
								if (fStart > start):
									df = df.append({'name': name, 'start': fStart, 'end' : fEnd}, ignore_index = True)
			
			# remove gene from df if in list
			df = df[df.name != config.__GENENAME__]
			
			if (df.empty): # if empty, continue searching
				increment+=increment #doubles increment
				df = df.drop(df.index, inplace=True) # clear data frame at end
			else: # else, finish searching
				check = 1
	
	if upDown == 'up':
		df = df.sort_values('end', ascending=False)
		df = df.reset_index(drop = True)
		config.__START__ = int((df['end'][0]) + 1)
	elif upDown == 'down':
		df = df.sort_values('start')
		df = df.reset_index(drop = True)
		config.__END__ = (df['start'][0]) - 1
	df = df.drop(df.index, inplace=True) # clear data frame at end

def main():
	print("*****STARTING FETCHING*****")
	try:
		errCode = getGeneRange() #gets range of gene
	except TypeError:
		return(-1)
	if (errCode == 0):
		length = getLength()
		getClosest('up', length) # gets upstream location
		getClosest('down', length)
		return(0)
	else:
		return(-1)
