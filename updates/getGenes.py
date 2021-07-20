# TO DO: read annotRelease from file and UPDATE once in a while
## TO DO: ask for email
#NOTE: changes start/end values so that 1st nt has value of 1 (vs. stored as 1st nt has index of 0)#
from Bio import Entrez
import xmltodict
import ujson
import os
from datetime import date, time 
import pandas as pd
from pathlib import Path

def getInfo(chrom):
	if chrom == 1:
		val = 'NC_000001'
		start = 1
		end = 248956422
	else:
		print("TO ADD")
	return val, start, end

def eSearch(chrom, start, end):
	Entrez.email = "anne-sophie.fratzscher@mail.mcgill.ca"
	paramEutils = {'usehistory': 'Y'} # use entrez search history to cache results
	
	db = 'gene'
	term =  '"Homo sapiens"[Organism] AND ("genetype protein coding"[Properties] AND alive[prop] \
		AND (' + chrom + '[nucl_accn] AND ' + str(start) + '[CHRPOS] : ' + str(end) + '[CHRPOS]))'

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

def countOverlap(df):
	# swap minus start/end so that larger number is end
	minus = df[df['strand'] == 'minus'].index.tolist()
	df.loc[minus,['start','end']] = df.loc[minus,['end','start']].values
			# swaps start/end for genes on minus strand (so that can find region on pos strand and not have error)

	#sort by start b/c inverted minus strand genes
	df = df.sort_values('start').reset_index(drop=True)
	df['num_overlap'] = 0

	for row in df.itertuples():
	    mask = (row.start <= df.end) & (row.end >= df.start)
	    df.loc[row.Index, 'num_overlap'] = sum(mask) - 1

	return df

def query(chrom):
	lst = []
	val, start, end = getInfo(chrom) # should store chrom val, start, end
	dsdocs = eSearch(val, start, end)
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

	# convert to df
	df = pd.DataFrame(lst)	
	df = countOverlap(df)
	df.insert(loc=0, column='gene_num', value=range(1, len(df)+1))
	df = df.sort_values('start').reset_index(drop=True)
	df = df.sort_values('end').reset_index(drop=True)
	df = df.sort_values('start').reset_index(drop=True)
	return df

# def toCSV(df, file):
# 	today = date.today().strftime("%B %d %Y")

# 	df = df.replace(False, "")
# 	name = file[:-5] + '.csv'
# 	f = open(name, "w", newline='')
# 	f.write("# Updated on: " + today + "\n")
# 	f.write("# NOTE: short arm = nt 10000-125184587; long arm = nt 143184588-248946422 \n")

# 	# remove false and replace with None for simplicity when reading for overlapping genes
# 	df.to_csv(f, index = None)
# 	f.close()

def rangeCalc(df, first, last, longFlag):
	startlist = []
	endlist = []

	# go through all genes and get range
	for idx in range(0, len(df.index)):
		rangeEnd = 0
		rangeStart = 0

		df = df.sort_values('start').reset_index(drop=True)
		gene = df.loc[idx]
			#sort by start (since doesnt work completely in getGenes.py b/c negative genes have flipped start/end)
		if gene['skip']:
			rangeEnd = None
		else:
			if idx == (len(df.index)-1): 
				rangeEnd = last  #CHANGE LATER!
			else: # if not last gene
				nextGene = df.loc[idx+1]
				# check if gene overlaps with next gene -> if so, increase idx
				i = idx
				while nextGene['skip']:
					i+=1
					if i > len(df.index)-1:
						endFlag = True
						break
					nextGene = df.loc[i]
				rangeEnd = nextGene['start'] -1

		df = df.sort_values('end').reset_index(drop=True)
		idx = df[df['gene'] == gene['gene']].index.values[0]
		if gene['skip']:
			rangeStart = None
		else:
			if idx == 0: # if not first gene
				rangeStart = first  
			else:
				prevGene = df.loc[idx-1]
				# check if gene overlaps with previous gene -> if so, decrease idx
				i = idx
				while prevGene['skip']:
					i-=1
					if i < 0:
						beginFlag = True
						break
					prevGene = df.loc[i]
				rangeStart = prevGene['end'] +1
		endlist.append(str(rangeEnd))	
		startlist.append(str(rangeStart))

	df = df.sort_values('start').reset_index(drop=True)
	df['range_start'] = startlist
	df['range_end'] = endlist
	return df


def getRange(df):
	# # swap minus start/end so that larger number is end
	# minus = df[df['strand'] == 'minus'].index.tolist()
	# df.loc[minus,['start','end']] = df.loc[minus,['end','start']].values
	# 		# swaps start/end for genes on minus strand (so that can find region on pos strand and not have error)

	#sort by start b/c inverted minus strand genes
	df = df.sort_values('start').reset_index(drop=True)

	# do for short arm first
	shortarm =  df[df['arm'] == 'short']
	shortarm = rangeCalc(shortarm, 10000, 125184587, False)
	file = direct + "/chr1_short.json"
	# shortarm.insert(loc=0, column='gene_num', value=range(1, len(shortarm)+1))
	shortarm.to_json(file, orient='records', indent=4)
	# toCSV(shortarm, file)

	# then do long arm
	longarm =  df[df['arm'] == 'long'].reset_index(drop=True)
	longarm = rangeCalc(longarm, 143184588, 248946422, True)
	file = direct + "/chr1_long.json"
	# longarm.insert(loc=0, column='gene_num', value=range(1, len(longarm)+1))
	longarm.to_json(file, orient='records', indent=4)
	# toCSV(longarm, file)

	#have file with both arms together
	# del longarm['gene_num']
	# del shortarm['gene_num']
	df = shortarm.append(longarm)
	df = df.sort_values('start').reset_index(drop=True)

	# check for overlap
	# df['start_overlap'] = df['range_start'].duplicated(keep=False) # True if have duplicate, False if no duplicate
	# df['end_overlap'] = df['range_end'].duplicated(keep=False)
	# # if True for start and end, then is overlapping gene
	# df['overlap'] = (df['start_overlap']) | (df['end_overlap'])
	# del df['start_overlap']
	# del df['end_overlap']

	file = direct + "/chr1.json"
	df.to_json(file, orient='records', indent=4)
	return df
	# toCSV(df, file)

def createFolder():
	folder = direct
	if not (os.path.isdir(direct)):
		try:
		    Path(direct).mkdir(parents=True, exist_ok=True)
		except OSError:
		    print ("Creation of the directory %s failed" % direct)

def selectOverlapGenes(df):
	overlap = df.loc[df['num_overlap'] > 0].reindex().sort_index(ascending=False)
	names = overlap['gene'].tolist()
	overlap['length'] = overlap['end'] - overlap['start']
		#overlapping genes have 'True' in overlap column
		#always pick longest one (NOTE: values already correct, dont need to be changed)
	prev = None
	flag = True
	maxgene = None
	loops = -1
	lst = []
	for index, row in overlap.iterrows():
		if flag: # if first item in df
			maxgene = row['gene']
			loops = row['num_overlap']
			prev = row
			flag = False
			continue
		if (prev['gene_num']-row['gene_num']) == 1:
			if row['length'] > prev['length']:
				maxgene = row['gene']
		else:
			lst.append(maxgene)
			maxgene = row['gene']
		prev = row
	lst.append(maxgene)

	toExclude = [x for x in names if x not in lst] 
	
	#remove range for genes that are excluded
	df['range_start'] = -1
	df['range_end'] = -1
	df.loc[df['gene'].isin(toExclude), 'range_start'] = None
	df.loc[df['gene'].isin(toExclude), 'range_end'] = None

	# True in skip if gene overlaps but is shorter (so dont use in analysis)
	df['skip']=df['gene'].apply(lambda x: True if x in toExclude else False)
	del df['num_overlap']
	return df

def main():
	print("*****STARTING QUERY*****")
	global annotRelease 
	global chrom
	global direct
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/src/data"

	createFolder()
	annotRelease = '109.20210514'

	# currently just for chrom 1
	chrom = 1
	df = query(chrom)
	df = selectOverlapGenes(df)
	df = getRange(df)

	file = direct + "/chr1.json"
	df.to_json(file, orient='records', indent=4)


if __name__ == '__main__':
	main()