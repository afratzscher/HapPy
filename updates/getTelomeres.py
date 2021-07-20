'''
FILE: getTelomeres.py
PURPOSE: get telomere and centromere sites
NOTE: FASTA starts with index 0, saves assuming index 0 for first nt
INPUT: df
OUTPUT: df with sequence added
CREATED BY: Anne-Sophie Fratzscher
'''

import os
from Bio import Entrez, SeqIO
import pandas as pd


def eFetch(chrom,start,end):
	Entrez.email = "afratzscher@yahoo.com"

	try:
		handle = Entrez.efetch(db="nucleotide", id=chrom,
			seq_start=start , seq_stop= end,
			rettype="fasta",retmode="text")
	except:
		return(None, -1)
	
	record = SeqIO.read(handle, "fasta")
	handle.close()
	return record, 0

def getFASTA(chrom):
	lst = []
	start = 0
	# start = 248000000
	# end = 249000000
	end = 1000000 # go in increments of 1 million
	flag = True
	while flag:
		# print(start, '-', end, 'vales')
		record, errCode= eFetch(chrom,start,end)
		if errCode == -1: # stops when done
			flag = False
			break
		seq = str(record.seq)
		if 'N' in seq: # have N
			# print(start, end)
			multipleFlag = True
			startidx = 0
			endidx = len(seq)
			# recursive answer
			while multipleFlag:
				tmp = seq[startidx:endidx]
				if (tmp == ''):
					break
				if (tmp.count('N') == len(tmp)): # if all N, break
					lst.append(str(startidx+start) + "-" + str(endidx-1+start))
					break

				nfirst = tmp.find("N")
				if nfirst == -1:
					startidx = endidx+1
					continue

				afirst = tmp[nfirst:].find("A")
				cfirst = tmp[nfirst:].find("C")
				tfirst = tmp[nfirst:].find("T")
				gfirst = tmp[nfirst:].find("G")
				first = nfirst
				try:
					last = min([val for val in [afirst, cfirst, tfirst, gfirst] if val >= 0]) + first
				except: #have N at end
					nlast = tmp.rfind("N")
					lst.append(str(nfirst+startidx+start) + "-" + str(endidx-1+start))
					break
				
				if last == end:
					lst.append(str(first+startidx+start) + "-" + str(endidx-1+start))
					multipleFlag = False
				else:
					lst.append(str(first+startidx+start) + "-" + str(last+startidx-1+start))
					startidx = last + startidx

			# print(lst)
			# first = seq.find('N')
			# last = seq.rfind('N')
			# print(first,last)
			# # try:
			# # 	print(seq[last-1], seq[last], seq[last+1])
			# # except:
			# # 	print(seq[first-1], seq[first], "last nt")
			# lst.append(str(first+start)+'-'+str(last+start))
		# start+=1000000
		# end+=1000000
		start+=1000000
		end+=1000000
		# print(lst)

	# print(lst)
	
	# then merge where N continues over multiple elements
	cleaned = []
	tmp_start = 0
	tmp_end = 0
	flag = False
	for i in range(0,len(lst)):
		first,last = lst[i].split("-")
		if ((int(last)%1000000) == 0) and not flag:
			tmp_start = first
			flag = True
			continue
		elif flag and ((int(last)%1000000)) == 0:
			continue
		elif flag:
			cleaned.append(tmp_start+"-"+last)
			flag = False
			continue
		else:
			cleaned.append(first+"-"+last)

	# get length of elements -> use to get centromere -> centromere = longest element
	length = []
	for i in cleaned:
		first,last = i.split("-")
		length.append(int(last)-int(first))

	idx = length.index(max(length)) #gets centromere based on longest stretch of "N"

	# first telomere in cleaned[0]
	# centromere in cleaned[idx]
	# second telomere in cleaned[-1] (last in list)

	file = open("GRCH38_chr_sites.txt", 'a')
	short = str(int(cleaned[0].split("-")[1])+1) + "-" + str(int(cleaned[idx].split("-")[0])-1) + "\t"
	long_val = str(int(cleaned[idx].split("-")[1])+1) + "-" + str(int(cleaned[-1].split("-")[0])-1) + "\t"
	telomere1 = cleaned[0] + "\t"
	centromere = cleaned[idx] + "\t"
	telomere2 = cleaned[-1] + "\n"

	num = chrom.split(".")[0].split("0")[-1]
	file.write(num+"\t"+chrom+"\t"+telomere1+short+centromere+long_val +telomere2)

def main():
	print('*****STARTING TELOMERE/CENTROMERE CALCULATION*****')
	name = os.getcwd().rsplit("/",1)[0] + "/src/GRCh38_chr_versions.txt"
	df = pd.read_csv(name, sep = '\t')

	for index,row in df.iterrows():
		chrom = row['version']
		print(chrom)
		seq = getFASTA(chrom)
	print('*****DONE WITH TELOMERE/CENTROMERE CALCULATION*****')

if __name__ == '__main__':
	main()
