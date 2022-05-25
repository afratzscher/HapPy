import pandas as pd
import os

def select(before, after, start_gene):
	i = before
	j = after
	direct = "/Users/afratzscher/Documents/GitHub/HapPy/src/"
	df = pd.read_json(direct + '/data/chr1.json')
	idx = df[df['gene'] == start_gene].index.tolist()[0]
	# print(idx)

	flag = True
	while flag:
		prevgenes = df.loc[idx-i:idx-1]
		prevgenes = prevgenes.loc[prevgenes['skip'] == False]
		# print(prevgenes)
		# print(len(prevgenes))
		# print(before)
		# exit()
		if len(prevgenes) < before:
			i+=1
		else:
			flag = False

	flag = True
	while flag:
		nextgenes = df.loc[idx+1:idx+j]
		nextgenes = nextgenes.loc[nextgenes['skip'] == False]
		if len(nextgenes) < after:
			j+=1
		else:
			flag = False

	lst = prevgenes['gene'].tolist()
	lst.append(start_gene)
	lst += (nextgenes['gene'].tolist())
	return lst[::-1]

def main():
	start_gene = 'H3-7'
	before = 0
	after = 832
	# 948 genes (including those to skip)
	# 833 genes (exclude those to skip, including pseudogenes)

	lst = select(before, after, start_gene)
	# print(lst)
	print(len(lst))
	print(lst[::-1])

if __name__ == '__main__':
	main()