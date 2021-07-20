import pandas as pd

def select(before, after):
	df = pd.read_json('data.json')
	idx = df[df['gene'] == 'ACKR1'].index.tolist()[0]
	genes = df.loc[idx-before:idx+after]
	print(genes)
	overlap = genes.loc[genes['overlap'] > 0].reindex().sort_index(ascending=False)
	overlap['length'] = overlap['end'] - overlap['start']
		#overlapping genes have 'True' in overlap column
		#always pick longest one (NOTE: values already correct, dont need to be changed)

	prev = None
	flag = True
	maxgene = None
	lst = []

	for index, row in overlap.iterrows():
		if flag:
			maxgene = row['gene']
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

	others = genes[~genes['overlap']]
	lst = lst + (others['gene'].tolist())
	return lst

def main():
	before = 19 # number genes before ACKR1
	after = 5 # number genes after ACKR1

	lst = select(before, after)
	print(lst)
	print(len(lst))


if __name__ == '__main__':
	main()