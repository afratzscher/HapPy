import pandas as pd
import os

def select(before, after):
	i = before
	j = after
	direct = os.getcwd()
	df = pd.read_json(direct + '/data/chr1.json')
	idx = df[df['gene'] == 'ACKR1'].index.tolist()[0]

	flag = True
	while flag:
		prevgenes = df.loc[idx-i:idx-1]
		prevgenes = prevgenes.loc[prevgenes['skip'] == False]
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
	lst.append('ACKR1')
	lst += (nextgenes['gene'].tolist())
	return lst[::-1]

def main():
	before = 50 # number genes before ACKR1
	after = 5 # number genes after ACKR1

	lst = select(before, after)
	print(lst)
	print(len(lst))
	print(lst[::-1])

if __name__ == '__main__':
	main()