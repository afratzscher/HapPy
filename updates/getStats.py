import pandas as pd
import os

def main():
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/src/data"
	file = direct +  "/data.json"
	df = pd.read_json(file)
	print('# genes on chr 1: ', df.shape[0])
	# longarm = df[df['arm'] == 'long']
	# shortarm = df[df['arm'] == 'short']
	# # tmp = df[df['arm'] != 'short']
	# # print(tmp)
	# # tmp2 = tmp[tmp['arm'] !=' long']
	# # print(tmp2)

	# for long arm
	longarm = df[df['arm'] == 'long']
	print('# long arm genes: ', longarm.shape[0])
	longminus = longarm[longarm['strand'] == 'minus']
	print('# long arm + minus strand genes: ', longminus.shape[0])
	longpos = longarm[longarm['strand'] == 'plus']
	print('# long arm + plus strand genes: ', longpos.shape[0])

	shortarm = df[df['arm'] == 'short']
	print('# short arm genes: ', shortarm.shape[0])
	shortminus = shortarm[shortarm['strand'] == 'minus']
	print('# short arm + minus strand genes: ', shortminus.shape[0])
	shortpos = shortarm[shortarm['strand'] == 'plus']
	print('# short arm + plus strand genes: ', shortpos.shape[0])

	
	# print(df)
	# df.to_csv(txtfile, index = False, mode = 'a')

if __name__ == '__main__':
	main()