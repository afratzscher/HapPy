import config
import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np

def main():
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results"
	file = config.__FILEPATH__ + "meta_" + config.__FILENAME__

	df = pd.read_csv(file, sep = '\t')
	df = df.sort_values(by=['length'], ascending=False)
	df = df[['start', 'end', 'length']]
	
	df['left'] =(df['start'] - config.__GENESTART__) / 1000
	df['right'] = (df['end'] - config.__GENESTART__) / 1000 

	# df['left'] =(config.__GENESTART__ - df['start']) / 1000 
	# df['right'] = (df['end'] - config.__GENEEND__) / 1000

	num_genes = range(1,len(df.index)+1)
	print(df)
	# # num_genes = df['start']
	# print(num_genes)
	# print(df)

	plt.rcParams["figure.figsize"] = (3,7)
	plt.rcParams['font.size'] = 8
	# ax.vlines(x=20, ymin = 1, ymax = len(df.index)+1, color='black', linewidth=.2)
	plt.hlines(y=num_genes, xmin=df['left'], xmax=df['right'], color='grey', linewidth=.2)
	plt.title("Length of haplotypes (" + config.__GENENAME__ + ")") 
	plt.xlabel('Length (nt x 10^3)')
	# plt.ylabel('Gene location (nt)')
	plt.ylabel('Genes')
	plt.yticks(np.arange(100, len(df.index)+1-100, 100))
	plt.yticks(list(plt.yticks()[0]) + [len(df.index)+1, 1])
	# plt.xticks(np.arange(0, len(df.index)+1-100, 1000))


	# plt.xticks(list(plt.xticks()[0]) + extraticks)
	plt.gca().invert_yaxis()

	# # Show the graph
	plt.savefig('graph_ACKR1.png', dpi=800)
	plt.show()



if __name__ == '__main__':
	main()