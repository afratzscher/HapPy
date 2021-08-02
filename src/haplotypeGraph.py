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

	mask = (df['end'] == config.__END__) & (df['start'] == config.__START__)

	df['left'] = (df['start'] - min(df['start'])) / 1000
	df['right'] = (df['end'] - min(df['start'])) / 1000

	# df['left'] = (df['start'] - config.__GENESTART__) / 1000
	# df['right'] = (df['end'] - config.__GENESTART__) / 1000 
	# df['left'] =(df['start']) / 1000
	# df['right'] = (df['end']) / 1000 

	full = df[mask]
	rest = df[~mask]
	print(full.shape)
	print("FULL***")

	num_genes = range(1,len(df.index)+1)
	num_full = range(1, len(full.index)+1)
	num_rest = range(len(full.index)+1, len(df.index)+1)

	plt.rcParams["figure.figsize"] = (3, 8)
	plt.rcParams['font.size'] = 6

	# adding area of gene
	fig, ax = plt.subplots()
	ax.axvspan((config.__GENESTART__ - min(df['start']))/1000, 
		(config.__GENEEND__ - min(df['start']))/1000, alpha=0.2, color='#FF9933')
	
	plt.hlines(y=num_full, xmin=full['left'], xmax=full['right'], color='#FF9933', linewidth=.2)
	plt.hlines(y=num_rest, xmin=rest['left'], xmax=rest['right'], color='grey', linewidth=.2)
	plt.title("Length of haplotypes (" + config.__GENENAME__ + ")") 
	plt.xlabel('Length (nt x 10^3)')
	plt.ylabel('Haplotypes')
	plt.yticks(np.arange(100, len(df.index)+1-100, 100))

	if len(df.index)+1 > 50:
		plt.yticks(list(plt.yticks()[0]) + [len(df.index)+1, 50])
	else:
		plt.yticks(list(plt.yticks()[0]) + [len(df.index)+1])
	plt.gca().invert_yaxis()

	# plt.xlim([0, max(df['right'])+10])
	# plt.ylim([-10, len(df.index)+2])

	# plt.show()
	# # Show the graph
	visualizationFolder = (config.__FOLDERPATH__ + "/results/" 
		+ config.__FOLDERNAME__ + "/visualization/")
	graphname = config.__GENENAME__ + '_haplotypes_graph_vrs02.png'
	plt.savefig(visualizationFolder + graphname, dpi=800, transparent=True)

if __name__ == '__main__':
	main()