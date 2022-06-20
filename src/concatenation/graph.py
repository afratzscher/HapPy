import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import pandas as pd
import json

# def getVersion():
# 	versions = pd.read_csv('GRCh38_chr_versions.txt', sep = '\t')
# 	config.__CHRVERSION__ = versions.loc[versions['chr'] == str(config.__CHR__)]['version'].values[0]


# def getRange():
# 	df = pd.read_json('updated_data.json')

# 	# swap start/end for minus
# 	if config.__GENENAME__ in df.gene.values:
# 		idx = df[df['gene'] == config.__GENENAME__].index.tolist()[0]
# 		gene = df.loc[idx]


def main():
	direct = os.path.split(os.getcwd())[0]
	ranges = pd.read_json(direct + '/src/data/chr1.json')
	genes = ['CD1D', 'CD1A', 'CD1B', 'CD1E', 'OR10T2', 
		'OR10K2', 'OR10K1', 'OR10R2', 'OR6Y1', 'OR6P1', 
		'OR10X1', 'SPTA1', 'OR6K2', 'OR6K3', 'OR6K6', 
		'OR6N1', 'PYHIN1', 'IFI16', 'AIM2', 'CADM3', 
		'ACKR1']
	start = []
	end = []
	for gene in genes:
		idx = ranges[ranges['gene'] == gene].index.tolist()[0]
		vals = ranges.loc[idx]
		start.append(vals['start'])
		end.append(vals['end'])

	num_genes = range(1,len(genes)+1)
	colours = ['green', 'purple', 'blue', 'cyan', 'grey', 'purple','green', 'olive', 'blue', 'cyan', 'green',
		'green', 'purple', 'blue', 'cyan', 'grey', 'purple','green', 'olive', 'blue', 'cyan', 'red']
	print(colours)

	minval = 158100263
	maxval = 159712288

	plt.hlines(y=0, xmin=minval, xmax=maxval, color='black')
	plt.hlines(y=0, xmin=minval, xmax=maxval, color='black', linewidth=1000, alpha=0.2)
	for i in range(0, len(num_genes)):
		if start[i] <= maxval:
			plt.hlines(y=0, xmin=start[i], xmax=end[i], color=colours[i], linewidth=20, label = genes[i])
		else:
			print('here')
			plt.hlines(y=0, xmin=start[i], xmax=end[i], color=colours[i], linewidth=20, label = genes[i])
	plt.tick_params(left = False, right = False , labelleft = False ,
                labelbottom = True, bottom = True)
	plt.xticks(list(plt.xticks()[0][2:]) + [minval, maxval], fontsize = 'small')
	plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')
	plt.ticklabel_format(useOffset=False, style='plain')
	plt.show()

	# # plt.show()
	# # # Show the graph
	# visualizationFolder = (config.__FOLDERPATH__ + "/results/" 
	# 	+ config.__FOLDERNAME__ + "/visualization/")
	# graphname = config.__GENENAME__ + '_haplotypes_graph_vrs02.png'
	# plt.savefig(visualizationFolder + graphname, dpi=800)

if __name__ == '__main__':
	main()