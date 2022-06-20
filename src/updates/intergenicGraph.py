import pandas as pd
import os
import matplotlib.pyplot as plt

def main():
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/src"
	file = direct +  "/chr1_long.json"
	df = pd.read_json(file)
	df = df[['gene', 'start', 'end', 'range_start', 'range_end']]
	df['left'] =(df['range_start'] - df['start']) / 1000 
	df['right'] = (df['range_end'] - df['end']) / 1000

	num_genes = range(1,len(df.index)+1)
	# num_genes = df['start']
	print(num_genes)
	print(df)
	plt.hlines(y=num_genes, xmin=df['left'], xmax=df['right'], color='black', linewidth=.2)
	plt.title("Length of intergenic regions, chromosome 1 long arm") 
	plt.xlabel('Length (nt) x 10^3')
	# plt.ylabel('Gene location (nt)')
	plt.ylabel('Genes')

	# Show the graph
	plt.savefig('intergenic_graph_vrs02.png', dpi=800)
	plt.show()

if __name__ == '__main__':
	main()