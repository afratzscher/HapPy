import gffpandas.gffpandas as gffpd

def main():
	# print('nef')
	# annotation = gffpd.read_gff3('GRCh38_latest_genomic.gff')
	# combined_df = annotation.filter_feature_of_type(['gene']).to_gff3('genes.gff')
	# gff_genes = open('genes.gff').read()
	# print(gff_genes)
	genes_annot = gffpd.read_gff3('genes.gff')
	# print(genes_annot)
	attr_to_columns = genes_annot.attributes_to_columns()
	print(attr_to_columns)
	print(attr_to_columns.source.unique()) # unique = BestRefSeq, Gnomon, BestRefSeq%2CGnomon, tRNAscan-SE, Curated Genomic, RefSeq -> use BestRefSeq
	# df = attr_to_columns[attr_to_columns['source'] == ('Curated Genomic' or 'BestRefSeq')]\
	df = attr_to_columns[attr_to_columns['source'].isin(['BestRefSeq', 'Curated Genomic', 'RefSeq', 'BestRefSeq%2CGnomon'])]
	# df[df['A'].isin([3, 6])]
	# df=attr_to_columns
	print(df)
	chr1 = df[df['seq_id'] == 'NC_000001.11']
	# print(chr1)
	lst = chr1.gene.unique()
	print('CADM3' in lst)
	# print(chr1.gene.unique())
	print(len(df.Name.unique()))
	print(len(df.Name))

if __name__ == '__main__':
	main()