# combines results for 11 genes into single DF
import os

genes = ['ACKR1', 'AIM2', 'APCS', 'CADM3',
		'FCER1A', 'IFI16', 'OR10J1', 'CRP', 
		'OR10J5', 'PYDC5', 'PYHIN1']
direct = '/'.join(os.getcwd().split('/')[:-1]) + "/results/"
# file = direct +  "/data.json"

for gene in genes:
	file = direct + gene + "/"
	print(file)