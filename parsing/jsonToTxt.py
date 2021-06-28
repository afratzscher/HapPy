import pandas as pd
import os
from datetime import date, time 

today = date.today()

def main():
	today = date.today().strftime("%B %d, %Y")
	direct = '/'.join(os.getcwd().split('/')[:-1]) + "/src"
	file = direct +  "/data.json"
	txtfile = open(direct + "/data.txt", "w")
	txtfile.write("# Updated on: " + today + "\n")
	df = pd.read_json(file)
	df.to_csv(txtfile, index = False, mode = 'a')

if __name__ == '__main__':
	main()