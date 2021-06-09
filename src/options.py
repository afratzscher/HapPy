'''
FILE: options.py
PURPOSE: parse options given by user (gene name, region, all flag)
INPUT: none
OUTPUT: none
'''
from optparse import OptionParser

# adding options
def main():
    parser = OptionParser()
    parser.add_option("-g", "--gene", dest="gene", default=False, type='string',
                      help="select by gene", metavar="GENE")    
    parser.add_option("-r", "--region", dest="region", default=False, type='string',
                      help="select by region (GRCh38)", metavar="chr#:#-#")
    parser.add_option("-a", "--all", dest="all", default=False, action="store_true",
                      help="run for all of chr 1")
    parser.add_option("-c", "--custom", dest="foldername", default=False, type='string',
                      help="select a custom folder name (otherwise, is chr#:#-#)")
    parser.add_option("-f", "--filter", dest="filter", default=False, action="store_true",
                      help="filters data to only include original 2504 individuals")
    
    (options, args) = parser.parse_args()

    return options



