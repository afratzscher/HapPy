'''
FILE: settings.py
PURPOSE: parse options given by user (gene name, region, all flag)
INPUT: none
OUTPUT: none
'''
import config
from optparse import OptionParser
import fetch
import pandas as pd

# adding options
def main():
    parser = OptionParser()
    parser.add_option("-g", "--gene", dest="gene", default=False, type='string',
                      help="select by gene", metavar="GENE")    
    parser.add_option("-r", "--region", dest="region", default=False, type='string',
                      help="select by region (GRCh38)", metavar="chr#:#-#")
    parser.add_option("-a", "--all", dest="all", default=False, action="store_true",
                      help="run for all of chr 1")
    parser.add_option("-f", "--folder", dest="foldername", default=False, type='string',
                      help="select a custom folder name (otherwise, is chr#:#-#)")
    (options, args) = parser.parse_args()

    option_dict = vars(options)

    if option_dict.get('all'):
        print('run for all of chromosome 1')
        print("TO IMPLEMENT (settings.py)")

    if option_dict.get('gene') and option_dict.get('region'):
        print("Select either gene or region (not both)")
        exit(0)

    if option_dict.get('foldername'):
        config.__FOLDERNAME__ = option_dict.get('foldername')

    # if given gene, search for region in gene database
    if option_dict.get('gene'):
        config.__GENENAME__ = option_dict.get('gene')
        if config.__FOLDERNAME__ == '':
            config.__FOLDERNAME__ = config.__GENENAME__
        errCode = fetch.main()
        if errCode == -1: # gene not found
            return(-1)

    elif option_dict.get('region'):
        print("TO FIX -> add gene AND whole region as vars...")
        print("TO DO: also fix haplotype (before/after DNE b/c same as gene)")
        if not 'chr' in option_dict.get('region'):
            print("Incorrect format. Use '-h' to get help")
            exit(0)
        string = option_dict.get('region').split(':')
        config.__CHR__ = string[0][3]
        string = string[1].split("-")
        if(string[0]=='' or string[1]==''):
            print("Incorrect format. Use '-h' to get help")
            exit(0)
        config.__START__ = int(string[0])
        config.__END__ = int(string[1])
        config.__GENESTART__ = int(string[0])
        config.__GENEEND__ = int(string[1])

        if config.__FOLDERNAME__ == '':
            config.__FOLDERNAME__ = 'chr'+config.__CHR__+":"+config.__START+"-"+config.__END__

        df = pd.read_csv('GRCh38_chr_versions.txt', sep='\t')
        df = df[df['chr'] == config.__CHR__]
        config.__CHRVERSION__ = df['version'][0]



