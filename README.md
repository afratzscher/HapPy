# HapPy

# Table of Contents
  * [Overview](#overview)
  * [Installation](#installation)
    + [Clone Github](#clone-github)
    + [Set up environment](#set-up-environment)
    + [Download files](#download-files)
  * [Running code](#running-code)
    + [Haplotypes](#haplotypes)
    + [Concatenation](#concatenation)
  * [Understanding Output](#understanding-output)
    + [Haplotypes](#haplotypes-1)
    + [Concatenation](#concatenation-1)
  * [Test Run](#test-run)
    + [Haplotypes](#haplotypes-2)
    + [Concatenation](#concatenation-2)
  * [TO DO](#to-do)
    + [Haplotypes](#haplotypes-3)
    + [Concatenation](#concatenation-3)

## Overview
This is package has two objectives:

1. Identifying haplotypes for genes from the 1000 Genomes Project
2. Concatenating haplotypes to construct consensus sequences

Code for `Objective 1` can be found in the [`src/getHaplotypes`](https://github.com/afratzscher/HapPy/tree/master/src/getHaplotypes) folder, while code for `Objective 2` can be found in the [`src/concatenation`](https://github.com/afratzscher/HapPy/tree/master/src/concatenation) folder.


## Installation
### Clone Github
First, clone the github using the following command:

```bash
git clone https://github.com/afratzscher/HapPy
```

### Set up environment
Next, set up the environment. You can do this using conda:

```bash
conda env create --name venv --file=environment.yml
```

*NOTE: make sure to activate the environment prior to running (`conda activate venv`)*

### Download files
Due to some issues with bcftools calling https, please download the following files into the `tmp` folder prior to executing scripts:

1. 1000 Genomes Project data for chromosome `x` (~ 2 GB) (example for chromosome 1):

```bash
curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz
curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz.tbi
```

*NOTE: make sure to download for other chromosomes if you are inputting genes not found on chromosome 1*

2. dbSNP files (~26 GB):

```bash
curl -O ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz
curl -O ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz.tbi
```

## Running code
### Haplotypes
All scripts to get haplotypes can be found in the `/src/getHaplotypes` folder.

To get haplotypes for a specific gene, you can run a command from Terminal:

```bash
python3 main.py -g ACKR1 -d /projects/molonc/roth_lab/afratz/projects/HapPy/tmp/
```

The `-g` flag is used for the name of the gene (in this example, the Duffy gene), and the `d` flag is the location of the `tmp` folder with the previously downloaded dbSNP and 1000 genome project files. 

It is also possible to run for a specific region instead of a gene using the `-r` flag:

```bash
python3 main.py -r chr11:33702782-33725810 -d /projects/molonc/roth_lab/afratz/projects/HapPy/tmp/
```

However, these commands only fetch the haplotypes for a **single gene/region**. To run for multiple genes (e.g. for all genes on the long arm of chromosome 1), you can edit the `names` variable in the `automated_run.py` file to be a list of genes:

```bash
TO EDIT:
data = /projects/molonc/roth_lab/afratz/projects/HapPy/tmp/
names = ['H3-2', 'FAM72C', 'PPIAL4E', 'NBPF15', 'PPIAL4F', 'SRGAP2B']
```

and then run the script using the following command:
```bash
python3 automated_run.py
```

Note: currently all 832 genes on the long arm of chromosome 1 are saved in the `names` variable. If you happen to lose them, you can run the following command and copy the outputed list:
```bash
python3 updateFilesRunAndGetStats.py 
```

### Concatenation
All scripts to run concatenation to get consensus sequences can be found in the `/src/concatenation` folder.

To run, first input the genes to concatenate in the `genes` variable and select a threshold for number of haplotypes to retain at each step:
```bash
TO EDIT:
threshold = 10000
genes = ['FAM72C', 'PPIAL4E', 'NBPF15', 'PPIAL4F']
```

Then run the following command from within the folder:
```bash
python3 combine_pairs.py
```

## Understanding Output

### Haplotypes
Input: 
1. Gene name using -g flag OR region using -r flag
2. Location of 1000 Genomes Project data (tmp folder)

Output (saved in `/results/geneName`)
1. `1000G_geneName.vcf` = raw data from the 1000 genomes project for the gene, outputed by `fetch.py`
2. `dbSNP_geneName.vcf` = raw dbSNP data for the gene, outputed by `fetch.py`
3. `raw_geneName.vcf` = combine dbSNP + 1000GP data, outputed by `fetch.py`
4. `cleaned_geneName.vcf` = data for samples with <= 1 heterozygous SNP, outputed by `cleaner.py`
5. `haplotypes_geneName.vcf = all haplotypes, outputed by `haplotypes.py`
6. `distinct_geneName.vcf` = distinct haplotypes (same as haplotypes_gene.vcf but duplicate haplotypes are removed), outputed by `distinct.py`
7. `mostfreq_geneName.vcf` = sequence of most frequent haplotype + population/superpopulation information, outputed by `popcounts.py`
8. `identical_geneName.vcf` = sequence of distinct haplotype + population/superpopulation information, outputed by `popcounts.py`
9. `full_length_haplotypes_geneName.vcf` = sequences + population/superpopulation information for full-length haplotypes only, outputed by `visualize.py`
10. `visualization/geneName_most_frequent1.png` = pie graph of populations for most frequent haplotype, outputed by `visualize.py`

### Concatenation
Input: 
1. `full_length_haplotypes_geneName.vcf` = full length haplotypes for a gene
OR
2. `distinct_geneName.vcf` = distinct haplotypes for a gene if that gene does not have any full-length haplotypes

Currently, the output is some information printed to the screen. It looks something like this
```bash
0 1 ADDING  ACKR1
11  initial
length:  533
```
It shows the number of long continuous stretch of homozygosity (LCSH, aka number of concatenated pairs) for the first two genes (e.g. 11 in this case), the gene that is trying to be concatenated to the existing sequence (e.g. ACKR1), and the number of LCSHs after adding the gene (e.g. 533). If no overlap is found, you get the following warning: `empty when adding: geneName`.


## Test Run
### Haplotypes
To test if the workflow is running properly, get haplotypes for the `ACKR1` gene. 
Run the following command:
```bash
python3 main.py -g ACKR1 -d /projects/molonc/roth_lab/afratz/projects/HapPy/tmp/
```
OR update the `gene` variable (and `data` variable) in `automated_run.py and run the following command:
```bash
python3 automated_run.py
```

The expected outputed can be found in the [`/results/test`](https://github.com/afratzscher/HapPy/tree/master/results/test) folder

### Concatenation
To test if the workflow is running properly, run the following script:
```bash
python3 test_concatenation.py
```
The expected output is the following
```bash
0 1 ADDING  ACKR1
11  initial
length:  533
1 2 ADDING  FCER1A
length:  201
2 3 ADDING  OR10J1
length:  21
3 4 ADDING  OR10J5
length:  31
4 5 ADDING  APCS
length:  12
12 LCSHs
```

## TO DO
### Haplotypes
1. Consider removing some outputs from getHaplotypes scripts (really only need to keep `full_length_haplotypes_geneName.vcf`, `distinct_geneName.vcf`, and the visualizations) to save space
2. Maybe change format of saved files (json or different type of file instead of .vcf (which is like .csv))?
3. Get haplotypes for all genes on long arm of chromosome 1 (run `automated_run.py`)

### Concatenation
4. Update `combine_pairs.py` so that it does not break when no haplotypes are found (currently, it tries to concatenate the next gene and breaks if there is no overlap, so you have to manually rerun with new set of genes)
5. Increase the threshold for `combine_pairs.py` so that there are fewer breaks
6. Consider using a different ranking OR a different overall approach to decrease number of breaks

