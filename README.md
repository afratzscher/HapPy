# HapPy

## Instructions
* You will need to download package dependencies. You can do so using conda:
     
     `conda env create --name venv --file=environment.yml`
     * NOTE: make sure to activate the environment prior to running (`conda activate venv`)

* Additionally, you will need to download the following large files to fetch gene locations using the following commands:

     1. 1000 Genomes Project data for chromosome `x` (~ 2 GB) (example for chromosome 1):
          
          `curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz`
     
     2. dbSNP files (~26 GB):
    
          `curl -O ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz`
     
          `curl -O ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz.tbi`
    
