Use the command "pip3 install -r requirements.txt" to install dependencies

Download files  into data folder

wget -P data/ "$url"

curl -O ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz
curl -O ftp://ftp.ncbi.nih.gov/snp/redesign/latest_release/VCF/GCF_000001405.39.gz

'bcftools view "$ftp" -I -r chr1:159203314-159283574 > /Users/afratzscher/Documents/GitHub/HapPy/results/ACKR1/raw_ACKR1.vcf', 
'bcftools query "$dbSNPFTP" -f "%POS\t%ID\t%REF\t%ALT\n" -r NC_000001.11:159203314-159283574> /Users/afratzscher/Documents/GitHub/HapPy/results/ACKR1/dbSNP_ACKR1.vcf', 
'rm CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz.tbi', 
'rm GCF_000001405.39.gz.tbi']