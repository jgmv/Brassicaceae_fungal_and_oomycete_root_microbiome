#!/usr/bin/env bash

### Generates a reference database for identification of Oomycota
###   based on the mitochondrial coxII gene


### prepare environment
# create folders for scripts and data
mkdir -p scripts


### download necessary scripts
# custom bash scripts
if test -f ./scripts/fetch_gb.py
then
  echo "scripts found"
else
  git clone https://github.com/jgmv/ncbi_data_analysis.git scripts
fi

# system-wide access to scripts
export PATH="$PATH:scripts"


### fetch data
# retrieve GIs for all oomycota cox2 from GenBank
#   (or search manually and download)
search_ncbi_gi_by_term.py -o oomycota.gi 'Stramenopiles[Organism] AND ((cytochrome c oxidase subunit II[Title]) OR (coxII[Title])) NOT chromosome NOT genome'

# fecth GB data
fetch_gb.py oomycota.gi -o oomycota.gb

# extract metadata from GB
get_metadata_from_gb.py oomycota.gb -o oomycota.csv

# extract sequences from GB
get_fasta_from_gb.py

# extract taxonomy from GB
sed -n -e '/ORGANISM/,/REFERENCE/ p' oomycota.gb > oomycota.tax
sed -i 's/  ORGANISM  //g' oomycota.tax
sed -i '/REFERENCE/d' oomycota.tax
sed -i ':a;N;$!ba;s/\n            Eukaryota/;Eukaryota/g' oomycota.tax
sed -i ':a;N;$!ba;s/\n            //g' oomycota.tax
sed -i 's/; /;/g' oomycota.tax
sed -i 's/ /_/g' oomycota.tax
sed -i ':a;N;$!ba;s/.\n/\n/g' oomycota.tax


### end
