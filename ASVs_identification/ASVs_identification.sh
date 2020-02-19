#!/usr/bin/env bash

### Identifies fungal ITS and oomycete coxII (forward) sequences from
###   MiSeq sequencing against reference databases (UNITE reference 
###   sequences for ITS, in-house database for coxII) using the
###   Naïve Bayesian Classifier (NBC) implemented in Mothur


### prepare environment
# create folders for scripts and data
mkdir -p scripts
mkdir -p ITS_identification
mkdir -p cox_identification


### download necessary scripts
# custom bash scripts
if test -f ./scripts/bash_seq_analysis/otuList.sh
then
  echo "scripts found"
else
  git clone https://github.com/jgmv/bash_seq_analysis.git scripts
fi

for i in $(ls scripts/*.sh)
do
  . $i
done


### ITS sequence identification
# download last version of UNITE ITS reference dataset, mothur release
# check https://unite.ut.ee/repository.php
if test -f ./data/UNITE_sh_dynamic.tax
then
  echo "UNITE database found"
else
  wget -P data https://files.plutof.ut.ee/public/orig/56/25/5625BDC830DC246F5B8C7004220089E032CC33EEF515C76CD0D92F25BDFA9F78.zip
  unzip data/*.zip
  rm data/*.zip
  mv data/UNITE*_sh_dynamic.fasta data/UNITE_sh_dynamic.fasta
  mv data/UNITE*_sh_dynamic.tax data/UNITE_sh_dynamic.tax
  rm data/UNITEv*
fi

# identify ITS sequences using mothur's NBC
mothur "#classify.seqs(fasta=fasta_files/ASV_seqs_its.fasta,\
  template=data/UNITE_sh_dynamic.fasta,\
  taxonomy=data/UNITE_sh_dynamic.tax, cutoff=60, probs=T)" 
mv fasta_files/ASV_seqs_its.UNITE* ITS_identification
mv mothur.* ITS_identification
removeTaxonTag \
  ITS_identification/ASV_seqs_its.UNITE_sh_dynamic.wang.taxonomy \
  ITS_identification/taxonomy_boot_its.csv

# create a copy of the taxonomy file without bootstrap values
cp ITS_identification/taxonomy_boot_its.csv ITS_identification/taxonomy_its.csv
sed -i 's/([^()]*)//g' ITS_identification/taxonomy_its.csv


### cox sequence identification
# identify co sequences using mothur's NBC
mothur "#classify.seqs(fasta=fasta_files/ASV_fwd_seqs_cox.fasta,\
  template=data/oomycota_refs.fasta,\
  taxonomy=data/oomycota_refs.tax, cutoff=60, probs=T)" 
mv fasta_files/ASV_fwd_seqs_cox.oomycota_refs* cox_identification
mv mothur.* cox_identification
removeTaxonTag \
  cox_identification/ASV_fwd_seqs_cox.oomycota_refs.wang.taxonomy \
  cox_identification/taxonomy_boot_cox.csv

# create a copy of the taxonomy file without bootstrap values
cp cox_identification/taxonomy_boot_cox.csv cox_identification/taxonomy_cox.csv
sed -i 's/([^()]*)//g' cox_identification/taxonomy_cox.csv


### end
