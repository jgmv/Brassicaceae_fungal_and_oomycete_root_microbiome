# Brassicaceae_fungal_and_oomycete_root_microbiome
Data and code for the process of sequence data and data analyses used in [Maciá-Vicente, Nam, Thines (2020) Root filtering, rather than host identity or age, determines the composition of root-associated fungi and oomycetes in three naturally co-occurring Brassicaceae. Soil Biology & Biochemistry, 146:107806](https://doi.org/10.1016/j.soilbio.2020.107806).

## Contents
### oomycota_taxonomy_database
Code to generate a database of *coxII* gene sequence data from NCBI GenBank, for its use as a reference in the identification of oomycete sequences. The file `oomycota_20190626.gi` provided contains GenBank GI numbers for the sequence records used in the paper (fetched on 2019-06-26).

### ASVs_identification
Data and code to identify MiSeq sequence data using the Naïve Bayesian Classifier tool [(Wang et al. 2007)](https://doi.org/10.1128/AEM.00062-07) implemented in mothur [(Schloss et al. 2009)](https://doi.org/10.1128/AEM.01541-09). Identifications are based on comparisons of fungal ITS sequences against records in the [UNITE](https://unite.ut.ee/) database of reference ITS sequences [(Kõljalg et al. 2005)](https://doi.org/10.1111/j.1469-8137.2005.01376.x), and of oomycete *coxII* sequences against the inhouse database described [above](https://github.com/jgmv/Brassicaceae_fungal_and_oomycete_root_microbiome/tree/master/oomycota_taxonomy_database).

### data_analysis
Data and code for the statistical analyses in R. The following R packages are required: `DESeq`,`emdbook`,`gplots`,`Hmisc`,`MASS`,`RColorBrewer`, and `vegan`.

## Additional files
* MiSeq sequence data is available at the NCBI Sequence Read Archive under BioProject number [PRJNA593383](https://www.ncbi.nlm.nih.gov/bioproject/593383).
* The code for processing raw MiSeq sequence data is available at [https://github.com/jgmv/MiSeq_process](https://github.com/jgmv/MiSeq_process).
