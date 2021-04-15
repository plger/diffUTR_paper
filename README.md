# Streamlining differential exon and 3' UTR analysis

This repository accompanies the [diffUTR paper](https://doi.org/10.1186/s12859-021-04114-7), providing a record of what was done and of the data used. For information about the diffUTR package and its usage, see the [package's repository](https://github.com/ETHZ-INS/diffUTR).

## Structure

* the [results](results/) folder contains the results of each method, applied to each dataset.
* the [data](data/) folder contains some necessary datasets (see about large data files, below).
* the [figures](figures/) folder contains the code used to generate the figures.
* the root directory contains files directly used in the manuscript.

## Large data files

Some files were two large to be included here. These include:

* The gencode gene annotation, which should be decompressed and placed in `data/large/`. It can be downloaded from `ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gtf.gz`.
* The raw data and aligned reads for the different datasets used. The reads from the published LTP dataset can be downloaded from the GEO/SRA ([GSE84643](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84643)). The simulated reads can be downloaded from [figshare](https://dx.doi.org/10.6084/m9.figshare.13726143).

