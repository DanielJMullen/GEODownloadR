# GEODownloadR
### R package for GEODownloadR to identify, download, and combine .fastq files associated with GSM IDs from the Gene Expression Omnibus (GEO)

##### GEODownloadR contains the following functions which take sample accession numbers, or GSM IDs (given as GSM#####) and:

1. identifySRRsFromGSMIDs: Identifies their associated SRR accession numbers, returns various information about the sample, links to each sample's associated SRR .fastq files on the European Nucleotide Archive (ENA), and size of the .fastq files ().

2. SRRFastqDownload: Downloads the .fastq files associated with each SRR accession.

3. combineSRRFastqsToGSMs: Combines the downloaded SRR .fastq files and combines them into .fastq file(s) associated with each of the original GSM IDs.

These functions are designed to be run in sequence, with each returning a data frame of information that is taken as input by the next function and then returned with additional info. The SRRFastqDownload function will also download SRR .fastq files to a user-specified directory, and the combineSRRFastqsToGSMs function will take the SRR .fastqs in this directory and combine them into ones associated with each GSM ID given by the user.

### Installing GEODownloadR

Currently, GEODownloadR can be installed using tools which allow users to install an R package from a github repository, such as the install_github function in the devtools package.

```r
## Install devtools if required.
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

## Use install_github to install the GEODownloadR package.
devtools::install_github("DanielJMullen/GEODownloadR")

## Load the library
library(GEODownloadR)
```
