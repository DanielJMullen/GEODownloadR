# GEODownloadR
### R package for GEODownloadR to identify, download, and combine .fastq files associated with GSM IDs from the Gene Expression Omnibus (GEO)

##### GEODownloadR contains functions which take sample accession numbers, or GSM IDs (given as GSM#####) and;

1. Identifies their associated SRR accession numbers, returns various information about the sample, links to each sample's associated SRR .fastq files on the European Nucleotide Archive (ENA), and size of the .fastq files.
2. Downloads the .fastq files associated with each SRR accession.
3. Combines the downloaded SRR .fastq files and combines them into .fastq file(s) associated with each of the original GSM IDs.


