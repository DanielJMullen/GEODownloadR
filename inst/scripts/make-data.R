## Code to prepare the `exampleDFGEODownloadR` object

## Successfully created 11/19/2025

## Specify some example GSMs to download for the examples.
exampleGSMVector <- c(
    "GSM4110236",
    "GSM2774971",
    "GSM2406903"
)

## Run the 'identifySRRsFromGSMIDs' first to get the SRR .fastq information.
exampleDFGEODownloadR <- GEODownloadR::identifySRRsFromGSMIDs(
    GSMIDVector = exampleGSMVector,
    pullFastqFileSizes = TRUE
)

## Create the data-raw directory if it doesn't exist.
if(!dir.exists("./data-raw")) {
    dir.create("./data-raw")
}

GSMIDValues = exampleDFGEODownloadR$GSM_IDs
SRRIDValues = exampleDFGEODownloadR$SRR_IDs
SRRSequencingType = exampleDFGEODownloadR$Sequencing_type
ftpSRRFileURLs = exampleDFGEODownloadR$Fastq_URLs
SRRFastqNames = exampleDFGEODownloadR$Fastq_file_names
SRRFastqFileSizes = exampleDFGEODownloadR$Fastq_file_sizes
internalDownloadSRRWithGSMPresent = TRUE
internalLocalDownloadDirectory = "./data-raw"
internalSRRDownloadNodeCount = 1
downloadIndices = 1

## Download the associated SRR .fastqs to the data-raw directory.
exampleDFGEODownloadR <- GEODownloadR::SRRFastqDownload(
    identifySRRsFromGSMIDsDF = exampleDFGEODownloadR,
    localDownloadDirectory = "./data-raw",
    downloadSRRWithGSMPresent = FALSE
)

## Check if the data subdirectory exists, and if it doesn't, create it.
if (!dir.exists("./data")) {
    dir.create("./data")
}

## Save the data.
save(
    exampleDFGEODownloadR,
    file = "./data/exampleDFGEODownloadR.rda",
    compress = "xz",
    compression_level = -9 ## Level -n = level n but with xz -e
)
