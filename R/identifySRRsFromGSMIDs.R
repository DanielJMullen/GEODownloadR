## Internal functions for the identifySRRsFromGSMIDs function:

## An internal function which takes individual GSM ID values and returns
## information on the GSE page and SRX page for them.
.GSEInformationFunction <- function(GSMIDValues) {

    ## Initialize a holding value for the GSMIndices, which will be used shortly
    ## when setting up parallelization using foreach().
    GSMIndices <- NULL

    ## Create a list of the nodes that will be used, and then register them for
    ## use with the foreach package. 60 nodes will be used by default.
    makeGSEClusterList <- parallel::makeCluster(60)
    doParallel::registerDoParallel(makeGSEClusterList)

    ## Use foreach() to iterate across the GSM ID values, and extract the
    ## relevant information from each, then combine this information into a
    ## dataframe using rbind. The relevant information here includes ...
    GSMInfoDF <- foreach::foreach(
        GSMIndices = seq_along(GSMIDValues),
        .combine = rbind
    ) %dopar% {

        ## Get the individual GSMID to use for this iteration, akin to a
        ## for loop.
        iterGSMIDValue <- GSMIDValues[GSMIndices]

        ## Set a URL to the GSM's page on the NCBI site.
        GSMPageURL <- paste0(
            "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",
            iterGSMIDValue
        )

        ## Set up an initial fail count - the function is designed to try
        ## accessing information 3 times. Also, create an empty GSMPage to load
        ## text into - if this function fails, we can also check if this
        ## variable remains NULL to see if information from the page was
        ## successfully loaded.
        GSEInfoFailCount <- 0
        GSMPage <- NULL

        while (GSEInfoFailCount < 3 && is.null(GSMPage)) {

            ## Establish a connection to the specified GSM's page.
            GSMPageURLConnection <- url(GSMPageURL)

            ## Try reading the webpage for the given GSM. If this fails, sleep
            ## 3 seconds before proceeding.
            tryCatch({
                GSMPage <- readLines(GSMPageURLConnection)
            }, error = function(cond) {
                Sys.sleep(3)
            })

            ## Close the connection to the GSM's page.
            close(GSMPageURLConnection)

            ## Increment the GSEInfoFailCount.
            GSEInfoFailCount <- (GSEInfoFailCount + 1)
        }

        ## There are two potential fail cases to account for. First, where a
        ## page is never loaded, and second where the GSM value is not in the
        ## database and a generic page is loaded instead (which will not have
        ## relevant information).
        if (is.null(GSMPage)) {
            return(
                c(
                    iterGSMIDValue,
                    GSMPageURL,
                    rep(NA, 4),
                    "FAIL"
                )
            )
        } else if (length(grep('GSE', GSMPage)) == 0) {
            return(
                c(
                    iterGSMIDValue,
                    GSMPageURL,
                    rep(NA, 4),
                    "FAIL"
                )
            )
        }

        ## Identify the GSE(s) the file is associated with:

        ## Get the strings on the webpage which contain "GSE".
        rawGSEStrings <- GSMPage[grepl("GSE", GSMPage, fixed=TRUE)]

        ## strsplit the rawGSEStrings based on the term "GSE".
        rawGSEStringsStrsplit <- strsplit(
            rawGSEStrings,
            "GSE"
        )

        ## The listing of each GSE associated with the GSM is listed at
        ## the start of the second elements of each line after splitting
        ## on "GSE".
        rawGSEStringsStrsplitGSEElements <- sapply(
            rawGSEStringsStrsplit,
            "[[",
            2
        )

        ## Get the parts of the string before the annotation '\"
        ## onmouseout' "GSE" also needs to be readded to the start of
        ## the values. All GSE values will also be collapsed into a
        ## single string with the GSEs delimited by commas.
        GSEValueString <- paste0(
            "GSE",
            sub(
                '\" onmouseout.*',
                '',
                rawGSEStringsStrsplitGSEElements
            ),
            collapse = ","
        )

        ## Occasionally there is an additional '\">"' at the end of the
        ## GSE string. If it is present, remove it.
        GSEValueString <- sub(
            '\".*',
            '',
            GSEValueString
        )

        ## Get the webpages associated with the given GSE entries.
        ## Again, these will also be collapsed into a single string with
        ## the URLs delimited by commas.
        GSEPageURL <- paste0(
            'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
            unlist(strsplit(GSEValueString, ",")),
            collapse = ","
        )

        ## This gets the string which contains the link to the SRA page
        ## from the overall GSM's webpage.
        rawSRXString <- GSMPage[grepl("term=SRX", GSMPage, fixed=TRUE)]

        ## Extract the SRX URL from the raw SRX string.
        ## The outer sub is needed to remove an escape character "\"
        ## left at the start of the link after the second level sub
        ## statement, as it is hard to quote properly.
        SRXURL <- sub(
            ".",
            "",
            sub(
                ".*href=",
                "",
                sub(
                    "\">.*",
                    "",
                    rawSRXString
                )
            )
        )

        ## Get the SRXID from the SRXURL. This is the part at the end of the
        ## URL string after 'term=SRX'. This is to better ensure a single match
        ## so we will readd the SRX back in after the removal.
        SRXID <- paste0(
            "SRX",
            sub(
                ".*term=SRX",
                "",
                SRXURL
            )
        )

        ## If successful, return the GSMPageURL, GSEValueString,
        ## GSEPageURL, and SRXURL. Additionally, return a seventh element
        ## which we set to "PASS" for now, so we know this step worked.
        return(
            c(
                iterGSMIDValue,
                GSMPageURL,
                GSEValueString,
                GSEPageURL,
                SRXID,
                SRXURL,
                "PASS"
            )
        )
    }

    ## Close the cluster list.
    ## suppressWarnings is used to silence a warning that relates to closing
    ## unused connections, which doesn't affect how the function works.
    suppressWarnings(parallel::stopCluster(makeGSEClusterList))

    ## Set the rownames of the dataframe to be output as the GSMIDValues.
    row.names(GSMInfoDF) <- GSMIDValues

    ## Set the colnames of the dataframe to describe the data contained in each
    ## column:
    colnames(GSMInfoDF) <- c(
        "GSM_IDs",
        "GSM_ID_URLs",
        "GSE_IDs",
        "GSE_ID_URLs",
        "SRX_IDs",
        "SRX_ID_URLs",
        "Passed_GSE_information"
    )

    ## Return the GSMInfoDF.
    return(GSMInfoDF)
}

## An internal function which takes the SRX page URLs derived for the original
## GSM values input by the user, and returns the SRR values associated with the
## GSM, links to those SRR files in the European Nucleotide Archive (ENA), and
## the type of sequencing that was done (single-end or paired).
.SRXInformationFunction <- function(GSMIDValues, SRXIDURLs) {

    ## Initialize a holding value for the SRXIDURLs, which will be used
    ## shortly when setting up parallelization using foreach().
    SRXURLIndices <- NULL

    ## Create a list of the nodes that will be used, and then register them for
    ## use with the foreach package. 60 nodes will be used by default.
    makeSRXClusterList <- parallel::makeCluster(60)
    doParallel::registerDoParallel(makeSRXClusterList)

    ## Use foreach() to iterate across the SRX URLs, and extract the
    ## relevant information from each, then combine this information into a
    ## dataframe using rbind. The relevant information here includes ...
    SRXInformationDF <- foreach::foreach(
        SRXURLIndices = seq_along(SRXIDURLs),
        .combine = rbind
    ) %dopar% {

        ## Get the associated GSM value, akin to a for loop.
        iterGSMIDValue <- GSMIDValues[SRXURLIndices]

        ## Get the individual SRX URL to use for this iteration, akin to a
        ## for loop.
        iterSRXURLValue <- SRXIDURLs[SRXURLIndices]

        ## Next, we need to check if the iterSRXURLValue is NA. If so, we will
        ## need to cancel out of this step to prevent an error from being
        ## induced as we can't load an NA URL.
        if(is.na(iterSRXURLValue)) {
            return(
                c(
                    iterGSMIDValue,
                    iterSRXURLValue,
                    rep(NA, 10),
                    "FAIL"
                )
            )
        }

        ## Set up an initial fail count - the function is designed to try
        ## accessing information 3 times. Also, create an empty SRXPage to load
        ## text into - if this function fails, we can also check if this
        ## variable remains NULL to see if information from the page was
        ## successfully loaded
        SRXInfoFailCount <- 0
        SRXPage <- NULL

        while (SRXInfoFailCount < 3 && is.null(SRXPage)) {

            ## Establish a connection to the specified SRX's page.
            SRXPageURLConnection <- url(iterSRXURLValue)

            ## Try reading the webpage for the given SRX. If this fails, sleep
            ## 3 seconds before proceeding. NOTE: supressWarnings is included to
            ## suppress an irrelevant warning message about an incomplete final
            ## line when reading the SRX page URL.
            tryCatch({
                SRXPage <- suppressWarnings(readLines(SRXPageURLConnection))
            }, error = function(cond) {
                Sys.sleep(3)
            })

            ## Close the connection to the SRX's page.
            close(SRXPageURLConnection)

            ## Increment the SRRInfoFailCount.
            SRXInfoFailCount <- (SRXInfoFailCount + 1)
        }

        ## There are two potential fail cases to account for. First, where a
        ## page is never loaded, and second where the SRX value is not in the
        ## database (this I imagine is hard to come across since we pull the
        ## value from a valid GSM ID, but just in case I've included the
        ## check here) and a generic page is loaded instead (which will not have
        ## relevant information).
        if (is.null(SRXPage)) {
            return(
                c(
                    iterGSMIDValue,
                    iterSRXURLValue,
                    rep(NA, 10),
                    "FAIL"
                )
            )
        } else if (length(grep('SRR', SRXPage)) == 0) {
            return(
                c(
                    iterGSMIDValue,
                    iterSRXURLValue,
                    rep(NA, 10),
                    "FAIL"
                )
            )
        }

        ## This gets the string which contains the link to the SRR page by
        ## pulling the string which contains "SRR".
        rawSRRString <- SRXPage[grepl("SRR", SRXPage, fixed=TRUE)]

        ## There may be multiple SRRs listed, so we will split the string
        ## containing the SRR IDs in a way that leaves the IDs
        ## themselves at the start of each element in the list, except the first
        ## one (which will contain no SRR IDs).
        rawSRRStringVector <- unlist(
            strsplit(rawSRRString, "Traces?run=", fixed = TRUE)
        )

        ## Remove the first element of the vector (it has no IDs in it).
        rawSRRStringVectorDropFirst <- rawSRRStringVector[-1]

        ## Get the SRR ID number(s) by getting the part of the string
        ## with the SRR IDs before ' "\"> '.
        SRRIDs <- sub("\">.*", "", rawSRRStringVectorDropFirst)

        ## Create links to the European Nucleotide Archive for each SRR.
        ENASRRURLs <- paste0(
            "http://www.ebi.ac.uk/ena/data/view/",
            SRRIDs
        )

        ## Get the information about the number of spots (reads).

        ## First, split each of the rawSRRStringVectorDropFirst strings on
        ## 'align=\"right\">', as the data we need for each SRR is not well
        ## annotated in the raw text.
        rawSRRStringVectorDropFirstAndSplit <- strsplit(
            rawSRRStringVectorDropFirst,
            "align=\"right\">",
            fixed = TRUE
        )

        ## The read counts in each SRR can be extracted from the second element
        ## of each SRR's list by getting the part of the string before
        ## '</td><td ', then the resulting string can have the commas extracted
        ## from it with gsub() and can be converted to a numeric value.
        ENASRRReadCounts <- as.numeric(
            gsub(
                ",",
                "",
                sub(
                    "</td><td .*",
                    "",
                    sapply(rawSRRStringVectorDropFirstAndSplit, "[[", 2)
                )
            )
        )

        ## The approximate total bases in each SRR can be extracted from the
        ## third element of each SRR's list by getting the part of the string
        ## before '</td><td '. This count is approximate, with a letter used for
        ## the approximate metric value, so we will keep that.
        ENASRRTotalBases <- sub(
            "</td><td .*",
            "",
            sapply(rawSRRStringVectorDropFirstAndSplit, "[[", 3)
        )

        ## The approximate file size of each SRR can be extracted from the
        ## fourth element of each SRR's list by getting the part of the string
        ## before '</td><td>'. This value is approximate, with letters used for
        ## the byte value, so we will keep that.
        ENASRRFileSize <- sub(
            "</td><td>.*",
            "",
            sapply(rawSRRStringVectorDropFirstAndSplit, "[[", 4)
        )

        ## The date each SRR was uploaded can also be extracted from the
        ## fourth element of each SRR's list by getting the part of the string
        ## after '</td><td>', but before '</td></tr></tbody> and sometimes
        ## '</td></tr><tr>'.
        ENASRRFileDate <- sub(
            "</td></tr><tr>.*",
            "",
            sub(
                "</td></tr></tbody>.*",
                "",
                sub(
                    ".*</td><td>",
                    "",
                    sapply(rawSRRStringVectorDropFirstAndSplit, "[[", 4)
                )
            )
        )

        ## Let's also get the organism studied, instrument used, as well as the
        ## experimental Strategy and Layout (sequencing format):

        ## Extract the organism the sequencing data was derived from.

        ## First, extract the part of the rawSRRString before
        ## '</a></span></div></div><div class=\"expand showed sra-full-data\">
        ## Library:' and after 'Organism: <span><a href=\"'.
        sequencingOrganism <- sub(
            '.*Organism: <span><a href=\"',
            '',
            sub(
                '</a></span></div></div><div class=\"expand showed sra-full-data\">Library:.*',
                '',
                rawSRRString
            )
        )

        ## Then, get the part of the resulting string after '' to remove the
        ## part of the string with the link to the organism's ID on the
        ## Taxonomy Browser, which is not needed.
        sequencingOrganism <- sub(
            '.*\">',
            '',
            sequencingOrganism
        )

        ## Extract the instrument used for the sequencing. To do this, extract
        ## the part of the rawSRRString before '</span></div><div>Strategy:' and
        ## after 'Instrument: <span>'.
        sequencingInstrument <- sub(
            '.*Instrument: <span>',
            '',
            sub(
                '</span></div><div>Strategy:.*',
                '',
                rawSRRString
            )
        )

        ## Extract the experiment strategy. To do this, extract the part of the
        ## rawSRRString before '</span></div><div>Source:' and after
        ## 'Strategy: <span>'.
        sequencingStrategy <- sub(
            '.*Strategy: <span>',
            '',
            sub(
                '</span></div><div>Source:.*',
                '',
                rawSRRString
            )
        )

        ## Also extract the type of sequencing which is also contained in the
        ## rawSRRString as well. To do this, extract the part of the
        ## rawSRRString before '</span></div><div>Construction protocol:' and
        ## after 'Layout: <span>'.
        sequencingFormat <- sub(
            '.*Layout: <span>',
            '',
            sub(
                '</span></div><div>Construction protocol:.*',
                '',
                rawSRRString
            )
        )

        ## For some page layout, an additional long string is left trailing. We
        ## can curtail it by keeping the start of the string before
        ## '</span></div></div>' (this won't affect anything when the string
        ## isn't present).
        sequencingFormat <- sub(
            '</span></div></div>.*',
            '',
            sequencingFormat
        )

        ## If successful, return the SRRIDs (comma-delimited), ENASRRURLs, the
        ## data on the SRR files, and the overall sequencing data. Additionally,
        ## return a eleventh element which we set to "PASS" for now, so we know
        ## this step worked. If multiple entries are noted for the dataset,
        ## this will be set up as a data.frame:
        return(
            cbind(
                rep(iterGSMIDValue, length(SRRIDs)),
                rep(iterSRXURLValue, length(SRRIDs)),
                SRRIDs,
                ENASRRURLs,
                ENASRRReadCounts,
                ENASRRTotalBases,
                ENASRRFileSize,
                ENASRRFileDate,
                rep(sequencingOrganism, length(SRRIDs)),
                rep(sequencingInstrument, length(SRRIDs)),
                rep(sequencingStrategy, length(SRRIDs)),
                rep(sequencingFormat, length(SRRIDs)),
                rep("PASS", length(SRRIDs))
            )
        )
    }

    ## Close the cluster list.
    ## suppressWarnings is used to silence a warning that relates to closing
    ## unused connections, which doesn't affect how the function works.
    suppressWarnings(parallel::stopCluster(makeSRXClusterList))

    ## Set the rownames of the dataframe to be output as the SRR values.
    row.names(SRXInformationDF) <- SRXInformationDF[,3]

    ## Set the colnames of the dataframe to describe the data contained in each
    ## column:
    colnames(SRXInformationDF) <- c(
        "GSM_IDs",
        "SRX_ID_URLs",
        "SRR_IDs",
        "SRR_ID_URLs",
        "Total_reads",
        "Approximate_total_bases",
        "Approximate_file_size",
        "Upload_date",
        "Species_name",
        "Sequencing_instrument",
        "Experiment_type",
        "Sequencing_type",
        "Passed_SRX_information"
    )

    ## Return the SRXInformationDF.
    return(SRXInformationDF)
}

## An internal function which takes the SRR page URLs derived for the original
## GSM values input by the user, and returns the individual SRR files and their
## ftp links for downloading.
.SRRFileIdentificationFunction <- function(
    GSMIDValues,
    SRRIDValues,
    SRRSequencingType
) {

    ## Initialize a holding value for the SRRIndices, which will be used shortly
    ## when setting up parallelization using foreach().
    SRRIndices <- NULL

    ## Create a list of the nodes that will be used, and then register them for
    ## use with the foreach package. 60 nodes will be used by default.
    makeSRRIDClusterList <- parallel::makeCluster(60)
    doParallel::registerDoParallel(makeSRRIDClusterList)

    ## Use foreach() to iterate across the SRR URLs, and extract the
    ## relevant information from each, then combine this information into a
    ## dataframe using rbind. The relevant information here includes ...
    SRRFileIDDF <- foreach::foreach(
        SRRIndices = seq_along(SRRIDValues),
        .combine = rbind
    ) %dopar% {

        ## Get the associated GSM ID value, akin to a for loop.
        iterGSMIDValue <- GSMIDValues[SRRIndices]

        ## Get the individual SRR ID to use for this iteration, akin to a
        ## for loop.
        iterSRRIDValue <- SRRIDValues[SRRIndices]

        ## Get the type of sequencing (single or paired end) that was done for
        ## the given SRR.
        iterSRRSequencingTypeValue <- SRRSequencingType[SRRIndices]

        ## Next, we need to check if the iterSRRURLValue is NA. If so, we will
        ## need to cancel out of this step to prevent an error from being
        ## induced as we can't load an NA URL.
        if(is.na(iterSRRIDValue)) {
            return(
                c(
                    iterGSMIDValue,
                    iterSRRIDValue,
                    rep(NA, 2),
                    "FAIL"
                )
            )
        }

        ## Set up an initial fail count - the function is designed to try
        ## accessing information 3 times. Also, create an empty SRRManifest to
        ## load text into - if this function fails, we can also check if this
        ## variable remains NULL to see if information from the page was
        ## successfully loaded
        SRRInfoFailCount <- 0
        SRRManifest <- NULL

        while (SRRInfoFailCount < 3 && is.null(SRRManifest)) {

            ## Try reading the webpage with info on SRRs.
            ## If this fails, sleep 3 seconds before proceeding.
            tryCatch({
                SRRManifest <- read.table(
                    paste0(
                        "https://www.ebi.ac.uk/ena/portal/api/filereport?",
                        "accession=",
                        iterSRRIDValue,
                        "&result=read_run&fields=experiment_accession,",
                        "run_accession,fastq_ftp&format=tsv&download=",
                        "true&limit=0"
                    ),
                    sep ="\t",
                    header = TRUE
                )
            }, error = function(cond) {
                Sys.sleep(3)
            })

            ## Increment the SRRInfoFailCount.
            SRRInfoFailCount <- (SRRInfoFailCount + 1)
        }

        ## There are three potential fail cases to account for. First, where a
        ## page is never loaded, the second where the 'fastq_ftp' column is not
        ## in the loaded dataframe, and the last where the given SRR value
        ## cannot be found in the ftp links of that column's data.
        ## (though I imagine these cases are hard to come across since we pull
        ## the value from a valid SRX ID, but just in case I've included the
        ## check here) and a generic page is loaded instead (which will not have
        ## relevant information).
        if (is.null(SRRManifest)) {
            return(
                c(
                    iterGSMIDValue,
                    iterSRRIDValue,
                    rep(NA, 2),
                    "FAIL"
                )
            )
        } else if (!("fastq_ftp" %in% colnames(SRRManifest))) {
            return(
                c(
                    iterGSMIDValue,
                    iterSRRIDValue,
                    rep(NA, 2),
                    "FAIL"
                )
            )
        } else if (length(grep(iterSRRIDValue, SRRManifest$fastq_ftp)) == 0) {
            return(
                c(
                    iterGSMIDValue,
                    iterSRRIDValue,
                    rep(NA, 2),
                    "FAIL"
                )
            )
        }

        ## The ftp links can be acquired by strsplitting the values in the
        ## 'fastq_ftp' column of the loaded SRRManifest, then pasting 'ftp://'
        ## to the start of them.
        ftpFastqSRRURLValues <- paste0(
            "ftp://",
            unlist(strsplit(SRRManifest$fastq_ftp, ";"))
        )

        ## Get the name of the SRR fastq.gz file that will be downloaded.
        SRRFastqFileNames <- basename(ftpFastqSRRURLValues)

        ## If successful, return the GSM IDs, SRR IDs, the URL to each
        ## constituent SRR file's ftp file, the name of the base fastq(s),
        ## and note whether or not each file succeeded.
        return(
            cbind(
                rep(iterGSMIDValue, length(ftpFastqSRRURLValues)),
                rep(iterSRRIDValue, length(ftpFastqSRRURLValues)),
                ftpFastqSRRURLValues,
                SRRFastqFileNames,
                rep("PASS", length(ftpFastqSRRURLValues))
            )
        )
    }

    ## Close the cluster list.
    ## suppressWarnings is used to silence a warning that relates to closing
    ## unused connections, which doesn't affect how the function works.
    suppressWarnings(parallel::stopCluster(makeSRRIDClusterList))

    ## Set the colnames of the dataframe to describe the data contained in each
    ## column:
    colnames(SRRFileIDDF) <- c(
        "GSM_IDs",
        "SRR_IDs",
        "Fastq_URLs",
        "Fastq_file_names",
        "Passed_SRR_file_identification"
    )

    ## Return the SSRRFileIDDF.
    return(SRRFileIDDF)
}

## An internal function which takes the ftp paths to files, pulls relevant
## information on them, particularly the individual file size.
.SRRInformationFunction <- function(
    GSMIDValues,
    SRRIDValues,
    SRRSequencingType,
    ftpSRRFileURLs,
    SRRFastqNames,
    errorWaitTime
) {
    ## Check if the ftpSRRFileURLs is NA. If so, we will need to cancel out of
    ## this step to prevent an error from being induced as we can't load an NA
    ## URL.
    if(is.na(ftpSRRFileURLs)) {
        return(
            c(
                GSMIDValues,
                SRRIDValues,
                SRRSequencingType,
                ftpSRRFileURLs,
                SRRFastqNames,
                rep(NA, 2),
                "FAIL"
            )
        )
    }

    ## Now load the text of the ftp pages. This will be used to extract the
    ## exact size of the files.

    ## If getSRRQuick is FALSE, don't try to load the page, instead
    ## proceed directly to the download section without checking the file size.

    ## Set up a fail count - the function is designed to try accessing
    ## information 3 times. Also, create an empty FTPSRRPage to load text
    ## into - if this function fails, we can also check if this variable
    ## remains NULL to see if information from the page was successfully
    ## loaded.
    FTPSRRInfoFailCount <- 0
    FTPSRRPage <- NULL

    while (FTPSRRInfoFailCount < 3 && is.null(FTPSRRPage)) {

        ## Try reading the data for the given SRR using the getURL()
        ## function from the RCurl package. To set this up, we need to
        ## strip the fastq file from the ftpSRRURLValue, and add back
        ## the trailing "/".If this fails, sleep 5 seconds before
        ## proceeding.
        tryCatch({
            FTPSRRPage <- RCurl::getURL(
                paste0(dirname(ftpSRRFileURLs),"/")
            )
        }, error = function(cond) {
            Sys.sleep(errorWaitTime)
        })

        ## Increment the ftpSRRInfoFailCount.
        FTPSRRInfoFailCount <- (FTPSRRInfoFailCount + 1)
    }

    ## There are two potential fail cases to account for. First, where a
    ## page is never loaded, and second where the SRR value is not in
    ## the database (this I imagine is hard to come across since we pull
    ## the value from a valid GSM/SRX ID, but just in case I've included
    ## the check here) and a generic page is loaded instead (which will
    ## not have relevant information).
    if (is.null(FTPSRRPage)) {
        return(
            c(
                GSMIDValues,
                SRRIDValues,
                SRRSequencingType,
                ftpSRRFileURLs,
                SRRFastqNames,
                rep(NA, 2),
                "FAIL"
            )
        )
    } else if (length(grep(SRRFastqNames, FTPSRRPage)) == 0) {
        return(
            c(
                GSMIDValues,
                SRRIDValues,
                SRRSequencingType,
                ftpSRRFileURLs,
                SRRFastqNames,
                rep(NA, 2),
                "FAIL"
            )
        )
    }

    ## Because two (or rarely, more) files may be listed under a single SRR ftp
    ## page if the experiment happens to be paired, we will need to determine
    ## which entry we are looking at.
    fileN <- grep(
        SRRFastqNames,
        unlist(
            strsplit(
                FTPSRRPage,
                "\n"
            )
        )
    )

    ## The page text is loaded as one big string, so we will need to
    ## strsplit it in order to extract the parts of interest.
    FTPSRRPageStrsplit <- strsplit(
        unlist(strsplit(FTPSRRPage, '\n')), ' +'
    )

    ## The size of the given file is in the 5th element of each list.
    SRRFastqFileSize <- as.numeric(
        vapply(FTPSRRPageStrsplit, '[', '', 5)[[fileN]]
    )

    ## The dates of first upload can be created by pasting elements
    ## 6-8.
    SRRFastqFileFirstDate <- paste(
        vapply(FTPSRRPageStrsplit, '[', '', 6)[[fileN]],
        vapply(FTPSRRPageStrsplit, '[', '', 7)[[fileN]],
        vapply(FTPSRRPageStrsplit, '[', '', 8)[[fileN]],
        sep = "-"
    )

    ## Do the return:
    return(
        c(
            GSMIDValues,
            SRRIDValues,
            SRRSequencingType,
            ftpSRRFileURLs,
            SRRFastqNames,
            SRRFastqFileSize,
            SRRFastqFileFirstDate,
            "PASS"
        )
    )
}

## Identify SRR files associated with GEO GSM IDs:

#' This function takes GSM IDs provided by the user and identifies the SRRs
#' associated with them. This function also provides URLs to their corresponding
#' .fastq files hosted by the European Nucleotide Archive (ENA). If selected,
#' this function will also query the ENA server to identify the size of these
#' .fastq files, which will take awhile to run but is essential for this
#' package's download pipeline. This function is the first step in the pipeline
#' to download and combine .fastq files associated with the user's GSM IDs of
#' interest.
#'
#' @param GSMIDVector Specify a vector of sample accession numbers to
#' investigate/download files for. The values in the GSMIDVector should be given
#' in the form "GSM" followed by the numbers associated with the sample, such as
#' "GSM2437766".
#' @param pullFastqFileSizes Set to TRUE to pull the size of the fastq file(s)
#' associated with each SRR, otherwise set to FALSE to not. This option will
#' take awhile to run but is essential to successfully download and combine the
#' fastq files later. Defaults to TRUE.
#' @return Returns a dataframe containing a dataframe with information about the
#' SRR files associated with the GSM IDs specified in the `GSMIDVector`,
#' including paths to their associated .fastq files on the ENA server.
#' @export
#'
#' @examplesIf interactive()
#' ## This example requires an internet connection, and will attempt to pull
#' ## information, including .fastq file size, for three example GSM IDs and
#' ## returns a data frame with that information.
#' returnValue <- identifySRRsFromGSMIDs(
#'     GSMIDVector = c("GSM4110236", "GSM2774971", "GSM2406903")
#' )

identifySRRsFromGSMIDs <- function(
    GSMIDVector,
    pullFastqFileSizes = TRUE
) {

    ## Verify the GSMIDVector argument.

    ## First ensure the GSMIDVector is a character vector.
    if (!.isCharacterVector(GSMIDVector)) {
        stop(
            "The `GSMIDVector` argument must be a character vector ",
            "consisting of GSM IDs each with 'GSM' followed by numerical ",
            "values."
        )
    }

    ## Then verify that the GSMIDs are valid.
    .isValidGSMIDs(GSMIDVector)

    ## Verify the pullFastqFileSizes argument.
    if (!.isSingularBoolean(pullFastqFileSizes)) {
        stop(
            "The `pullFastqFileSizes` argument must be a singular boolean ",
            "value."
        )
    }

    ## Run the .GSEInfoFunction() to pull infomration on the GSE and SRX pages
    ## associated with each GSMID.
    GSEInformationDF <- as.data.frame(.GSEInformationFunction(GSMIDVector))

    ## Do a second pass at downloading by taking the GSMs which failed the
    ## first time and re-running those through the .GSEInformationFunction.
    failedGSEInformationDFIters <- unique(
        which(GSEInformationDF$Passed_GSE_information == "FAIL")
    )

    ## If there are any GSMs for which information was not obtained properly,
    ## rerun the .GSEInformationFunction on those samples and replace the info
    ## for those GSMs in the output DF.
    if (length(failedGSEInformationDFIters) > 0) {

        GSEInformationDFSecondPass <- as.data.frame(
            .GSEInformationFunction(GSMIDVector[failedGSEInformationDFIters])
        )

        GSEInformationDF[
            failedGSEInformationDFIters,
        ] <- GSEInformationDFSecondPass
    }

    ## Now run the .SRXInformationFunction using the results obtained from the
    ## previous .GSEInformationFunction to pull information about the SRX, such
    ## as the total reads, approximate bases and file sizes, species,
    ## sequencing instrument, etc.
    SRXInformationDF <- as.data.frame(
        .SRXInformationFunction(
            GSMIDValues = GSMIDVector,
            SRXIDURLs = GSEInformationDF$SRX_ID_URLs
        )
    )

    ## Do a second pass at downloading by taking the GSMs which failed the
    ## first time and re-running those through the .SRXInformationFunction.
    ## Note, since this function expands out, with more rows than the GSMs
    ## as multiple SRRs can be associated with a single GSM, we will need to
    ## identify which GSMs invocated a failure
    failedSRXInformationDFIters <- unique(
        which(SRXInformationDF$Passed_SRX_information == "FAIL")
    )

    ## If there are any SRXs for which information was not obtained properly,
    ## rerun the .SRXInformationFunction on the unique GSMs and SRX URLs which
    ## failed at least once, and replace all rows associated with those entries
    ## with the new data.
    if (length(failedSRXInformationDFIters) > 0) {

        ## These should be the same length, and essentially paired with each
        ## other, as each GSM has one corresponding SRX URL.
        uniqueGSMFailedSRXInformation <- unique(
            SRXInformationDF$GSM_IDs[failedSRXInformationDFIters]
        )

        uniqueSRXURLFailedSRXInformation <- unique(
            SRXInformationDF$SRX_ID_URLs[failedSRXInformationDFIters]
        )

        ## Rerun the results for the unique values that failed.
        SRXInformationDFSecondPass <- as.data.frame(
            .SRXInformationFunction(
                GSMIDValues = uniqueGSMFailedSRXInformation,
                SRXIDURLs = uniqueSRXURLFailedSRXInformation
            )
        )

        ## Replace all entries with at least one failure with the new results.
        SRXInformationDF[
            which(SRXInformationDF$GSM_IDs %in% uniqueGSMFailedSRXInformation),
        ] <- SRXInformationDFSecondPass
    }

    ## Do a right merge with the SRXInformationDF into the GSEInformationDF.
    informationDFMerge <- merge(
        GSEInformationDF,
        SRXInformationDF[
            , - which(colnames(SRXInformationDF) %in% "SRX_ID_URLs")
        ],
        by = "GSM_IDs",
        all.y = TRUE
    )

    ## Run the .SRRFileIdentificationFunction to get the paths to the SRR file
    ## URLs and the name of the SRR fastq files on the ENA server.
    SRRIdentificationDF <- as.data.frame(
        .SRRFileIdentificationFunction(
            GSMIDValues = informationDFMerge$GSM_IDs,
            SRRIDValues = informationDFMerge$SRR_IDs,
            SRRSequencingType = informationDFMerge$Sequencing_type
        )
    )

    ## Do a second pass at downloading by taking the samples which failed the
    ## first time and re-running those through the
    ## .SRRFileIdentificationFunction.
    ## Note, since this function expands out, with more rows than the GSMs
    ## as multiple SRRs can be associated with a single GSM, we will need to
    ## identify which GSMs invoked a failure.
    failedSRRIdentificationDFIters <- unique(
        which(SRRIdentificationDF$Passed_SRR_file_identification == "FAIL")
    )

    ## If there are any files which did not download information properly, do
    ## another pass at downloading them, and replace the information for those
    ## files
    if (length(failedSRRIdentificationDFIters) > 0) {

        ## Get the unique SRR IDs that failed.
        uniqueSRRFailedSRRIdentification <- unique(
            SRRIdentificationDF$SRR_IDs[failedSRRIdentificationDFIters]
        )

        ## Get the GSM IDs that are associated with those unique SRR IDs which
        ## failed.
        GSMFailedSRRIdentification <- informationDFMerge[
            which(
                informationDFMerge$SRR_IDs %in% uniqueSRRFailedSRRIdentification
            ),
            "GSM_IDs"
        ]

        ## Also get the sequencing types for each of the unique SRR IDs which
        ## failed.
        sequencingTypeFailedSRRIdentification <- informationDFMerge[
            which(
                informationDFMerge$SRR_IDs %in% uniqueSRRFailedSRRIdentification
            ),
            "Sequencing_type"
        ]

        ## Rerun the results for the unique SRR values (and associated data)
        ## that failed.
        SRRIdentificationDFSecondPass <- .SRRFileIdentificationFunction(
            GSMIDValues = GSMFailedSRRIdentification,
            SRRIDValues = uniqueSRRFailedSRRIdentification,
            SRRSequencingType = sequencingTypeFailedSRRIdentification
        )

        ## Replace all entries with at least one failure with the new results.
        SRRIdentificationDF[
            which(
                SRRIdentificationDF$SRR_IDs %in%
                    uniqueSRRFailedSRRIdentification
            ),
        ] <- SRRIdentificationDFSecondPass
    }

    ## Do a right merge with the SRRIdentificationDF into the
    ## informationDFMerge.
    SRRDFMerge <- merge(
        informationDFMerge,
        SRRIdentificationDF,
        by = "SRR_IDs",
        all.y = TRUE
    )

    ## Do some cleaning of duplicate columns caused by the merge.
    SRRDFMerge$GSM_IDs.y <- NULL

    colnames(SRRDFMerge)[
        which(colnames(SRRDFMerge) == "GSM_IDs.x")
    ] <- "GSM_IDs"

    ## If the user has opted, run the internal .SRRInformationFunction, which
    ## will pull the fastq sizes and the date they were first uploaded, as well
    ## as generate a 'Passed_SRR_file_information' column noting the status of
    ## the information pull.
    if (pullFastqFileSizes) {

        ## Run the .SRRInformationFunction to get the fastq file sizes and their
        ## date of upload.
        SRRInformationDF <- mapply(
            FUN = .SRRInformationFunction,
            GSMIDValues = SRRDFMerge$GSM_IDs,
            SRRIDValues = SRRDFMerge$SRR_IDs,
            SRRSequencingType = SRRDFMerge$Sequencing_type,
            ftpSRRFileURLs = SRRDFMerge$Fastq_URLs,
            SRRFastqNames = SRRDFMerge$Fastq_file_names,
            MoreArgs = list("errorWaitTime" = 30)
        )

        SRRInformationDF <- as.data.frame(t(SRRInformationDF))
        colnames(SRRInformationDF) <- c(
            "GSM_IDs",
            "SRR_IDs",
            "Sequencing_type",
            "Fastq_URLs",
            "Fastq_file_names",
            "Fastq_file_sizes",
            "Fastq_file_first_upload_date",
            "Passed_SRR_file_information"
        )

        ## Do a second pass at getting fastq data by taking the fastqs which
        ## failed the first time.
        failedSRRInformationDFIters <- unique(
            which(SRRInformationDF$Passed_SRR_file_information == "FAIL")
        )

        ## If there are any files which did not pull information properly, do
        ## another pass at using the .SRRInformationFunction function on them,
        ## and replace the information for those files.
        if (length(failedSRRInformationDFIters) > 0) {

            SRRInformationDFSecondPass <- mapply(
                FUN = .SRRInformationFunction,
                GSMIDValues = SRRDFMerge$GSM_IDs[
                    failedSRRInformationDFIters
                ],
                SRRIDValues = SRRDFMerge$SRR_IDs[
                    failedSRRInformationDFIters
                ],
                SRRSequencingType = SRRDFMerge$Sequencing_type[
                    failedSRRInformationDFIters
                ],
                ftpSRRFileURLs = SRRDFMerge$Fastq_URLs[
                    failedSRRInformationDFIters
                ],
                SRRFastqNames = SRRDFMerge$Fastq_file_names[
                    failedSRRInformationDFIters
                ],
                MoreArgs = list("errorWaitTime" = 30)
            )

            SRRInformationDFSecondPass <- as.data.frame(
                t(SRRInformationDFSecondPass)
            )

            colnames(SRRInformationDFSecondPass) <- c(
                "GSM_IDs",
                "SRR_IDs",
                "Sequencing_type",
                "Fastq_URLs",
                "Fastq_file_names",
                "Fastq_file_sizes",
                "Fastq_file_first_upload_date",
                "Passed_SRR_file_information"
            )

            SRRInformationDF[
                failedSRRInformationDFIters,
            ] <- SRRInformationDFSecondPass
        }

        ## Add columns which aren't duplicates from the SRRInformationDF to the
        ## SRRDFMerge
        SRRDFMerge <- cbind(
            SRRDFMerge,
            SRRInformationDF[,-c(1:5)]
        )

        ## Convert the Fastq_file_sizes column to numeric (as it gets saved as a
        ## character).
        SRRDFMerge$Fastq_file_sizes <- as.numeric(
            SRRDFMerge$Fastq_file_sizes
        )

    } else {

        ## If the user hasn't elected to do the information pull, create the
        ## columns which would be pulled, and assign NA values to them.
        SRRDFMerge$Fastq_file_sizes <- NA
        SRRDFMerge$Fastq_file_first_upload_date <- NA
        SRRDFMerge$Passed_SRR_file_information <- NA
    }

    ## Rearrange the columns so the ones associated with the function checks
    ## are at the end.
    SRRDFMerge <- SRRDFMerge[
        c(2:7,1,9:17,19:20,22:23,8,18,21,24)
    ]

    ## Return the dataframe to the user.
    return(SRRDFMerge)
}
