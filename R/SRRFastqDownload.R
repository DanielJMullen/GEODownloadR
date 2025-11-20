## Internal functions for the SRRFastqDownload function:

## An internal function which takes the ftp paths to .fastq files on the ENA
## server, along with other data including the file sizes, and downloads them
## to a specified directory.
.internalSRRDownloadFunction <- function(
    GSMIDValues,
    SRRIDValues,
    SRRSequencingType,
    ftpSRRFileURLs,
    SRRFastqNames,
    SRRFastqFileSizes,
    downloadErrorWaitTime, # This exsits to easily modulate in the future
    internalDownloadSRRWithGSMPresent = FALSE,
    internalLocalDownloadDirectory = NULL
) {

    ## Initialize a holding value for the downloadIndices, which will be used
    ## shortly when setting up parallelization using foreach().
    downloadIndices <- NULL

    ## Create a list of the nodes that will be used, and then register them for
    ## use with the foreach package. 12 nodes will be used by default.
    makeDownloadClusterList <- parallel::makeCluster(12)
    doParallel::registerDoParallel(makeDownloadClusterList)

    ## Use foreach() to iterate across the SRR URLs, and download the file to
    ## each, then return information on the success (or failure) of the
    ## download.
    downloadFileVec <- foreach::foreach(
        downloadIndices = seq_along(ftpSRRFileURLs),
        .combine = 'c'
    ) %dopar% {

        ## Get the associated GSM ID value, akin to a for loop.
        iterGSMIDValue <- GSMIDValues[downloadIndices]

        ## Get the individual SRR ID to use for this iteration, akin to a
        ## for loop.
        iterSRRIDValue <- SRRIDValues[downloadIndices]

        ## Get the type of sequencing (single or paired end) that was done for
        ## the given SRR.
        iterSRRSequencingTypeValue <- SRRSequencingType[downloadIndices]

        ## Get the individual fastq file URL.
        iterFTPSRRFileURL <- ftpSRRFileURLs[downloadIndices]

        ## Get the individual file fastq name.
        iterSRRFastqName <- SRRFastqNames[downloadIndices]

        ## Finally, get the individual fastq file size.
        iterSRRFastqFileSize <- SRRFastqFileSizes[downloadIndices]

        ## Check if the iterFTPSRRFileURL is NA. If so, we will need to cancel
        ## out of this step to prevent an error from being induced as we can't
        ## load an NA URL.
        if(is.na(iterFTPSRRFileURL)) {
            return("NO_FTP_URL_GIVEN")
        }

        ## Also check if the iterSRRFastqFileSize is NA. If so, we will need to
        ## cancel out of this step to prevent an error from being induced as we
        ## can't verify file is properly downloaded without knowing it's size.
        if(is.na(iterSRRFastqFileSize)) {
            return("FILE_SIZE_NOT_GIVEN")
        }

        ## Create a variable noting where the SRR file will be downloaded to.
        iterSRRFileDownloadPath <- file.path(
            internalLocalDownloadDirectory,
            iterSRRFastqName
        )

        ## Check that the given SRR file doesn't already exist, and that if it
        ## does it is of the correct length. If it is not either of these, try
        ## downloading the file.
        if (
            !file.exists(iterSRRFileDownloadPath) ||
            file.info(
                iterSRRFileDownloadPath
            )[1,1] != iterSRRFastqFileSize
        ) {

            ## If the overall GSM is present and
            ## internalDownloadSRRWithGSMPresent is FALSE, also check to make
            ## sure the given GSM file isn't present
            if (!internalDownloadSRRWithGSMPresent) {

                ## Determine the paired extension for the GSM.
                if (iterSRRSequencingTypeValue == "PAIRED") {

                    ## This check exists because there seems to be rare
                    ## repositories like GSM2856807 - SRR6290078 where there are
                    ## two normal read end SRAs listed (SRR6290078_1.fastq.gz
                    ## and SRR6290078_2.fastq.gz) but there is an additional SRA
                    ## (SRR6290078.fastq.gz) which doesn't seem to have a pair.
                    ## This will help to correct for this to ensure the extra
                    ## file is not downloaded. when either paired GSM is
                    ## present. NOTE: suppressWarnings is used here to prempt a
                    ## rare cirumstance where NAs may be induced by coercion.
                    ## This does not affect the functionality of the code
                    ## however.
                    iterFileNameNoExtension <- suppressWarnings(
                        sub(
                            '\\.fastq.gz.*',
                            '',
                            iterSRRFastqName
                        )
                    )

                    if (
                        substr(
                            iterFileNameNoExtension,
                            (nchar(iterFileNameNoExtension)-1),
                            (nchar(iterFileNameNoExtension)-1)
                        ) == "_"
                    ) {

                        ## These are files where the paired file extension (_1
                        ## or _2) is listed as the last two characters of the
                        ## name before the .fastq.gz at the end.
                        iterPairedFilePortion <- substr(
                            iterFileNameNoExtension,
                            (nchar(iterFileNameNoExtension)-1),
                            iterFileNameNoExtension
                        )
                    } else {

                        ## These are examples of rare files where they do not
                        ## have a paired file extension. In these cases, we will
                        ## not download the file if either paired GSM file
                        ## exists.
                        iterPairedFilePortion <- c("_1", "_2")
                    }

                    iterGSMFileName <- paste0(
                        iterGSMIDValue,
                        iterPairedFilePortion,
                        ".fastq.gz"
                    )
                } else {

                    ## This is the name of the single GSM file.
                    iterGSMFileName <- paste0(iterGSMIDValue, ".fastq.gz")
                }

                ## Check if the GSM file associated with the given SRR has been
                ## downloaded already, if it hasn't, download the SRR. If it
                ## has, or the size of theexisting file doesn't equal the size
                ## the file should be, skip the download of this SRR with a note
                ## why.
                if (
                    all(!file.exists(
                        file.path(
                            internalLocalDownloadDirectory,
                            iterGSMFileName
                        )
                    ))
                ) {

                    iterSRRDownloadFailCount <- 0

                    while (
                        iterSRRDownloadFailCount < 3 &&
                        (
                            !file.exists(iterSRRFileDownloadPath) ||
                            file.info(
                                iterSRRFileDownloadPath
                            )[1,1] != iterSRRFastqFileSize
                        )
                    ) {

                        ## Try downloading the SRR file from the ftp link. If it
                        ## doesn't download properly, wait 5 seconds before
                        ## attempting again.
                        ## supressMessages is used to suppress the output from
                        ## the curl implementation of download.file() as it
                        ## relates only to the status of the download.
                        tryCatch ({
                            # supressMessages(
                            download.file(
                                iterFTPSRRFileURL,
                                destfile = iterSRRFileDownloadPath,
                                method = "curl",
                                mode = "wb"
                            )
                            # )
                        }, error = function(cond) {
                            Sys.sleep(downloadErrorWaitTime)
                        })

                        ## Increment the iterSRRDownloadFailCount.
                        iterSRRDownloadFailCount <- (
                            iterSRRDownloadFailCount + 1
                        )
                    }

                    ## On success or error out, wait 5 seconds standard.
                    Sys.sleep(5)

                    ## Do a final check to see that the file was downloaded
                    ## and is of the proper size:
                    if (!file.exists(iterSRRFileDownloadPath)) {

                        iterFileCheckReturnValue <- "FAILED_DOWNLOAD"
                    } else {

                        if (
                            file.info(
                                iterSRRFileDownloadPath
                            )[1,1] != iterSRRFastqFileSize
                        ) {
                            iterFileCheckReturnValue <-
                                "FAILED_DOWNLOAD_FILE_TRUNCATED"
                        } else{
                            iterFileCheckReturnValue <- "PASS"
                        }
                    }

                    ## Do the return:
                    return(iterFileCheckReturnValue)

                } else {
                    ## Do the return noting the related GSM file is already
                    ## present as a reason for not downloading the SRR again.
                    return("GSM_FILE_ALREADY_DOWNLOADED")
                }

            } else{

                ## Proceed with the download regardless of if the related GSM
                ## file is present or not.
                iterSRRDownloadFailCount <- 0

                while (
                    iterSRRDownloadFailCount < 3 &&
                    (
                        !file.exists(iterSRRFileDownloadPath) ||
                        file.info(
                            iterSRRFileDownloadPath
                        )[1,1] != iterSRRFastqFileSize
                    )
                ) {
                    ## Try downloading the SRR file from the ftp link. If it
                    ## doesn't download properly, wait 5 seconds before
                    ## attempting again.
                    ## supressMessages is used to suppress the output from the
                    ## curl implementation of download.file() as it relates
                    ## only to the status of the download.
                    tryCatch({
                        # suppressMessages(
                        download.file(
                            iterFTPSRRFileURL,
                            destfile = iterSRRFileDownloadPath,
                            method = "curl",
                            mode = "wb"
                        )
                        # )
                    }, error = function(cond) {
                        Sys.sleep(downloadErrorWaitTime)
                    })

                    ## Increment the iterSRRDownloadFailCount.
                    iterSRRDownloadFailCount <- (iterSRRDownloadFailCount + 1)
                }

                ## On success or error out, wait 5 seconds standard.
                Sys.sleep(5)

                ## Do a final check to see that the file was downloaded
                ## and is of the proper size:
                if (!file.exists(iterSRRFileDownloadPath)) {

                    iterFileCheckReturnValue <- "FAILED_DOWNLOAD"
                } else {

                    if (
                        file.info(
                            iterSRRFileDownloadPath
                        )[1,1] != iterSRRFastqFileSize
                    ) {
                        iterFileCheckReturnValue <-
                            "FAILED_DOWNLOAD_FILE_TRUNCATED"
                    } else{
                        iterFileCheckReturnValue <- "PASS"
                    }
                }

                ## Do the return:
                return(iterFileCheckReturnValue)
            }

        } else{

            ## Create a return noting the SRR has already been properly
            ## downloaded.
            return("SRR_FILE_ALREADY_DOWNLOADED")
        }
    }

    ## Close the cluster list.
    ## suppressWarnings is used to silence a warning that relates to closing
    ## unused connections, which doesn't affect how the function works.
    suppressWarnings(parallel::stopCluster(makeDownloadClusterList))

    ## Return the final vector.
    return(downloadFileVec)
}

## Download .fastq files associated with SRR IDs:

#' This function takes a data frame with information on SRR IDs and their
#' associated .fastq files from the European Nucleotide Archive, such as one
#' generated by the 'combineSRRFastqsToGSMs' function, and downloads the .fastq
#' files associated with them to a user-selected directory.
#'
#' @param identifySRRsFromGSMIDsDF Specify a dataframe with relevant information
#' about the SRRs and their associated .fastq files, such as one output by the
#' identifySRRsFromGSMIDs function. This dataframe should contain columns with
#' the names 'GSM_IDs', 'SRR_IDs', 'Sequencing_type', 'Fastq_URLs',
#' 'Fastq_file_names', and 'Fastq_file_sizes'.
#' @param localDownloadDirectory Set a path to a directory where the user wishes
#' to download the .fastq files associated with the SRRs detailed in the
#' `SRRIdentificationDF`. Defaults to the user's working directory.
#' @param downloadSRRWithGSMPresent Set to TRUE to try downloading SRR .fastq
#' files even when the final GSM .fastq files associated with those SRRs
#' (compiled by the 'combineSRRFastqsToGSMs' function) is present, otherwise set
#' to FALSE to not. Defaults to FALSE.
#' @return Downloads and checks the .fastq files to the directory specified as
#' the `localDownloadDirectory`. The function will also return the dataframe
#' given as the `identifySRRsFromGSMIDsDF` argument with an additional column
#' named "SRR_file_download_status" noting the status of the download of each
#' .fastq file.
#' @export
#'
#' @examplesIf interactive()
#' ## This example requires an internet connection, and, if run, will attempt to
#' ## download the .fastq files for SRRs associated with three example GSM IDs
#' ## using an example dataset with results from the previous
#' ## 'identifySRRsFromGSMIDs' function to the user's working directory.
#'
#' ## Load the example dataset.
#' exampleDFGEODownloadR <- utils::data(
#'     "exampleDFGEODownloadR",
#'     package = "GEODownloadR"
#' )
#'
#' ## Download the SRR .fastqs. There should be a total of 11 .fastq files
#' ## donwloaded.
#' returnValue <- SRRFastqDownload(
#'     identifySRRsFromGSMIDsDF = exampleDFGEODownloadR
#' )

SRRFastqDownload <- function(
    identifySRRsFromGSMIDsDF,
    localDownloadDirectory = getwd(),
    downloadSRRWithGSMPresent = FALSE
) {

    ## Verify the identifySRRsFromGSMIDsDF.

    ## First ensure the object that is given to the identifySRRsFromGSMIDsDF
    ## argument is a dataframe (or matrix could work).
    .isDataFrameOrMatrix(identifySRRsFromGSMIDsDF)

    ## Then, ensure all the relevant columns are present in the data frame or
    ## matrix given as the identifySRRsFromGSMIDsDF.
    .verifyColumnsPresent(
        c(
            "GSM_IDs",
            "SRR_IDs",
            "Sequencing_type",
            "Fastq_URLs",
            "Fastq_file_names",
            "Fastq_file_sizes"
        ),
        identifySRRsFromGSMIDsDF
    )

    ## Make sure the "Fastq_file_sizes" column is numeric.
    identifySRRsFromGSMIDsDF$Fastq_file_sizes <- as.numeric(
        identifySRRsFromGSMIDsDF$Fastq_file_sizes
    )

    ## Ensure that there are still valid fastq file sizes left, and issue an
    ## error to the user if not.
    if(all(is.na(identifySRRsFromGSMIDsDF$Fastq_file_sizes))) {
        stop(
            "The values in the 'Fastq_file_sizes' column of the ",
            "dataframe given as the `identifySRRsFromGSMIDsDF` argument are ",
            "not numeric. Please ensure numeric values are given in the ",
            "'Fastq_file_sizes' column and rerun this function."
        )
    }

    ## Check that the localDownloadDirectory is a valid path.

    ## First ensure the localDownloadDirectory is a character vector.
    if (!.isSingularCharacter(localDownloadDirectory)) {
        stop(
            "The `localDownloadDirectory` argument must be a single character ",
            "element noting a path where the .fastq files associated with the ",
            "SRRs in the `identifySRRsFromGSMIDsDF` dataframe should be ",
            "downloaded to."
        )
    }

    ## Next ensure that the path exists.
    if (!file.exists(localDownloadDirectory)) {
        stop(
            "The path given as the `localDownloadDirectory` argument doesn't ",
            "seem to exist. Please give a valid path to download .fastq files ",
            "associated with the SRRs in the `identifySRRsFromGSMIDsDF` ",
            "dataframe to."
        )
    }

    ## Finally, make sure there is no trailing "/" in the path, and remove it
    ## if there is.
    localDownloadDirectory <- .trailingFSlashRemover(localDownloadDirectory)

    ## Verify the downloadSRRWithGSMPresent argument.
    if (!.isSingularBoolean(downloadSRRWithGSMPresent)) {
        stop(
            "The `downloadSRRWithGSMPresent` argument must be a singular ",
            "boolean value."
        )
    }

    ## Get the user's options:
    original_options <- options()

    ## Reset the timeout value so we don't lose the file download:
    options(timeout = max(12000, getOption("timeout")))

    ## Set a trigger to restore the user's options when the function ends.
    on.exit(options(original_options))

    ## Download the fastq files associated with the SRR IDs.
    SRRDownloadVec <- .internalSRRDownloadFunction(
        GSMIDValues = identifySRRsFromGSMIDsDF$GSM_IDs,
        SRRIDValues = identifySRRsFromGSMIDsDF$SRR_IDs,
        SRRSequencingType = identifySRRsFromGSMIDsDF$Sequencing_type,
        ftpSRRFileURLs = identifySRRsFromGSMIDsDF$Fastq_URLs,
        SRRFastqNames = identifySRRsFromGSMIDsDF$Fastq_file_names,
        SRRFastqFileSizes = identifySRRsFromGSMIDsDF$Fastq_file_sizes,
        downloadErrorWaitTime = 30,
        internalDownloadSRRWithGSMPresent = downloadSRRWithGSMPresent,
        internalLocalDownloadDirectory = localDownloadDirectory
    )

    ## Do a second check and run for files that didn't download properly.
    failedSRRDownloadVecIters <- unique(
        which(!(SRRDownloadVec %in% c("PASS", "SRR_FILE_ALREADY_DOWNLOADED")))
    )

    if (length(failedSRRDownloadVecIters) > 0) {

        ## If there are files that didn't download properly, redownload just
        ## those files.
        SRRDownloadVecSecondPass <- .internalSRRDownloadFunction(
            GSMIDValues = identifySRRsFromGSMIDsDF$GSM_IDs[
                failedSRRDownloadVecIters
            ],
            SRRIDValues = identifySRRsFromGSMIDsDF$SRR_IDs[
                failedSRRDownloadVecIters
            ],
            SRRSequencingType = identifySRRsFromGSMIDsDF$Sequencing_type[
                failedSRRDownloadVecIters],
            ftpSRRFileURLs = identifySRRsFromGSMIDsDF$Fastq_URLs[
                failedSRRDownloadVecIters
            ],
            SRRFastqNames = identifySRRsFromGSMIDsDF$Fastq_file_names[
                failedSRRDownloadVecIters
            ],
            SRRFastqFileSizes = identifySRRsFromGSMIDsDF$Fastq_file_sizes[
                failedSRRDownloadVecIters
            ],
            downloadErrorWaitTime = 30,
            internalDownloadSRRWithGSMPresent = downloadSRRWithGSMPresent,
            internalLocalDownloadDirectory = localDownloadDirectory
        )

        ## Add the results for the redownloaded files back into the main
        ## SRRDownloadVec
        SRRDownloadVec[failedSRRDownloadVecIters] <- SRRDownloadVecSecondPass
    }

    ## Add the results of the download to the identifySRRsFromGSMIDsDF
    ## dataframe.
    identifySRRsFromGSMIDsDF$SRR_file_download_status <- SRRDownloadVec

    ## Return the data frame back to the user.
    return(identifySRRsFromGSMIDsDF)
}
