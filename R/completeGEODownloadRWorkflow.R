#' A complete workflow of GEODownloadR functions
#'
#' This function acts as a wrapper which combines the three GEODownloadR
#' functions (identifySRRsFromGSMIDs, SRRFastqDownload, and optionally
#' combineSRRFastqsToGSMs) into a single streamlined function. This function
#' will also include automatic rerunning on GSM IDs that fail to download or
#' combine properly, and will issue a message to the user (along with the usual
#' dataframe of information output by either the SRRFastqDownload or
#' combineSRRFastqsToGSMs functions) noting the GSM IDs which have still failed
#' after this function's checks, so the user can more closely inspect those IDs.
#'
#' @param GSMIDVector Specify a vector of sample accession numbers to
#' investigate/download files for. The values in the GSMIDVector should be given
#' in the form "GSM" followed by the numbers associated with the sample, such as
#' "GSM2437766".
#' @param useCombineSRRFastqsToGSMs Set to TRUE to use the
#' combineSRRFastqsToGSMs function to combine the downloaded SRR .fastq files
#' back into .fastq files associated with each GSM ID directly. Defaults to
#' TRUE.
#' @param localDownloadDirectory Set a path to a directory where the user wishes
#' to download the .fastq files associated with the SRRs. Defaults to the user's
#' working directory.
#' @param outputGSMFastqDirectory Set a path to a directory where the user wants
#' the combined .fastq files for the GSM IDs to be output to. Defaults to the
#' user's working directory.
#' @param downloadSRRWithGSMPresent Set to TRUE to try downloading SRR .fastq
#' files even when the final GSM .fastq files associated with those SRRs
#' (compiled by the 'combineSRRFastqsToGSMs' function) is present, otherwise set
#' to FALSE to not. Defaults to FALSE.
#' @param SRRInformationNodeCount Specify a maximum number of nodes to use when
#' getting SRR information (particularly the SRR file sizes). This will only
#' have an effect if `pullFastqFileSizes` is TRUE. Try setting to a relatively
#' low number (like 1) if the function is having difficulty in acquiring the
#' fastq file size (and there are many "FAILED" values in the
#' 'Passed_SRR_file_information' column of the returned data frame). Defaults to
#' 12.
#' @param SRRDownloadNodeCount Specify a maximum number of nodes to use when
#' downloading SRR .fastqs. Try setting to a relatively
#' low number (like 1) if the function is having difficulty in acquiring the
#' fastq files (and there are many "FAILED" values in the
#' 'SRR_file_download_status' column of the returned data frame). Defaults to
#' 12.
#' @return First, this function finds and downloads the SRR .fastq files
#' associated with the user's GSM IDs given as the `GSMIDVector` to the
#' directory specified as the `localDownloadDirectory`, then, if
#' `useCombineSRRFastqsToGSMs` is set to TRUE, this function combines these into
#' .fastq files for each GSM ID as a whole. These GSM IDs will be deposited in
#' the directory specified as the `outputGSMFastqDirectory` argument. The
#' function will also return a data frame with information related to the GSMs,
#' SRR files, and the status of each function that was run.
#' @export
#'
#' @examplesIf interactive()
#' ## This example will attempt to pull, download, and combine information for
#' ## three example GSM IDs and returns a data frame with that information to
#' ## the user. Running this example requires an internet connection, and it
#' ## will attempt to download and combine files to the user's working
#' ## directory. It will also use 8 nodes when pulling information and
#' ## downloading files.
#'
#' returnValue <- identifySRRsFromGSMIDs(
#'     GSMIDVector = c("GSM4110236", "GSM2774971", "GSM2406903"),
#'     useCombineSRRFastqsToGSMs = TRUE,
#'     SRRInformationNodeCount = 8
#'     SRRDownloadNodeCount = 8
#' )

completeGEODownloadRWorkflow <- function(
    GSMIDVector,
    useCombineSRRFastqsToGSMs = TRUE,
    localDownloadDirectory = getwd(),
    outputGSMFastqDirectory = getwd(),
    downloadSRRWithGSMPresent = FALSE,
    SRRInformationNodeCount = 12,
    SRRDownloadNodeCount = 12
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

    ## Verify the useCombineSRRFastqsToGSMs argument.
    if (!.isSingularBoolean(useCombineSRRFastqsToGSMs)) {
        stop(
            "The `useCombineSRRFastqsToGSMs` argument must be a single ",
            "boolean value."
        )
    }

    ## Check that the localDownloadDirectory and outputGSMFastqDirectory are
    ## valid paths.

    ## First ensure the paths are single character values.
    if (!.isSingularCharacter(localDownloadDirectory)) {
        stop(
            "The `localDownloadDirectory` argument must be a single character ",
            "element noting a path where the .fastq files associated with the ",
            "SRRs should be downloaded to."
        )
    }

    if (!.isSingularCharacter(outputGSMFastqDirectory)) {
        stop(
            "The `outputGSMFastqDirectory` argument must be a single character ",
            "element noting a path where the .fastq files associated with the ",
            "GSMs created by combining their component SRR .fastq files ",
            "should be saved to."
        )
    }

    ## Next ensure that the paths exist.
    if (!dir.exists(localDownloadDirectory)) {
        stop(
            "The path given as the `localDownloadDirectory` argument doesn't ",
            "seem to exist. Please give a valid path to download .fastq files ",
            "associated with the SRRs to."
        )
    }

    if (!dir.exists(outputGSMFastqDirectory)) {
        stop(
            "The path given as the `outputGSMFastqDirectory` argument doesn't ",
            "seem to exist. Please give a valid path to save the combined ",
            "GSM .fastq files to."
        )
    }

    ## Finally, make sure there is no trailing "/" in the paths, and remove them
    ## if there are.
    localDownloadDirectory <- .trailingFSlashRemover(localDownloadDirectory)

    outputGSMFastqDirectory <- .trailingFSlashRemover(outputGSMFastqDirectory)

    ## Verify the downloadSRRWithGSMPresent argument.
    if (!.isSingularBoolean(downloadSRRWithGSMPresent)) {
        stop(
            "The `downloadSRRWithGSMPresent` argument must be a single ",
            "boolean value."
        )
    }

    ## Finaly, verify the SRRInformationNodeCount and SRRDownloadNodeCount
    ## arguments.
    if (!.isSinglePositiveInteger(SRRInformationNodeCount)) {
        stop(
            "The `SRRInformationNodeCount` argument must be a singular ",
            "positive integer value."
        )
    }

    if (!.isSinglePositiveInteger(SRRDownloadNodeCount)) {
        stop(
            "The `SRRDownloadNodeCount` argument must be a singular ",
            "positive integer value."
        )
    }

    ## Now start with the identifySRRsFromGSMIDs function.
    resultsDF <- GEODownloadR::identifySRRsFromGSMIDs(
        GSMIDVector = GSMIDVector,
        pullFastqFileSizes = TRUE,
        SRRInformationNodeCount = SRRInformationNodeCount
    )

    ## Use the data from the identifySRRsFromGSMIDs function to then download
    ## the SRR fastqs with the SRRFastqDownload function.
    resultsDF <- GEODownloadR::SRRFastqDownload(
        resultsDF,
        localDownloadDirectory = localDownloadDirectory,
        downloadSRRWithGSMPresent = downloadSRRWithGSMPresent,
        SRRDownloadNodeCount = SRRDownloadNodeCount
    )

    ## If the user has set useCombineSRRFastqsToGSMs to TRUE, run the
    ## useCombineSRRFastqsToGSMs function to combine SRR fastqs into those per
    ## GSM ID, then identify the unique GSM IDs which failed to either download
    ## SRRs or combine back into GSMs.
    if (useCombineSRRFastqsToGSMs) {

        ## Run the combineSRRFastqsToGSMs function.
        resultsDF <- GEODownloadR::combineSRRFastqsToGSMs(
            SRRFastqDownloadDF = resultsDF,
            localDownloadDirectory = localDownloadDirectory,
            outputGSMFastqDirectory = outputGSMFastqDirectory
        )

        ## Identify the SRR fastq file names which failed to combine.
        failedSRRFiles <- resultsDF[
            !(
                resultsDF$GSMCombinationStatus %in% c(
                    "PASS", "GSM_FASTQ_FILE_ALREADY_PRESENT"
                )
            ),
            "Fastq_file_names"
        ]

        ## Identify the unique GSM IDs which failed to combine.
        failedGSMIDs <- unique(
            resultsDF[
                !(
                    resultsDF$GSMCombinationStatus %in% c(
                        "PASS", "GSM_FASTQ_FILE_ALREADY_PRESENT"
                    )
                ),
                "GSM_IDs"
            ]
        )

        ## Issue a message to the user noting the number and names of the SRR
        ## .fastq files and the GSM IDs which failed and will be rerun.
        if (length(failedGSMIDs) > 0) {
            message(
                paste0(
                    "As many as ",
                    length(failedSRRFiles),
                    " SRR .fastq files, and ",
                    length(failedGSMIDs),
                    " unique GSMs failed to properly combine. We will now ",
                    "attempt to redownload and recombine these files."
                )
            )
        }
    } else {

        ## If useCombineSRRFastqsToGSMs is FALSE, only get the SRR .fastq names
        ## and GSM IDs which failed the SRR download step instead.
        failedSRRFiles <- resultsDF[
            !(
                resultsDF$SRR_file_download_status %in% c(
                    "PASS",
                    "SRR_FILE_ALREADY_DOWNLOADED"
                )
            ),
            "Fastq_file_names"
        ]

        failedGSMIDs <- unique(
            resultsDF[
                !(
                    resultsDF$SRR_file_download_status %in% c(
                        "PASS",
                        "SRR_FILE_ALREADY_DOWNLOADED"
                    )
                ),
                "GSM_IDs"
            ]
        )

        ## Issue a message to the user noting the number and names of the SRR
        ## .fastq files and the GSM IDs which failed and will be rerun.
        if (length(failedGSMIDs) > 0) {
            message(
                paste0(
                    "As many as ",
                    length(failedSRRFiles),
                    " SRR .fastq files from ",
                    length(failedGSMIDs),
                    " unique GSMs failed to properly download. We will now ",
                    "attempt to redownload these files."
                )
            )
        }
    }

    ## If there are some GSM IDs that failed, rerun the steps on those IDs,
    ## and replace the results for them with the rerun ones.

    ## Remove the previous results for the unique GSM IDs that failed
    resultsDF <- resultsDF[
        !(resultsDF$GSM_IDs %in% failedGSMIDs),
    ]

    ## Rerun the identifySRRsFromGSMIDs function on the failed IDs, then
    ## redownload those samples with the SRRFastqDownload function.
    if (length(failedGSMIDs) > 0) {

        failedGSMIDsResultsDF <- GEODownloadR::identifySRRsFromGSMIDs(
            GSMIDVector = failedGSMIDs,
            pullFastqFileSizes = TRUE,
            SRRInformationNodeCount = SRRInformationNodeCount
        )

        failedGSMIDsResultsDF <- SRRFastqDownload(
            failedGSMIDsResultsDF,
            localDownloadDirectory = localDownloadDirectory,
            downloadSRRWithGSMPresent = downloadSRRWithGSMPresent,
            SRRDownloadNodeCount = SRRDownloadNodeCount
        )

        ## If the user has set useCombineSRRFastqsToGSMs to TRUE, run the
        ## useCombineSRRFastqsToGSMs function again.
        if (useCombineSRRFastqsToGSMs) {

            ## Run the combineSRRFastqsToGSMs function.
            failedGSMIDsResultsDF <- GEODownloadR::combineSRRFastqsToGSMs(
                SRRFastqDownloadDF = failedGSMIDsResultsDF,
                localDownloadDirectory = localDownloadDirectory,
                outputGSMFastqDirectory = outputGSMFastqDirectory
            )

            ## Identify the SRR fastq file names which failed to combine.
            failedSRRFiles2 <- failedGSMIDsResultsDF[
                !(
                    failedGSMIDsResultsDF$GSMCombinationStatus %in% c(
                        "PASS", "GSM_FASTQ_FILE_ALREADY_PRESENT"
                    )
                ),
                "Fastq_file_names"
            ]

            ## Identify the unique GSM IDs which failed to combine.
            failedGSMIDs2 <- unique(
                failedGSMIDsResultsDF[
                    !(
                        failedGSMIDsResultsDF$GSMCombinationStatus %in% c(
                            "PASS", "GSM_FASTQ_FILE_ALREADY_PRESENT"
                        )
                    ),
                    "GSM_IDs"
                ]
            )

            ## Issue a message to the user noting the number and names of the SRR
            ## .fastq files and the GSM IDs which failed and will be rerun.
            if (length(failedGSMIDs2) > 0) {
                message(
                    paste0(
                        "As many as ",
                        length(failedSRRFiles2),
                        " SRR .fastq files, and ",
                        length(failedGSMIDs2),
                        " unique GSMs still failed to properly combine. ",
                        "Please check these entries - it is possible they may ",
                        "not exist or there may be an issue with the data ",
                        "that is saved for these GSMs. It may also be ",
                        "possible the server connection is compromised. ",
                        "Consider rerunning this function perhaps with lower ",
                        "node counts."
                    )
                )
            }

        } else {

            ## If useCombineSRRFastqsToGSMs is FALSE, only get the SRR .fastq
            ## names and GSM IDs which failed the SRR download step instead.
            failedSRRFiles2 <- failedGSMIDsResultsDF[
                !(
                    failedGSMIDsResultsDF$SRR_file_download_status %in% c(
                        "PASS",
                        "SRR_FILE_ALREADY_DOWNLOADED"
                    )
                ),
                "Fastq_file_names"
            ]

            failedGSMIDs2 <- unique(
                failedGSMIDsResultsDF[
                    !(
                        failedGSMIDsResultsDF$SRR_file_download_status %in% c(
                            "PASS",
                            "SRR_FILE_ALREADY_DOWNLOADED"
                        )
                    ),
                    "GSM_IDs"
                ]
            )

            ## Issue a message to the user noting the number and names of the SRR
            ## .fastq files and the GSM IDs which failed and will be rerun.
            if (length(failedGSMIDs2) > 0) {
                message(
                    paste0(
                        "As many as ",
                        length(failedSRRFiles2),
                        " SRR .fastq files, and ",
                        length(failedGSMIDs2),
                        " unique GSMs still failed to properly download. ",
                        "Please check these entries - it is possible they may ",
                        "not exist or there may be an issue with the data ",
                        "that is saved for these GSMs. It may also be ",
                        "possible the server connection is compromised. ",
                        "Consider rerunning this function perhaps with lower ",
                        "node counts."
                    )
                )
            }
        }

        ## Add the results for the reprocessed files in with the previous
        ## resultsDF.
        resultsDFOutput <- rbind(
            resultsDF,
            failedGSMIDsResultsDF
        )

        ## Sort the dataset by GSM ID then SRR .fastq name.
        resultsDFOutput <- resultsDFOutput[
            order(resultsDFOutput$GSM_IDs, resultsDFOutput$Fastq_file_names),
        ]

    } else {

        ## Message the user to inform them that no errors were detected.
        message("None of the specified GSM IDs failed to download or combine")

        ## Rename the resultsDF since there were no errors.
        resultsDFOutput <- resultsDF
    }

    ## Return the final resultsDFOutput
    return(resultsDFOutput)
}
