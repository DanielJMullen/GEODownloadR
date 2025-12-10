## Internal functions for the combineSRRFastqsToGSMs function

## An internal function which combines downloaded SRR files into files assigned
## per GSM.
## Note: This function does assume the SRR files it is trying to combined are
## gzipped, and thus the resulting GSM fastq file that is made with have the
## .gz extension.
.combineSRRToGSMFunction <- function(
        SRRInfoAndDownloadFunctionDF,
        internalLocalDownloadDirectory,
        internalOutputGSMFastqDirectory
) {

    ## Get the unique GSM values that appear in the dataframe output by the
    ## .SRRInfoAndDownloadFunction.
    uniqueGSMIDValues <- unique(SRRInfoAndDownloadFunctionDF$GSM_IDs)

    ## Initialize a holding value for the uniqueGSMIndices, which will be used
    ## shortly when setting up parallelization using foreach().
    uniqueGSMIndices <- NULL

    ## Create a list of the nodes that will be used, and then register them for
    ## use with the foreach package. Up to 60 nodes will be used by default.
    combineGSMClusterList <- parallel::makeCluster(
        min(length(uniqueGSMIDValues), 60)
    )
    doParallel::registerDoParallel(combineGSMClusterList)

    ## Use foreach() to iterate across the GSM ID values, and extract the
    ## relevant information from each, then combine this information into a
    ## dataframe using c(). The relevant information here includes ...
    combineGSMVec <- foreach::foreach(
        uniqueGSMIndices = seq_along(uniqueGSMIDValues),
        .combine = 'c'
    ) %dopar% {

        ## Get the individual GSMID to use for this iteration, akin to a
        ## for loop.
        iterUniqueGSMID <- uniqueGSMIDValues[uniqueGSMIndices]

        ## Get the local SRR fastqs that should exist for the given GSM.
        localSRRFilesPaths <- file.path(
            internalLocalDownloadDirectory,
            paste0(
                SRRInfoAndDownloadFunctionDF[
                    SRRInfoAndDownloadFunctionDF$GSM_IDs == iterUniqueGSMID,
                    "Fastq_file_names"
                ]
            )
        )

        ## Get the expected SRR file sizes as well.
        localSRRFilesSizes <- SRRInfoAndDownloadFunctionDF[
            SRRInfoAndDownloadFunctionDF$GSM_IDs == iterUniqueGSMID,
            "Fastq_file_sizes"
        ]

        ## Get the sequencing type used for the given experiment. This should be
        ## a single unique value across all the SRR values annotated to the
        ## given unique GSM value, and if it is not, return a FAILED message
        ## early.
        sequencingFormat <- unique(
            SRRInfoAndDownloadFunctionDF[
                SRRInfoAndDownloadFunctionDF$GSM_IDs == iterUniqueGSMID,
                "Sequencing_type"
            ]
        )

        if (length(sequencingFormat) != 1 || is.na(sequencingFormat)) {
            return(
                "PROPER_SEQUENCING_FORMAT_NOT_RECOGNIZED"
            )
        }

        ## Create a vector noting for each of the SRR files that compose the
        ## unique GSM, if it is present and the correct size.
        SRRFilesReadyVector <- c()

        ## First let's check to ensure all the relevant files have been
        ## downloaded, and if they haven't let's try redownloading them.
        for (SRRFileIter in seq_along(localSRRFilesPaths)) {

            ## Check if the individual file's given size is known.
            if (is.na(localSRRFilesSizes[SRRFileIter])) {

                ## If the file size is not known, we cannot be sure that the
                ## given file, if it is even present, is correctly downloaded
                ## so we will return FALSE for that file.
                SRRFilesReadyVector <- c(
                    SRRFilesReadyVector,
                    FALSE
                )

            } else {

                ## Do a final check if the file is now present and of the
                ## correct size. If it is, return TRUE for the
                ## SRRFilesReadyVector, else, FALSE.
                if (
                    !file.exists(localSRRFilesPaths[SRRFileIter]) ||
                    file.info(
                        localSRRFilesPaths[SRRFileIter]
                    )[1,1] != localSRRFilesSizes[SRRFileIter]
                ) {

                    SRRFilesReadyVector <- c(
                        SRRFilesReadyVector,
                        FALSE
                    )
                } else {
                    SRRFilesReadyVector <- c(
                        SRRFilesReadyVector,
                        TRUE
                    )
                }
            }
        }

        ## Now, check that all the SRR files that comprise the given GSM are
        ## present and ready to be combined. If they are not, return
        ## a message noting the failure.
        if (!all(SRRFilesReadyVector)) {
            return(
                "ALL_SRR_FILES_NOT_PRESENT"
            )
        } else {

            ## Since all the files are present, start trying to combine them
            ## into either one or two fastq files for the given GSM, depending
            ## on if single or paired end sequencing was performed.
            if (sequencingFormat == "PAIRED") {

                ## Determine the fastq files that will be created for the GSM.
                GSMFastqPath1 <- file.path(
                    internalOutputGSMFastqDirectory,
                    paste0(
                        iterUniqueGSMID,
                        "_1.fastq.gz"
                    )
                )

                GSMFastqPath2 <- file.path(
                    internalOutputGSMFastqDirectory,
                    paste0(
                        iterUniqueGSMID,
                        "_2.fastq.gz"
                    )
                )

                ## Do a check to see if either of the GSM fastq files already
                ## exist and return a message if they do.
                if (all(
                    file.exists(GSMFastqPath1), file.exists(GSMFastqPath2)
                )) {
                    return(
                        "GSM_FASTQ_FILE_ALREADY_PRESENT"
                    )
                }

                ## There is a rare issue seen in some repos, such as that for
                ## GSM2856807 - SRR6290078, where there are only three files
                ## instead of the usual pair of 2. When this is the case,
                ## return a failed message noting the issue, and don't attempt
                ## to combine files.
                if (length(localSRRFilesPaths) %% 2 != 0) {
                    return(
                        "IMPROPER_FILE_PAIRING_DETECTED"
                    )
                }

                ## All files of a given read direction should be annotated
                ## with "_1" before the file extension, and the files of
                ## the other read direction should be annotated with "_2".
                ## Get the read_1, read_2, and other read file(s)
                filesListedRead1 <- localSRRFilesPaths[
                    substring(
                        sub('\\.fastq.*', '', basename(localSRRFilesPaths)),
                        nchar(
                            sub('\\.fastq.*', '', basename(localSRRFilesPaths))
                        ) -1,
                        nchar(
                            sub('\\.fastq.*', '', basename(localSRRFilesPaths))
                        )
                    ) == "_1"
                ]

                filesListedRead2 <- localSRRFilesPaths[
                    substring(
                        sub('\\.fastq.*', '', basename(localSRRFilesPaths)),
                        nchar(
                            sub('\\.fastq.*', '', basename(localSRRFilesPaths))
                        ) -1,
                        nchar(
                            sub('\\.fastq.*', '', basename(localSRRFilesPaths))
                        )
                    ) == "_2"
                ]

                filesListedReadOther <- localSRRFilesPaths[
                    !(substring(
                        sub('\\.fastq.*', '', basename(localSRRFilesPaths)),
                        nchar(
                            sub('\\.fastq.*', '', basename(localSRRFilesPaths))
                        ) -1,
                        nchar(
                            sub('\\.fastq.*', '', basename(localSRRFilesPaths))
                        )
                    ) %in% c("_1", "_2"))
                ]

                ## If reads without _1 or _2 before the file extensions,
                ## return a FAILED message regarding improper pairing.
                if (length(filesListedReadOther) > 0) {
                    return(
                        "IMPROPER_FILE_PAIRING_DETECTED"
                    )
                }

                ## Also do another check to ensure the same number of files
                ## are listed for both the _1 and _2.
                if (length(filesListedRead1) != length(filesListedRead2)) {
                    return(
                        "IMPROPER_FILE_PAIRING_DETECTED"
                    )
                }

                ## Create a vector noting whether each of the paired end files needed to
                ## be created (to better craft return message).
                pairedEndReturnVector <- c()

                ## There may be cases where a paired end file already exists but the
                ## other doesn't. In that case, note which has been done here.
                if (!file.exists(GSMFastqPath1)) {

                    ## If only a single SRR file is listed for each read direction
                    ## Simply rename the existing SRR files to match the GSM (since
                    ## no files need to be combined).
                    if (length(filesListedRead1) == 1) {

                        ## The SRR files are wanted, so copy them and rename the
                        ## copies after the GSM.
                        file.copy(
                            from = filesListedRead1,
                            to = GSMFastqPath1,
                            overwrite = TRUE
                        )

                    } else {

                        ## If multiple SRR files are included for each read
                        ## direction, combine them into a single fastq for each
                        ## direction (named after the GSM) using the
                        ## Rfastp::catfastq function.
                        Rfastp::catfastq(
                            output = GSMFastqPath1,
                            inputFiles = filesListedRead1
                        )
                    }

                    ## Add to the return vector.
                    pairedEndReturnVector <- c("end1")
                }

                if (!file.exists(GSMFastqPath2)) {

                    ## If only a single SRR file is listed for each read direction
                    ## Simply rename the existing SRR files to match the GSM (since
                    ## no files need to be combined).
                    if (length(filesListedRead2) == 1) {

                        ## The SRR files are wanted, so copy them and rename the
                        ## copies after the GSM.
                        file.copy(
                            from = filesListedRead2,
                            to = GSMFastqPath2,
                            overwrite = TRUE
                        )

                    } else {

                        ## If multiple SRR files are included for each read
                        ## direction, combine them into a single fastq for each
                        ## direction (named after the GSM) using the
                        ## Rfastp::catfastq function.
                        Rfastp::catfastq(
                            output = GSMFastqPath2,
                            inputFiles = filesListedRead2
                        )
                    }

                    ## Add to the return vector.
                    pairedEndReturnVector <- c(
                        pairedEndReturnVector,
                        "end2"
                    )
                }

                ## Do a final check to ensure the two fastqs refecting the
                ## two read directions were created for the given GSM.
                if (
                    file.exists(GSMFastqPath1) && file.exists(GSMFastqPath2)
                ) {
                    return("PASS")
                } else if ("end1" %in% pairedEndReturnVector) {

                    ## Make sure the first paired end fastq was created.
                    if (file.exists(GSMFastqPath1)) {
                        return("PASS_PAIRED_END_1_CREATED")
                    } else {
                        return("GSM_FASTQ_NOT_CREATED")
                    }

                } else if ("end2" %in% pairedEndReturnVector) {

                    ## Make sure the second paired end fastq was created.
                    if (file.exists(GSMFastqPath1)) {
                        return("PASS_PAIRED_END_2_CREATED")
                    } else {
                        return("GSM_FASTQ_NOT_CREATED")
                    }

                } else {
                    return("GSM_FASTQ_NOT_CREATED")
                }

            } else if (sequencingFormat == "SINGLE") {

                ## Determine the fastq file that will be created for the GSM.
                GSMFastqPath <- file.path(
                    internalOutputGSMFastqDirectory,
                    paste0(
                        iterUniqueGSMID,
                        ".fastq.gz"
                    )
                )

                ## Do a check to see if the GSM fastq file already exists and
                ## return a message it does.
                if (file.exists(GSMFastqPath)) {
                    return(
                        "GSM_FASTQ_FILE_ALREADY_PRESENT"
                    )
                }

                ## If only a single SRR file is listed, simply rename the
                ## existing SRR files to match the GSM (since no files need to
                ## be combined).
                if (length(localSRRFilesPaths) == 1) {

                    ## The SRR files are wanted, so copy them and rename the
                    ## copies after the GSM.
                    file.copy(
                        from = localSRRFilesPaths,
                        to = GSMFastqPath,
                        overwrite = FALSE
                    )
                } else {

                    ## If multiple SRR files are included, combine them into a
                    ## single fastq (named after the GSM) using the
                    ## Rfastp::catfastq function.
                    Rfastp::catfastq(
                        output = GSMFastqPath,
                        inputFiles = localSRRFilesPaths
                    )
                }

                ## Do a final check to ensure the fastq was created for the
                ## given GSM.
                if (
                    file.exists(GSMFastqPath)
                ) {
                    return("PASS")
                } else {
                    return("GSM_FASTQ_NOT_CREATED")
                }
            } else {

                ## An improper sequencingFormat value is detected, return a
                ## FAILED message reflecting this.
                return(
                    "PROPER_SEQUENCING_FORMAT_NOT_RECOGNIZED"
                )
            }
        }
    }

    ## Close the cluster list.
    ## suppressWarnings is used to silence a warning that relates to closing
    ## unused connections, which doesn't affect how the function works.
    suppressWarnings(parallel::stopCluster(combineGSMClusterList))

    ## Set the names of the combineGSMVec to be the unique GSM IDs.
    names(combineGSMVec) <- uniqueGSMIDValues

    ## Return the combineGSMVec.
    return(combineGSMVec)
}

#' Combine downloaded SRR .fastq files into .fastq files assigned per GSM ID.
#'
#' This function takes a data frame with information on SRR IDs and the state of
#' the download of their associated .fastq files from the European Nucleotide
#' Archive, such as one generated by the 'SRRFastqDownload' function, and
#' combines the .fastq files associated with the SRRs in a user-selected
#' directory into larger .fastq files associated with each GSM ID that called
#' the individual SRRs.
#'
#' @param SRRFastqDownloadDF Specify a dataframe with relevant information
#' about the SRRs and their associated .fastq files, such as one output by the
#' SRRFastqDownload function. This dataframe should contain columns with
#' the names 'GSM_IDs', 'SRR_IDs', 'Sequencing_type', and 'Fastq_file_sizes'.
#' @param localDownloadDirectory Set a path to a directory where the user has
#' downloaded the .fastq files associated with the SRRs detailed in the
#' `SRRFastqDownloadDF`. Defaults to the user's working directory.
#' @param outputGSMFastqDirectory Set a path to a directory where the user wants
#' the combined .fastq files for the GSM IDs to be output to. Defaults to the
#' user's working directory.
#' @return Combines the previously downloaded SRR .fastq files in the directory
#' specified as the `localDownloadDirectory` into .fastq files for each GSM ID.
#' These GSM IDs will be deposited in the directory specified as the
#' `outputGSMFastqDirectory`. The function will also return the data frame given
#' as the `SRRFastqDownloadDF` with an additional column named
#' "GSMCombinationStatus" noting the status of the combination of the .fastq
#' files for each GSM ID present.
#' @export
#'
#' @examplesIf interactive()
#' ## This example will attempt to combine the 6? .fastq files for SRRs
#' ## associated with three example GSM IDs using an example dataset with
#' ## results from the previous 'SRRFastqDownload' function.
#' ## NOTE: TO RUN THIS EXAMPLE, THE EXAMPLE FROM THE 'SRRFastqDownload' MUST
#' ## HAVE BEEN PREVIOUSLY RUN TO DOWNLOAD SRR .fastq FILES TO THE USER'S
#' ## WORKING DIRECTORY.
#'
#' ## Load the example dataset.
#' utils::data(
#'     "exampleDFGEODownloadR",
#'     package = "GEODownloadR"
#' )
#'
#' ## Combine the SRR .fastqs in the user's working directory into .fastqs
#' ##  associated with the original 3 GSM IDs.
#' returnValue <- combineSRRFastqsToGSMs(
#'     SRRFastqDownloadDF = exampleDFGEODownloadR
#' )

combineSRRFastqsToGSMs <- function(
    SRRFastqDownloadDF,
    localDownloadDirectory = getwd(),
    outputGSMFastqDirectory = getwd()
) {

    ## Verify the SRRFastqDownloadDF.

    ## First ensure the object that is given to the SRRFastqDownloadDF
    ## argument is a dataframe (or matrix could work).
    .isDataFrameOrMatrix(SRRFastqDownloadDF)

    ## Then, ensure all the relevant columns are present in the data frame or
    ## matrix given as the SRRFastqDownloadDF.
    .verifyColumnsPresent(
        c(
            "GSM_IDs",
            "SRR_IDs",
            "Sequencing_type",
            "Fastq_file_sizes",
            "Fastq_file_names"
        ),
        SRRFastqDownloadDF
    )

    ## Make sure the "Fastq_file_sizes" column is numeric.
    SRRFastqDownloadDF$Fastq_file_sizes <- as.numeric(
        SRRFastqDownloadDF$Fastq_file_sizes
    )

    ## Ensure that there are still valid fastq file sizes left, and issue an
    ## error to the user if not.
    if(all(is.na(SRRFastqDownloadDF$Fastq_file_sizes))) {
        stop(
            "The values in the 'Fastq_file_sizes' column of the ",
            "dataframe given as the `SRRFastqDownloadDF` argument are ",
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
    if (!dir.exists(localDownloadDirectory)) {
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

    ## Check that the outputGSMFastqDirectory is a valid path.

    ## First ensure the outputGSMFastqDirectory is a character vector.
    if (!.isSingularCharacter(outputGSMFastqDirectory)) {
        stop(
            "The `outputGSMFastqDirectory` argument must be a single character ",
            "element noting a path where the .fastq files associated with the ",
            "SRRs in the `identifySRRsFromGSMIDsDF` dataframe should be ",
            "downloaded to."
        )
    }

    ## Next ensure that the path exists.
    if (!dir.exists(outputGSMFastqDirectory)) {
        stop(
            "The path given as the `outputGSMFastqDirectory` argument doesn't ",
            "seem to exist. Please give a valid path to download .fastq files ",
            "associated with the SRRs in the `identifySRRsFromGSMIDsDF` ",
            "dataframe to."
        )
    }

    ## Finally, make sure there is no trailing "/" in the path, and remove it
    ## if there is.
    outputGSMFastqDirectory <- .trailingFSlashRemover(outputGSMFastqDirectory)

    ## Run the internal .combineSRRToGSMFunction to combine the SRR .fastqs into
    ## those for the GSM IDs themselves. Save the output vector as a new column
    ## in the `SRRFastqDownloadDF` data frame.
    GSMCombinationStatusVec <- .combineSRRToGSMFunction(
        SRRInfoAndDownloadFunctionDF = SRRFastqDownloadDF,
        internalLocalDownloadDirectory = localDownloadDirectory,
        internalOutputGSMFastqDirectory = outputGSMFastqDirectory
    )

    SRRFastqDownloadDF$GSMCombinationStatus <- GSMCombinationStatusVec[
        SRRFastqDownloadDF$GSM_IDs
    ]

    ## Return the data frame with the combination status to the user.
    return(SRRFastqDownloadDF)
}
