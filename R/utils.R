## Define %dopar% from the foreach package so it can be seen and used by
## package functions.
`%dopar%` <- foreach::`%dopar%`

## Write an internal function to check if a given object is a single boolean in
## order to verify arguments of package functions. Note: the !is.matrix()
## clause prevents a 1x1 boolean from being registered as a singular object.
.isSingularBoolean <- function (x) {
    return((is.logical(x) && length(x) == 1) && !is.matrix(x))
}

## Write an internal function to check if a given object is a single character
## element in order to verify arguments of package functions. Note: the
## !is.matrix() clause prevents a 1x1 vector from being registered as a
## character still.
.isSingularCharacter <- function (x) {
    return((is.character(x) && length(x) == 1) && !is.matrix(x))
}

## Write an internal function to check if a given object is a character vector
## in order to verify arguments of package functions. Note: the !is.matrix()
## clause prevents a 1x1 vector from being registered as a character still.
.isCharacterVector <- function (x) {
    return(is.character(x) && !is.matrix(x))
}

## Write an internal function to check if a given object is a data frame or
## matrix in order to verify arguments of package functions.
.isDataFrameOrMatrix <- function (x) {
    return((is.data.frame(x) || is.matrix(x)))
}

## Write an internal function to check if a vector of GSM IDs are valid IDs,
## starting with 'GSM' and containing numeric values after that. If at least one
## isn't, stop the function and issue a message to the user, listing which of
## the given GSM IDs weren't valid.
.isValidGSMIDs <- function (x) {
    improperGSMIDIters <- unique(
        sort(
            c(
                which((substr(x, 1, 3) != "GSM")),
                which(
                    is.na(suppressWarnings(as.numeric(substr(x, 4, nchar(x)))))
                )
            )
        )
    )

    if(length(improperGSMIDIters > 0)) {
        stop(
            "The following GSMIDs are invalid: ",
            paste(
                x[improperGSMIDIters],
                collapse = ", "
            ),
            ". Please ensure these GSM IDs are properly formatted. Each ",
            "should start with 'GSM' and end with numeric values.",
            call. = FALSE
        )
    }
}

## Write an internal function, which checks if the last character in a path is
## a "/", and removes it if it is.
.trailingFSlashRemover <- function (x) {
    if (substr(x, nchar(x), nchar(x)) == "/") {
        return(substr(x, 1, (nchar(x)-1)))
    } else {
        return(x)
    }
}

## Write a function to verify columns are present in the given data frame. If
## at least one column isn't present, stop the function and issue a message to
## the user, noting which columns aren't present.
.verifyColumnsPresent <- function (columnVector, dataFrame) {

    ## Find the column names in the columnVector that aren't present in the
    ## dataframe.
    columnsNotPresent <- setdiff(columnVector, colnames(dataFrame))

    # ## If any of the specified columns aren't present, stop the function and
    # ## tell the user which columns aren't present.
    # dataFrameArgument <- deparse(substitute(paste0("`", dataFrame, "`")))

    if (length(columnsNotPresent) > 0) {
        stop(
            "The following columns were not found in the `",
            deparse(substitute(dataFrame)),
            "`: ",
            paste(
                columnsNotPresent,
                collapse = ", "
            ),
            ". Please ensure these columns are present in the data frame ",
            "given as the `",
            deparse(substitute(dataFrameArgument)),
            "` argument and try again.",
            call. = FALSE
        )
    }
}
