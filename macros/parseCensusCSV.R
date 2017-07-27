# Retrieve general statistics from input Census CSV file (with labels as values)
# inputFile: input CSV file with Census results
# outputFile: output TSV file with filtered sample (see format below)
# yCatLabel: label of category that you want to extract (e.g. 'virgin')
# yCatValues: array of values for yCatLabel to assign to y-category 0 (others will be assigned to category 1)
# yCatExclude: array of values for yCatLabel to exclude from analysis
#
# Output file is a tab-separated file (easily readable in a ROOT TTree) with the following columns:
# 'id', 'address_letter', 'address_hour', yCatLabel, 'weightnerds'
#
# A filtering of events based on valid streets for analysis as well as valid entries and weight is performed.
#
parseCensusCSV <- function(inputFile, outputFile='', yCatLabel='virgin', yCat0Values=c("yes"), yCatExclude=c()) {
    print('Loading input file')
    censusTable <- read.csv(inputFile, header=TRUE)

    print('Extracing statistics')
    listOfStreets <- c('Esplanade', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L')
    nEntries <- nrow(censusTable[,])
    validEntries <- (!is.na(censusTable$weightnerds) & censusTable[,'address_letter'] %in% listOfStreets & !is.na(censusTable[,'address_hour']))
    nValidEntries <- nrow(censusTable[validEntries,])

    yCat0Entries <- (censusTable[,yCatLabel] %in% yCat0Values & validEntries)
    nYCat0 <- nrow(censusTable[yCat0Entries,])
    fracYCat0 <- nYCat0 / nValidEntries
    weightTotal <- sum(censusTable$weightnerds[validEntries])
    weightYCat0 <- sum(censusTable$weightnerds[yCat0Entries])
    fracYCat0Weights <- weightYCat0 / weightTotal    

    # Print stats
    print('Statistics:')
    print(paste(' Total entries: ', nEntries))
    print(paste(' Valid entries: ', nValidEntries, '(', nValidEntries / nEntries * 100, '%)'))
    print(paste(' Fraction of ', yCatLabel, ' in [', yCat0Values, '] excluding [', yCatExclude, '] (weighted):', fracYCat0, ' (', fracYCat0Weights, ')'))

    # Save data to a file easily readable in ROOT
    if (! outputFile == '') {
        print(paste0('Saving filtered table to file ', outputFile))
        censusTable[,yCatLabel] <- ifelse(censusTable[,yCatLabel] %in% yCat0Values, 1, 0)
        interestingDataCols <- c('id', 'address_letter', 'address_hour', yCatLabel, 'weightnerds')
        outTable <- censusTable[validEntries,interestingDataCols]
        write.table(outTable, file=outputFile, sep='\t', row.names=FALSE, quote=FALSE)
    }
    print('All Done.')
}
