# Retrieve general statistics from input Census CSV file (with labels as valies)

parseCensusCSV <- function(inputFile, outputFile='') {
    print('Loading input file')
    censusTable <- read.csv(inputFile, header=TRUE)

    print('Extracing statistics')
    listOfStreets <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L')
    nEntries <- nrow(censusTable[,])
    validEntries <- (!is.na(censusTable$weightnerds) & censusTable[,'address_letter'] %in% listOfStreets & !is.na(censusTable[,'address_hour']))
    nValidEntries <- nrow(censusTable[validEntries,])
    virginEntries <- (censusTable$virgin == 'yes' & validEntries)
    nVirgins <- nrow(censusTable[virginEntries,])
    fracVirgins <- nVirgins / nValidEntries
    weightTotal <- sum(censusTable$weightnerds[validEntries])
    weightVirgins <- sum(censusTable$weightnerds[virginEntries])
    fracVirginsWeights <- weightVirgins / weightTotal    

    # Print stats
    print('Statistics:')
    print(paste(' Total entries: ', nEntries))
    print(paste(' Valid entries: ', nValidEntries, '(', nValidEntries / nEntries * 100, '%)'))
    print(paste(' Fraction virgins (weighted):', fracVirgins, ' (', fracVirginsWeights, ')'))

    # Save data to a file easily readable in ROOT
    if (! outputFile == '') {
        print('Saving filtered table to file.')
        censusTable[,'virgin_bool'] <- ifelse(censusTable$virgin == 'yes', 1, 0)
        interestingDataCols <- c('id', 'address_letter', 'address_hour', 'virgin_bool', 'weightnerds')
        outTable <- censusTable[validEntries,interestingDataCols]
        write.table(outTable, file=outputFile, sep='\t', row.names=FALSE, quote=FALSE)
    }
    print('All Done.')
}
