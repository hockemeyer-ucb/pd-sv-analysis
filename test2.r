get_random_location_within_sub_subject <- function(subject, nItems = 25, seed = NA) {
  ## Randomly select GRangesList
  if (!is.na(seed)) set.seed(seed)
  targets <- sample(seq_along(subject)[elementNROWS(subject) > 0], nItems)

  ## For each selected items find a location within a sub-item
  if (!is.na(seed)) set.seed(seed)
  randomLocations <- do.call(c, lapply(targets, function(x) {
    ## Pick one of the items in the selected entry
    y <- sample(length(subject[x]), 1)

    ## Get the chromsome of the random sub-item
    chr <- seqnames(subject[[x]][y])

    ## Return a single random postion within the selected sub item
    set.seed(1234)
    loc <- sample(seq(start(subject[[x]][y]), end(subject[[x]][y])), 1)

    ## Return a GRanges for that random item
    GRanges(chr, IRanges(loc))
  }))
  return(randomLocations)

}


###############################################################################
###############################################################################
## Get the 5' UTR from Txdb
subject <- fiveUTRsByTranscript(txdb)

## Get random locations within 5' UTR
randomLocations <- get_random_location_within_sub_subject(subject, seed = 77458)

## Map the locations of the random GRanges to the subjects
map <- mapToTranscripts(randomLocations, subject)

## Are the mapped items within the ranges of the GRangesList
all(randomLocations[map$xHits] %within% unlist(range(subject[map$transcriptsHits])))

###############################################################################
###############################################################################

## Get the introns from Txdb
subject <- intronsByTranscript(txdb)

## Get random locations
randomLocations <- get_random_location_within_sub_subject(subject)

## Map the locations of the random GRanges to the subjects
map <- mapToTranscripts(randomLocations, subject)

## Are the mapped items within the ranges of the GRangesList
all(randomLocations[map$xHits] %within% unlist(range(subject[map$transcriptsHits])))
