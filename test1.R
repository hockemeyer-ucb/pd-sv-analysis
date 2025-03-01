get_random_location_within_sub_subject <- function(subject, nItems = 25, seed = NA) {
  ## Randomly select 5 subjects with items in them (not all tx have introns or UTR)
  if (!is.na(seed)) set.seed(seed)
  targets <- sample(seq_along(subject)[elementNROWS(subject) > 0], nItems)

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
}

test_locateVariants <- function(txdb, subjectFunc) {
  ## the genes() and locateVariants() throws warnings and messages, just drop them globally for this function
  suppressWarnings({
    suppressMessages({
      subject <- subjectFunc(txdb)
      ## Get the location of all genes
      genes <- genes(txdb)

      ## Get some random location within sub items of the GRangesList subject object
      randomLocations <- get_random_location_within_sub_subject(subject)

      ## Run locateVariants on the random location
      varLocs <- locateVariants(randomLocations, txdb, AllVariants())

      ## test to see if the return location fall within the subjects
      f <- !is.na(varLocs$GENEID) # Remove entries with NA as GENEIDs
      randomLocations[varLocs$QUERYID[f]] %within% genes[varLocs$GENEID[f]]
    })
  })
}

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(VariantAnnotation)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- keepStandardChromosomes(txdb)

## Test the different classes of item searched by locateVariants

## Test locateVarians within 5' UTR subsets
res <- test_locateVariants(txdb, fiveUTRsByTranscript)
all(res)

## Test locateVarians within 3' UTR subsets
res <- test_locateVariants(txdb, threeUTRsByTranscript)
all(res)

## Test locateVarians within cds subsets
res <- test_locateVariants(txdb, cdsBy)
all(res)

## Test locateVarians within introns subsets
res <- test_locateVariants(txdb, intronsByTranscript)
all(res)
