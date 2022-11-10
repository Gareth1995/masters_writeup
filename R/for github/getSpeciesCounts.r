#' Returns counts for a specific species at a specified location
#'
#' @param loc_code Numeric or character string corresponding to a CWAC site code
#' @param id Numeric or character string corresponding to a species
#'
#' @return A tibble containing the counts of a species
#' @export
#'
#' @examples
#' getSpeciesCounts(23312919, 269)
#' getSpeciesCounts("23312919", "269")

getSpeciesCounts <- function(loc_code, speciesID){
  
  cards <- CWAC::listCwacCards(loc_code)
  
  # setting up the data frame to be populated
  recs <- data.frame(matrix(nrow=dim(cards)[1], ncol = length(speciesID)))
  
  # populate the data frame
  for (i in 1:dim(cards)[1]) {
    # get records from each card in the specified site
    cd <- get_cwac_records(card=cards$Card[i])
    
    for (j in 1:length(species)) {
      # get counts of specified species in that card
      x <- as.numeric(as.character(cd$count[which(cd$spp == species[j])]))
      if (length(x)==1) recs[i,j] <- x
    }
  }
  # order data dframe by date
  spec_counts <- cbind(recs, cards$startDate, cards$Season)
  colnames(spec_counts) <- c(paste("X",species, sep=''), "startDate", "season")
  spec_counts$startDate <- lubridate::ymd(spec_counts$startDate)
  
  
  return(dplyr::arrange(spec_counts, startDate))
}

getSpeciesCounts(26352535, 269)
