#' Retrives all species recorded from a given site
#'
#' @param loc_code Numeric or character string corresponding to a CWAC site code
#'
#' @return A tibble with the species in the given site
#' @export
#'
#' @examples
#' getLocalSpecies(23312919)
#' getLocalSpecies("23312919")

getLocalSpecies <- function(loc_code){
  
  url = paste0('http://api.adu.org.za/cwac/cards/single/list?locationcode=',loc_code)
  
  myfile <- getURL(url, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)
  
  json_file <- rjson::fromJSON(myfile)
  
  # changing nulls in the Json files to NAs
  json_file <- lapply(json_file, function(x) {
    x[sapply(x, is.null)] <- NA
    unlist(x)
  })
  
  df <- as.data.frame(do.call("rbind", json_file))
  
  sdf <- stack(df)
  sdf$ind <- as.character(sdf$ind)
  
  species <- as.data.frame(sdf[which(startsWith(sdf$ind, 'list.common')),'values'])
  specIds <- as.data.frame(sdf[which(startsWith(sdf$ind, 'list.ref')),'values'])
  taxonomic <- as.data.frame(sdf[which(startsWith(sdf$ind, 'list.taxonomic')),'values'])
  
  colnames(species) <- "species"
  colnames(specIds) <- "id"
  colnames(taxonomic) <- "Scientific name"
  
  specList <- cbind(species, specIds, taxonomic)
  
  return(specList)
  
}

getLocalSpecies(26352535)
