# Functions for CWAC data download using the CWAC API
# 11 Oct 2019

# from the CWAC website: http://cwac.adu.org.za/api.php


require(RCurl)
#library(digest)
require(rjson)


# Retrieve a full list of all the active CWAC sites registered in the Western Cape
# *******************************************************************************
get_cwac_sites <- function(province="western%20cape"){
  # Options: Northern Province, Mpumalanga, North West, Gauteng, KwaZulu-Natal, 
  # Free State, Northern Cape, Western Cape, Eastern Cape, Kenya, Angola, Tanzania, Limpopo
  
url=paste('http://api.adu.org.za/cwac/sites/list?province=',province, sep='')

myfile <- getURL(url, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

json_file <- rjson::fromJSON(myfile)

json_file <- lapply(json_file, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

return(as.data.frame(do.call("rbind", json_file)))
}


# List of all cards submitted for a CWAC site
# *******************************************************************************

get_cwac_cards <- function(site="33561832"){
  
url=paste('http://api.adu.org.za/cwac/site/cards/list?locationCode=',site, sep='')

myfile <- getURL(url, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

json_file <- rjson::fromJSON(myfile)

json_file <- lapply(json_file$cards, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

return(as.data.frame(do.call("rbind", json_file)))
}



# download full data for a single list
# *******************************************************************************

get_cwac_records <- function(card="505786"){
  
  url=paste('http://api.adu.org.za/cwac/cards/single/get?card=',card,sep='')

myfile <- getURL(url, ssl.verifyhost=FALSE, ssl.verifypeer=FALSE)

json_file <- rjson::fromJSON(myfile)

json_file <- lapply(json_file$records, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

return(as.data.frame(do.call("rbind", json_file)))
}







