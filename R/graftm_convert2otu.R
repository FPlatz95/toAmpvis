#' A function for converting GraftM outputs
#' @name graftm_convert2otu
#' @description This function allows you to convert GraftM into Ampvis2 compatible format
#' @details something
#' @title GraftM
#' @param input Somethin
#' @keywords GraftM Ampvis2
#' @export
#' @examples
#' graftm_convert2otu(graftm_input)

library(dplyr)
library(stringr)

graftm_inserttax=function(input,taxsplit){
  input%>%
    mutate(Kingdom=sapply(taxsplit,function(x) x[2]))%>%
    mutate(Phylum=sapply(taxsplit,function(x) x[3]))%>%
    mutate(Class=sapply(taxsplit,function(x) x[4]))%>%
    mutate(Order=sapply(taxsplit,function(x) x[5]))%>%
    mutate(Family=sapply(taxsplit,function(x) x[6]))%>%
    mutate(Genus=sapply(taxsplit,function(x) x[7]))%>%
    mutate(Species=sapply(taxsplit,function(x) x[8]))%>%
    filter(Kingdom!="Eukaryota")
}

graftm_convert2otu=function(input){
  taxsplit=strsplit(as.character(input$ConsensusLineage),";")
  data=graftm_inserttax(input,taxsplit)
  output=data%>%
    mutate(OTU=paste0("OTU",1:nrow(data),"_graftm"))%>%
    select(c(11,2,4:10))
}