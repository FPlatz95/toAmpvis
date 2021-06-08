#' A function for converting SingleM outputs
#' @name singlem_convert2otu
#' @description This function allows you to convert SingleM into Ampvis2 compatible format
#' @details something
#' @title SingleM
#' @param input Something
#' @keywords SingleM Ampvis2
#' @export
#' @examples
#' singlem_convert2otu(singlem_input)

library(dplyr)
library(stringr)


inserttax = function(input){
  taxsplit = strsplit(as.character(input$Group.1),";")
  input %>%
    mutate(Kingdom=sapply(taxsplit,function(x) x[2])) %>%
    mutate(Phylum=sapply(taxsplit,function(x) x[3]))%>%
    mutate(Class=sapply(taxsplit,function(x) x[4]))%>%
    mutate(Order=sapply(taxsplit,function(x) x[5]))%>%
    mutate(Family=sapply(taxsplit,function(x) x[6]))%>%
    mutate(Genus=sapply(taxsplit,function(x) x[7]))%>%
    mutate(Species=sapply(taxsplit,function(x) x[8]))
}

singlem_convertotu=function(input){
  input %>%
    mutate(OTU=paste0("OTU",1:nrow(input),"_singlem")) %>%
    select(c(10,2,3:9)) %>%
    filter(Kingdom == c("Bacteria", "Archaea"))
}