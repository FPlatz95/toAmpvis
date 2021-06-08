#' A function for converting Kaiju outputs
#' @name kaiju_convert2otu
#' @description This function allows you to convert Kaiju into Ampvis2 compatible format
#' @details something
#' @title Kaiju
#' @param input Something
#' @keywords GraftM Ampvis2
#' @export
#' @examples
#' kaiju_convert2otu(kaiju_input)

library(dplyr)
library(stringr)

######## FUNCTIONS ############

insert_tax=function(input,taxsplit){
  input%>%
    mutate(Kingdom=sapply(taxsplit,function(x) x[1]))%>%
    mutate(Phylum=sapply(taxsplit,function(x) x[2]))%>%
    mutate(Class=sapply(taxsplit,function(x) x[3]))%>%
    mutate(Order=sapply(taxsplit,function(x) x[4]))%>%
    mutate(Family=sapply(taxsplit,function(x) x[5]))%>%
    mutate(Genus=sapply(taxsplit,function(x) x[6]))%>%
    mutate(Species=sapply(taxsplit,function(x) x[7]))
}


kaijutootutable = function(input){
  taxsplit=strsplit(as.character(input$taxon_name),";")
  data = insert_tax(input, taxsplit)
  output = data %>%
    mutate(suffix = "_kaiju") %>%
    mutate(OTU = paste0("OTU_", taxon_id, suffix)) %>%
    select(OTU,reads,Kingdom,Phylum,Class,Order,Family,Genus,Species) %>%
    filter(OTU != paste0("OTU_", "NA_","kaiju"))
}