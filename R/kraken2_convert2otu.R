#' A function for converting Kraken2 outputs
#' @name kraken2_convert2otu
#' @description This function allows you to convert Kraken2 into Ampvis2 compatible format
#' @details something
#' @title Kraken2
#' @param input Something
#' @keywords Kraken2 Ampvis2
#' @export
#' @examples
#' kraken2_convert2otu(kraken2_input)

library(dplyr)
library(stringr)

kraken_fillempty=function(data){
  maxl=max(length(data))
  for (i in 7:(maxl-1)){
    while(length(ind <- which(data[i] == "")) > 0){
      data[ind,i] <- data[ind-1,i]
    }
  }
  return(data)
}

remove1s = function(data) {
  m1 = c()
  for(i in 1:nrow(data)) {
    if(data[i,4] == data[max(nrow(data)),4]){
      m1[i] = data[i,3]
    }
    else if(data[i+1,4] == "S1" || data[i+1,4] == "S2" || data[i+1,4] == "G1" || data[i+1,4] == "G2" || data[i+1,4] == "F1" || data[i+1,4] == "F2" || data[i+1,4] == "O1" || data[i+1,4] == "O2" || data[i+1,4] == "C1" || data[i+1,4] == "C2" || data[i+1,4] == "P1" || data[i+1,4] == "P2" || data[i+1,4] == "K1" || data[i+1,4] == "K2" || data[i+1,4] == "D1" || data[i+1,4] == "D2" || data[i+1,4] == "R1" || data[i+1,4] == "R2") {
      m1[i] = data[i+1,3] + data[i,3]
    }
    else if (data[i+1,4] != "S1" || data[i+1,4] != "S2" || data[i+1,4] != "G1" || data[i+1,4] != "G2" || data[i+1,4] != "F1" || data[i+1,4] != "F2" || data[i+1,4] != "O1" || data[i+1,4] != "O2" || data[i+1,4] != "C1" || data[i+1,4] != "C2" || data[i+1,4] != "P1" || data[i+1,4] != "P2" || data[i+1,4] != "K1" || data[i+1,4] != "K2" || data[i+1,4] != "D1" || data[i+1,4] != "D2" || data[i+1,4] != "R1" || data[i+1,4] != "R2") {
      m1[i] = data[i,3]
    }
    else if(data[i,4] == data[min(nrow(data)),4]){
      m1[i] = data[i,3]
    }
    else{
      NULL
    }
  }
  data$V3 = m1
  return(data)
}

kraken_convertotu=function(input){
  output=input%>%
    mutate(OTU=paste0("OTU",1:nrow(input),"_kraken"))%>%
    select(c(14,3,7:13)) %>%
    filter(V3 > 0) %>%
    filter(!is.na(Kingdom))
  return(output)
}

kraken_convert = function(data){
  taxtab = unique(data$V4)
  taxrm=c("1","2","3","4","5","6","7","8","9","K")
  taxtab2 = taxtab[grepl(paste0(taxrm,collapse="|"),taxtab)]
  data = remove1s(data)
  for (i in taxtab2){
    data[grepl(i,data[,4]),4]=NA
  }
  tax=c("D","P","C","O","F","G","S")
  taxcol=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  filter(data,is.na(V4)==F)
  for (i in 1:length(tax)){
    colname=taxcol[i]
    data=data%>%
      mutate("{colname}":=as.character(V6))
    if(i==length(tax)){
      data[!grepl(tax[i],data[,4]),(i+6)]=NA
    }
    else{
      data[!grepl(tax[i],data[,4]),(i+6)]=""
    }
  }
  
  data=data%>%
    filter(is.na(V4)==F)
  for (i in 1:(length(tax)-1)){
    data[grepl(tax[i],data[,4]),(i+7):13]=NA
  }
  data[grepl("U|R",data[,4]),7:13]=NA
  data = kraken_fillempty(data)
  data = kraken_convertotu(data)
  return(data)
}