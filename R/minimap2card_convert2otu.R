#' A function for converting minimap2card outputs
#' @name minimap2card_convert2otu
#' @description This function allows you to convert minimap2card into Ampvis2 compatible format
#' @details something
#' @title minimap2card
#' @param input Somethin
#' @keywords minimap2 CARD Ampvis2
#' @export
#' @examples
#' minimap2card_convert2otu(minimap2card_input)




library(ggplot2)
library(pafr)
library(dplyr)

get_aroparentid=function(parentobj){
  parentid=gsub(".*/ontology/(.+?)'>.*","\\1",parentobj)
  return(parentid)
}
get_parentname=function(parentobj){
  parentname=gsub(".*<a href='/ontology/.*'>(.+?)</a>.*","\\1",parentobj)
  return(parentname)
}
get_urldata=function(aroid){
  url=read.delim(paste0("https://card.mcmaster.ca/ontology/",aroid),stringsAsFactors = F)
  urldata=url[grepl("<tr><td>",url[,1]),1]
  return(urldata)
}
get_parentobj=function(urldata){
  parentobj=gsub(".*<div id=show-parents class=collapse>(.+?)</div>.*","\\1",urldata)
  parentobj=unlist(strsplit(parentobj,"<br>"))
  if(sum(ifelse(grepl("\\[AMR Gene Family\\]",parentobj),T,F))!=0){
    parentobj=parentobj[grepl("\\[AMR Gene Family\\]",parentobj)==T]
    parentobj=parentobj[1]
    return(parentobj)
  } else if(sum(ifelse(grepl("\\[Efflux Component\\]",parentobj),T,F))!=0){
    parentobj=parentobj[grepl("\\[Efflux Component\\]",parentobj)==T]
    parentobj=parentobj[1]
    return(parentobj)
  } else if(sum(ifelse(grepl("\\[Efflux Regulator\\]",parentobj),T,F))!=0){
    parentobj=parentobj[grepl("\\[Efflux Regulator\\]",parentobj)==T]
    parentobj=parentobj[1]
    return(parentobj)
  } else{
    parentobj=parentobj[grepl("\\[.*\\]",parentobj)==F]
    parentobj=parentobj[grepl("participates_in",parentobj)==F]
    parentobj=parentobj[grepl("part_of",parentobj)==F]
    parentobj=parentobj[grepl("derives_from",parentobj)==F]
    return(parentobj)
  }
  return(parentobj)
}
get_parentterm=function(urldata){
  parentterm=c()
  for(i in 1:5){
    if(i==1){
      parentobj=gsub(".*<td>AMR Gene Family</td><td>(.+?)</td>.*","\\1",urldata)
      print(parentobj)
      parentname=get_parentname(parentobj)
      print(parentname)
      parentid=get_aroparentid(parentobj)
      print(parentid)
      urldata=get_urldata(parentid)
    } else{
      if(grepl("<div id=show-parents class=collapse>",urldata)==T && parentid!=36006){ #&& parentid!=36005 && parentid!=36003){
        parentobj=get_parentobj(urldata)
        print(parentobj)
        parentname=get_parentname(parentobj)
        print(parentname)
        parentid=get_aroparentid(parentobj)
        print(parentid)
        urldata=get_urldata(parentid)
      } else{
        parentname=NA
      }
    }
    parentterm[i]=parentname
  }
  return(parentterm)
}

arg_convert_minimap_card=function(input,name){
  print(name)
  card_ontology=read.delim("/srv/MA/Projects/microflora_danica/analysis/projects/MFD_seges/src/huy/aro.tsv",stringsAsFactors = F)
  data=input
  data$tname=gsub(".*ARO_Name:","",data$tname)
  output=data.frame(table(data$tname))
  output=cbind(output,output[,1])
  output$Var1=gsub(".ARO:.*","",output$Var1)
  output[,3]=gsub(".*ARO:","",output[,3])
  output[,3]=gsub(".Detection.*$","",output[,3])
  output[,3]=gsub("^[a-zA-Z]",NA,output[,3])
  parentmatrix=matrix(,nrow=lengths(output),ncol=5)
  for (i in 1:lengths(output)){
    if (is.na(output[i,3])==F){
      cardid=card_ontology[grepl(paste0("ARO:",output[i,3]),card_ontology$Accession),4]
      url=read.delim(paste0("https://card.mcmaster.ca/ontology/",cardid),stringsAsFactors = F)
      urldata=url[grepl("<tr><td>",url[,1]),1]
      parentmatrix[i,1:5]=get_parentterm(urldata)
    } else{
      parentmatrix[i,1:5]=c(NA,NA,NA,NA,NA)
    }
  }
  output=cbind(output,parentmatrix)
  colnames(output)=c("Gene",name,"ARO_ID","Parent1","Parent2","Parent3","Parent4","Parent5")
  return(output)
}
arg_reorder=function(input){
  output=data.frame(lapply(input,as.character),stringsAsFactors = F)
  for(i in 4:8){
    output[grepl("process or component of antibiotic biology or chemistry",output[,i]),i]=NA
    output[grepl("determinant of antibiotic resistance",output[,i]),i]=NA
  }
  output=output%>%
    select(c(3:2,8:4,1))
  
  maxl=length(output)
  for(i in (maxl-1):3){
    output[is.na(output[,i]),(i):(maxl-1)]=output[is.na(output[,i]),(i+1):(maxl)]
  }
  output=output%>%
    mutate(Species=Gene)%>%
    filter(is.na(ARO_ID)==F)%>%
    mutate(OTU=paste0("OTU_",ARO_ID,"_card"))%>%
    select(c(10,2:9))
  output[,2]=as.numeric(output[,2])  
  colnames(output)[3:9]=c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  return(output)
}