#' A function for converting DeepARG outputs
#' @name deeparg_convert2otu
#' @description This function allows you to convert DeepARG into Ampvis2 compatible format
#' @details something
#' @title DeepARG
#' @param input Something
#' @keywords DeepARG Ampvis2
#' @export
#' @examples
#' deeparg_convert2otu(deeparg_input)


library(dplyr)

convert_deeparg=function(input,name){
  print(nrow(input))
  input2=input%>%
    filter(input[9]>=60)
  print(nrow(input2))
  strings=strsplit(input[,6],split="\\|")
  data=input2%>%
    mutate(Gene=sapply(strings,function(x) x[5]))
  datafr=as.data.frame(table(data[13]),stringsAsFactors = F)
  class=data[c(5,13)]
  class2=class[!duplicated(class),]
  temp=datafr%>%
    left_join(class2,by=c("Var1"="Gene"))
  output=temp%>%
    mutate(OTU=paste0("OTU",1:nrow(temp),"_deeparg"))%>%
    select(c(4,2,3,1))%>%
    mutate(Var2=Var1)%>%
    mutate(Var3=Var1)%>%
    mutate(Var4=Var1)%>%
    mutate(Var5=Var1)%>%
    mutate(Var6=Var1)
  colnames(output)=c("OTU",name,"Kingdom","Phylum","Class","Order","Family","Genus","Species")
  return(output)
}