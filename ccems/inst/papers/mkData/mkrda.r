rm(list=ls(all=TRUE))  # clean up left overs from previous run, though shouldn't need to.
#library(RODBC)
setwd("/users/radivot/case/active/ccems/ccems/inst/papers/mkData")  # directory where the next two files live
##ch <- odbcConnectExcel(paste(system.file(package = "ccems"),"/papers/mkData/RNR.xls",sep=""))
#ch <- odbcConnectExcel("/users/radivot/case/active/ccems/ccems/inst/papers/mkData/RNR.xls")
#RNR <- sqlFetch(ch, "RNR")
#odbcCloseAll()
RNR=read.table(file="RNR.txt",header=T)
RNR
#save(RNR,file="RNR.RData")
save(RNR,      file="/users/radivot/case/active/ccems/ccems/data/RNR.rda")  

