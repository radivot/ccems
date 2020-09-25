rm(list=ls(all=TRUE))  # clean up left overs from any previous run
setwd("/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/data") 
#load("TK2.rda")
d=read.table("TK2inSec.txt",header=T)
TK2=d
save(TK2,file="TK2.rda")
#TK2BBA77=read.table("TK2patientsBBA77.txt",header=T)
#save(TK2BBA77,file="TK2BBA77.rda")


