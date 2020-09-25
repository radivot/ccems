library(ccems)  # loads the ccems package into R
setwd("~")
topology <- list(
    heads=c("R1X0","R2X2","R4X4","R6X6"),   # X = ATP and R = R1
    sites=list(   # s-sites are already filled only in (j>1)-mers 
        a=list(       #a-site                                                   
            m=c("R1X1"),                               # monomer
            d=c("R2X3","R2X4"),                        # dimers
            t=c("R4X5","R4X6","R4X7","R4X8"),          # tetramers
            h=c("R6X7","R6X8","R6X9","R6X10", "R6X11", "R6X12") # hexamers
        ), # tails of a-site threads are heads of h-site threads
        h=list(   # h-site
            m=c("R1X2"),                               # 1-mer  
            d=c("R2X5", "R2X6"),                       # 2-mers
            t=c("R4X9", "R4X10","R4X11", "R4X12"),     # 4-mers
            h=c("R6X13", "R6X14", "R6X15","R6X16", "R6X17", "R6X18") #6-mers
        )
    )
)
g=mkg(topology) 
dd=subset(RNR,(year==2002)&(fg==1)&(X>0),select=c(R,X,m,year)) # get data
names(dd)[1:2]=c("RT","XT")
cpus=c("localhost"=4,"compute-0-0"=4,"compute-0-1"=4,"compute-0-2"=4,
     "compute-0-3"=4, "compute-0-4"=4,"compute-0-5"=4,"compute-0-6"=4)
top10=ems(dd,g,cpusPerHost=cpus,maxTotalPs=3,ptype="SOCK",KIC=100,topN=10,transform="none") 
