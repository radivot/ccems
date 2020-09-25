library(RODBC)
setwd("/Users/radivot/case/active/ccems/ccems/inst/papers/TK2/data") 
TK2=read.table("TK2.txt",header=T)
TK2

# 1 unit= 1 nmol/min = 
# v is in units/mg
# 29 mg/umole  and 2 ng placed in 50 uL
# 1000 nmol = 29 mg => 1 mg = 1000/29 = 34.4 nmoles
# nmole/29 ug * 0.002 ug/50 uL = 0.002/29/50= 1.37e-6 nmoles/uL = 1.37 nmoles/L = 1.37 nM 
d=subset(TK2,(dT>0)&(year==1991))
d1=transform(d,E=0.00137,k=VdT/60/34.4)
d=subset(TK2,(dC>0)&(year==1991))
d2=transform(d,E=0.00137,k=VdC/60/34.4)
dd1=rbind(d1,d2)

# 1999 FEBS letter
# .25 mg/mL   46 fold dilution
# .25/46 = 5.4 ug/mL
# 5.43 ug/ml * 1 uL used in an assay = 5.43 ng enzyme in volume 50 uL
# 0.0054/29/50
# 3.75 nM
d=subset(TK2,(dT>0)&(year==1999))
d1=transform(d,E=0.00375,k=VdT/60/34.4)
d=subset(TK2,(dC>0)&(year==1999))
d2=transform(d,E=0.00375,k=VdC/60/34.4)
dd2=rbind(dd1,rbind(d1,d2))

# here one unit is defined as a nmol/min/mg
# from above 1 mg = 1000/29 = 34.4 nmoles
d=subset(TK2,(!is.na(VdT))&(year==2003)&(frstAut=="Wang"))
d1=transform(d,E=0.00375,k=VdT/60/34.4)
d=subset(TK2,(!is.na(VdC))&(year==2003)&(frstAut=="Wang"))
d2=transform(d,E=0.00375,k=VdC/60/34.4)
dd3=rbind(dd2,rbind(d1,d2))
dd3

# 50 uL reaction vessel with 1 ug/mL TK2 (1 nmole/29 ug) = 1000/29 = 34 nmoles/L, i.e. close to the lowest [S] of 50 nM
d=subset(TK2,(!is.na(VdT))&(year==2003)&(frstAut=="Barroso"))
d1=transform(d,E=0.0344,k=VdT/60/34.4)
d=subset(TK2,(!is.na(VdC))&(year==2003)&(frstAut=="Barroso"))
d2=transform(d,E=0.0344,k=VdC/60/34.4)
dd4=rbind(dd3,rbind(d1,d2))
dd4

# reduced [E] to .2 ug/mL for dT and .5 for dC, in profiles, so divide by 5 and 2 from above 
d=subset(TK2,(!is.na(VdT))&(year==2005)&(frstAut=="Barroso"))
d1=transform(d,E=0.0344/5,k=VdT/60/34.4)
d=subset(TK2,(!is.na(VdC))&(year==2005)&(frstAut=="Barroso"))
d2=transform(d,E=0.0344/2,k=VdC/60/34.4)
dd5=rbind(dd4,rbind(d1,d2))
dd5

write.table(dd5,file="TK2inSec.txt",sep="\t",row.names=F)
