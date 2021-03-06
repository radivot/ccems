`mkg` <-
    function(strct,TCC=TRUE,activity=FALSE,free=FALSE)  
{ 
  lets=c(LETTERS,letters)
  nums=paste(0:9)
#    browser()
  mkgC <- function(g) {
#      setwd(paste(g$wDir,"models",sep="/"))
    if (sum(dir()=="models")==0) system("mkdir models")
    setwd("models")
#		print("in mkgC")
    ofile=file(fn<-paste(g$id,".c",sep=""),"wt")
#    cat(sprintf("/* compile with R CMD SHLIB %s.c or (on Windows) R CMD SHLIB %s.c */\n",g$id,g$id),file=ofile)
    cat(sprintf("/* compile with R CMD SHLIB %s.c */\n",g$id),file=ofile)
    cat("#include <R.h> /* gives F77_CALL through R_ext/RS.h */\n\n",file=ofile)
    pnames=names(g$parmsTCC)
    cat(sprintf("static double parms[%d];\n",length(pnames)),file=ofile)
    for (i in 1:length(pnames))
      cat(sprintf("#define %25s  parms[%d]\n",pnames[i],i-1),file=ofile)
    cat(sprintf("void %s(void (* odeparms)(int *, double *))\n",id),file=ofile)
    cat(sprintf("{ int N=%d;\n",length(pnames)),file=ofile)   
    cat(" odeparms(&N, parms);",file=ofile); cat("}\n",file=ofile)
    cat("\nvoid myderivs(int *neq, double *curtime, double *statevec, double *ODERHSvec)\n",file=ofile)
    cat("{",file=ofile); 
#		print(g$specieS)
    cat(sprintf(" double %s;\n",paste(g$specieS,sep='',collapse=',')),file=ofile)
#		print(g$atomS)
    if (!g$free) for (i in 1:g$nAtomS) cat(sprintf("  %s = statevec[%d];",g$atomS[i],i-1),file=ofile) else {
      for (i in 1:1) cat(sprintf("  %s = statevec[%d];",g$atomS[i],i-1),file=ofile)
      for (i in 2:g$nAtomS) cat(sprintf("  %s = %s;",g$atomS[i],pnames[i]),file=ofile)
    }
    print(g$Z)
    cat("\n",file=ofile)
    for (Zj in g$Z) cat(sprintf("%10s = %s/Kj_%s;\n",Zj,paste(g$reactantS[[Zj]],collapse="*"),Zj),file=ofile)
    for (i in 1:ifelse(g$free,1,g$nAtomS))
    {
      if (i==1)      cat(sprintf("ODERHSvec[%d] = p*%sT ",i-1,g$atomS[i]),file=ofile) else
        cat(sprintf("ODERHSvec[%d] =     %sT ",i-1,g$atomS[i]),file=ofile)
      for (j in 1:g$nSpecieS) 
        if (j<=g$nAtomS) 
        {if (g$W[g$specieS[j],g$atomS[i]]==1) 
            cat(sprintf("-%s",g$specieS[j]),file=ofile)}
        else         
        if (g$W[g$specieS[j],g$atomS[i]]==0)
          cat("              ",file=ofile)  else
        if (g$W[g$specieS[j],g$atomS[i]]==1)
          cat(sprintf(" - %s  ",g$specieS[j]),file=ofile) else
          cat(sprintf(" - %d*%s",g$W[g$specieS[j],g$atomS[i]],g$specieS[j]),file=ofile) 
      cat(";\n",file=ofile)
    }
    cat("}\n",file=ofile)
    close(ofile)
# note: cpu time savings of putting in blanks for multiplications by 0 and 1 was neglibible => something else is making this 2x slower than 
# hard coding without masks.  
# system(sprintf("R CMD SHLIB %s.c\n",g$id))
 system(sprintf("%s CMD SHLIB %s.c\n",file.path(R.home("bin"), "R"), g$id))
#    if (.Platform$OS.type=="windows") system(sprintf("Rcmd SHLIB %s.c\n",g$id)) else  system(sprintf("R CMD SHLIB %s.c\n",g$id))
#    system(sprintf("mv %s%s  bin\n",id,.Platform$dynlib.ext))
    unlink(sprintf("%s.o",g$id))
    unlink(sprintf("%s.d",g$id))
    g$code=readLines(fn)
#    setwd(g$wDir)
    setwd("..")
    g
  }
  
  mkZ <- function(g) {  # Adds the R function fback to the g object. 
    strn="fback<-function(atoms, parmsTCC) {\n"
    for (i in 1:g$nAtomS)
      strn=paste(strn,sprintf("%10s = atoms[%d];\n",g$atomS[i],i),sep="")
    for (Zj in g$Z)
      strn=paste(strn,sprintf("%10s = %s/parmsTCC[\"Kj_%s\"];\n",Zj,paste(g$reactantS[[Zj]],collapse="*"),Zj),sep="")
    strn=paste(strn,"c(",sep="")
    for (i in 1:g$nSpecieS)
      strn=paste(strn,sprintf("%s=%s%s",g$specieS[i],g$specieS[i],ifelse(i<length(g$specieS),",",")\n}\n\n")),sep="")
    g$fback=eval(parse(text=strn))
    g
  }
  
  mkRP <- function(g) {
    strn="frp<-function(parmsTCC,kis) {\n"
    for (i in 2:g$nAtomS)
      strn=paste(strn,sprintf("%10s = parmsTCC[\"%sT\"]\n",g$atomS[i],g$atomS[i]),sep="")
    for (i in 1:g$nZ)
      {Stothen=paste(rep("S",g$W[g$Z[i],"S"]),collapse="*")
      strn=paste(strn,sprintf("%10s = %s/parmsTCC[\"Kj_%s\"];\n",g$Z[i],Stothen,g$Z[i]),sep="")
    }
    strn=paste(strn,"denom=1",sep="")
    for (i in 1:g$nZ)
      strn=paste(strn,sprintf("+%s%s",g$Z[i],ifelse(i<g$nZ,"","\n")),sep="")
    strn=paste(strn,"num=(1/",max(g$W[,"S"]),")*(",sep="")
    for (i in 1:nZ) # 
      strn=paste(strn,sprintf("%d*kis[\"k%s\"]*%s%s",g$W[g$Z[i],"S"],g$Z[i],g$Z[i],ifelse(i<g$nZ,"+",")\n")),sep="")
    strn=paste(strn,"EY=num/denom\nEY\n}\n\n",sep="")
    g$frp=eval(parse(text=strn))
    g
  }

  mkPoly <- function(g) {
    strn="fpoly<-function(parmsTCC) {\n"
    strn=paste(strn,sprintf("%10sT = parmsTCC[\"%sT\"]\n",g$atomS[1],g$atomS[1]),sep="")
    for (i in 2:g$nAtomS)
      strn=paste(strn,sprintf("%10s = parmsTCC[\"%sF\"]\n",g$atomS[i],g$atomS[i]),sep="")
    strn=paste(strn,sprintf("%10s = 1\n",g$hubChar),sep="")
    for (i in 1:g$nZ)
    {if (g$W[g$Z[i],2]>0) Stothen=paste(rep(g$atomS[2],g$W[g$Z[i],2]),collapse="*") else Stothen="1"
      strn=paste(strn,sprintf("%10s = %s/parmsTCC[\"Kj_%s\"];\n",g$Z[i],Stothen,g$Z[i]),sep="")
    }
    uniW=unique(g$W[,1])
    newW=g$W[-2,1,drop=FALSE]
    rnms=row.names(newW)
    pows=1:max(uniW)
#    print(rnms)
#    print(newW)
    strn=paste(strn,"Rpoly=c(",g$hubChar,"T,",sep="")
    for (j in pows) { print(j)
      str<-paste(rnms[which(newW==j)],collapse="+")
      if (str=="") str="0"
      strn=paste(strn,str,sep=ifelse((j<=max(uniW))&(j>1),",",""))
    }
    strn=paste(strn,")*c(-1,", 1,":",max(uniW) ,")\n",sep="")
    strn=paste(strn,"Rpoly=as.polynom(Rpoly)\n",sep="")
    strn=paste(strn,"rts=solve(Rpoly)\n",sep="")
#    strn=paste(strn,"rts=polyroot(Rpoly)\n",sep="")
#    strn=paste(strn,"print(rts)\n",sep="")
    strn=paste(strn,"out<-Re(rts[(Mod(Im(rts))<1e-6)&(Re(rts)>0)])\n",sep="")
#    strn=paste(strn,"out<-Re(rts[Im(rts)==0])\n",sep="")
#    strn=paste(strn,"print(out)\n",sep="")
    strn=paste(strn,"out\n}\n\n",sep="")
    g$fpoly=eval(parse(text=strn))
    g
  }
  
  
  
  conv2long<-function(react) {
    tmp=""
    i=1
    while (i < length(react)) {
      if(react[i]%in%lets) curLet=react[i]
      i=i+1
      tmpn=""
      while (react[i]%in%nums) {tmpn=paste(tmpn,react[i],sep=""); i=i+1}
      tmp=paste(tmp,paste(rep(curLet,as.numeric(tmpn)),collapse=""),sep="")
    }
    tmp
  }
  
  getnums<-function(react) {
    vals=NULL
    nms=NULL
    tmp=""
    i=1
    while (i < length(react)) {
      if(react[i]%in%lets) nms=c(nms,react[i])
      i=i+1
      tmpn=""
      while (react[i]%in%nums) {tmpn=paste(tmpn,react[i],sep=""); i=i+1}
      vals=c(vals,as.numeric(tmpn))
    }
    names(vals)<-nms
    vals
  }
  
  oneLess<-function(react){
    i = length(react)
    tmpn=""
    while (react[i]%in%nums) {tmpn=paste(react[i],tmpn,sep=""); i=i-1}
    strn=paste(react,collapse="")
    basK=substr(strn, 1, i)
    basn=as.numeric(tmpn)-1
    endK=substr(strn, i,i)
    paste(basK,basn,"_",endK,sep="")
  }
  
  mapStrct<-function(g) {
    # this function maps the topology in strct to a set of easier to use lists and vectors
    mylets=c("D","F","L","M","N","Q","V","Y",paste(2:9),"e","f","m","n","q","v","w","y") # leave I and J for inf and free as before
    jj=length(g$hds) # jj points to nodes
    iThread=1 # count up all the blocks and index them in order
    threads=list(NULL)
    threadsWithinSites=list(NULL)
    nodesWithinSites=list(NULL)
    usedLets=character(0)
    nSites=length(g$strct$sites)
    dfThreads=data.frame(NULL)
    for (iSite in 1:nSites) {
      currSite=g$strct$sites[[iSite]]
      threadsWithinSites[[iSite]]=as.numeric(NULL)
      nodesWithinSites[[iSite]]=as.numeric(NULL)
      for (iOligo in 1:length(currSite)) {
        currOligo=currSite[[iOligo]]
        threadSize=length(currOligo)
        nodes=(jj+1):(jj+threadSize)
        names(nodes)<-currOligo
        threads[[iThread]]=list(site=iSite,let=mylets[iThread],nodes=nodes)
        usedLets=c(usedLets,mylets[iThread])
        threadsWithinSites[[iSite]]=c(threadsWithinSites[[iSite]],iThread)
        nodesWithinSites[[iSite]]=c(nodesWithinSites[[iSite]],nodes)
        iThread=iThread+1
        jj=jj+threadSize
      } # iOligo loop through oligos, within sites, i.e. things that can be equal
      dfThreads=rbind(dfThreads,threadsWithinSites[[iSite]])
#    iSite=iSite+1
    } # currSite loop over things that cannot be equal
    nThreads=iThread-1
    names(dfThreads)<-names(g$strct$sites[[1]])
    rownames(dfThreads)<-names(g$strct$sites)
    names(threads)<-paste("t",1:nThreads,sep="")
    g$singleThread=(length(g$hds)==0)
    g$dfThreads=dfThreads
    g$threads=threads
    g$threadsWithinSites=threadsWithinSites
    g$nodesWithinSites=nodesWithinSites
    if(g$singleThread) usedLets= mylets[g$nodesWithinSites[[1]]]
    g$usedLets=usedLets
    g$nThreads=nThreads
    g$mylets=mylets
    g$nSites=nSites
    g
  }
  
  #  ###########  function definitions above
  
   
  
  if (length(strct$heads)==length(strct$sites[[1]])) {
    # check to see if first head is root and if so, eliminate it
    head0=strsplit(strct$heads[1],NULL)[[1]];
#    print(head0)
    head0N=getnums(head0)
#    print(head0N)
    if ((head0N[1]==1)&(head0N[2]==0)) strct$heads=strct$heads[-1]
  }
  Z=unlist(strct)
  nZ=length(Z);
  reacts=strsplit(Z,NULL);
  hubChar=reacts[[1]][1]
#  cat("\nhubChar is ",hubChar,"\n")
  names(reacts)<-Z
  atomS=union(hubChar,setdiff(unlist(reacts),paste(0:9,sep=""))) # union makes sure central/hub protein is first
  nAtomS=length(atomS)
#  id=paste(atomS,ifelse(free,"F",""),sep="")
  id=paste(paste(atomS,collapse=""),ifelse(free,"F",""),sep="")
#	print(id)
  reactantS=strsplit(sapply(reacts,conv2long),NULL)
  specieS=c(atomS,Z)
  nSpecieS=length(specieS)
  names(specieS)<-specieS
#	print(specieS)  
  reactaNts=sapply(reacts,getnums)
#	print(reactaNts)
  eye=diag(nAtomS)
  rownames(eye)=atomS
  colnames(eye)=atomS
  W=rbind(eye,t(reactaNts))
  W=as.data.frame(W)
#	print(W)
#	print(class(W))
  hdS=strct$heads
#	print(hdS)
  tmp=NULL
  for (i in 1:length(hdS)) tmp=c(tmp,which(Z==hdS[i])) # convert node names to indices
  hds=tmp
#  print(hds)
  KdS=sapply(reacts,oneLess) # binary Kds 
  KdS[hds]=Z[hds] #of non-head nodes
  initialStateTCC=rep(0,nAtomS) 
  names(initialStateTCC)<-atomS
  if (free) initialStateTCC=initialStateTCC[1]
#  parmsTCC=c(rep(1,nAtomS),1,rep(100,nZ))  # this won't work because spurs matrix ICs overwrite this
  parmsTCC=c(rep(1,nAtomS),1,rep(1,nZ))
  pnamesTCC=c(paste(atomS[1],"T",sep=""),paste(atomS[-1],ifelse(free,"F","T"),sep=""),"p",paste("Kj",Z,sep="_"))
  names(parmsTCC)<-pnamesTCC
  print(parmsTCC)
  gObj=list(id=id,
      hubChar=hubChar, 
      Z=Z,nZ=nZ,
      atomS=atomS, nAtomS=nAtomS,
      specieS=specieS, nSpecieS=nSpecieS,
      reactantS=reactantS,
      strct=strct,
      W=W,
      KdS=KdS,
      KS=KdS,
      kS=paste("k",Z,sep=""),
      hdS=hdS,hds=hds,
#			wDir=sub("C:","",getwd()),
      TCC=TCC,
      activity=activity,
      free=free,
      sstime = 1e6,rtol=1e-5,atol=1e-7,
      parmsTCC=parmsTCC,
      initialStateTCC=initialStateTCC
  ); 
  gObj=mapStrct(gObj)  # this was previously done in mkGrids
  gObj=mkZ(gObj)
# library(odesolve); 
#  gObj=mkgC(gObj)
  if (TCC&!free) { gObj=mkgC(gObj); testgC(gObj)} else {
    if ((activity)&(gObj$singleThread)) gObj=mkRP(gObj)
  }
  if (free) { library(PolynomF); gObj=mkPoly(gObj)}
#  if (free) { library(polynom); gObj=mkPoly(gObj)}
  gObj
}


