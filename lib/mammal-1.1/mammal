#!/usr/bin/env Rscript
bindir <- ""
cseries.dir <- ""
##------------------------------------------------------------------------------
HclustCenters <- function(x,nclass=21,dmethod="manhattan",
                          linkage="average",hclust.type="hclust"){
  if(hclust.type == "hclust"){
    d <- dist(x,method=dmethod)
    h <- hclust(d, method=linkage)
  }else{
    h <- Rclusterpp.hclust(x,method=linkage,distance=dmethod)
  }
  centers <- tapply(x,
                    list(rep(cutree(h,nclass),ncol(x)), col(x)), mean)
  return(centers)
}
DataFrequencies <- function(seqfile,clean=TRUE,
                            bindir="~/bio/huaichun/bordor/"){
  system(paste(bindir,"mult-data -s -i ",seqfile," > tmp.fr",sep=""))
  fr <- matrix(scan("tmp.fr",quiet=TRUE),ncol=22,byrow=TRUE)
  fr <- fr[,1:20]/apply(fr[,1:20],1,sum)
  if(clean) file.remove("tmp.fr")
  return(fr)
}
OptimizeWeights <- function(Sigma){
  ## min(-d^T b + 1/2 b^T D b) with the constraints A^T b >= b_0
  ntaxa <- nrow(Sigma)
  Amat <- rep(1,ntaxa)
  Amat <- cbind(Amat, diag(rep(1,ntaxa)))
  bvec <- c(1,rep(0,ntaxa))
  l <- solve.QP(Dmat=Sigma*2, dvec=rep(0,ntaxa), Amat=Amat, bvec=bvec)
  return(list(w=l$solution,v=l$value))
}
CalculateWeightsE <- function(seqfile,normalize=FALSE,clean=TRUE,DNA=FALSE,
                              bindir="~/bio/mix-freq/comp-lik-variance/"){
  cmdline <- paste(bindir,"mammal-sigma",sep="")
  if(DNA) cmdline <- paste(cmdline,"-D")
  cmdline <- paste(cmdline,"<",seqfile,">tmp.Sigma")
  ntaxa <- scan(seqfile,quiet=TRUE,n=1)
  ## cat(cmdline,"\n")
  system(cmdline)
  Sigma <- matrix(scan("tmp.Sigma",quiet=TRUE),ncol=ntaxa)
  l <- OptimizeWeights(Sigma)
  if(clean) file.remove("tmp.Sigma")
  if(normalize){
    l$w[l$w<0] <- 0
    l$w <- ntaxa*l$w/sum(l$w)
  }
  return(l)
}
CreateIQFreqfile <- function(fr,fname){
  cat("#nexus\nbegin models;\n",file=fname)
  for(ic in 1:dim(fr)[1]){
    cat("frequency ESclass",ic," = ",sep="",file=fname,append=TRUE)
    cat(fr[ic,],";\n",file=fname,append=TRUE)
  }
  cat("\n    frequency ESmodel = FMIX{",sep="",file=fname,append=TRUE)
  for(ic in 1:(dim(fr)[1]-1))
    cat("ESclass",ic,",",sep="",file=fname,append=TRUE)
  cat("ESclass",dim(fr)[1],"};\nend;\n",sep="",file=fname,append=TRUE)
}
##------------------------------------------------------------------------------
## MAMMaL Version 1.1
require("quadprog",quietly=TRUE)

## External programs
dgpe <- paste(bindir,"dgpe",sep="")
mult.data <- paste(bindir,"mult-data",sep="")

FnameE <- function(fname){
  if(!file.exists(fname)) stop(paste(fname,"does not exist"))
  return(fname)
}

## Process Command Line
args <- commandArgs()
iarg <- length(args)

iqtreefile <- nclass <- seqfile <- treefile <- NULL
q <- 0.75 # quantile for rate estimation
iq.freqs <- lwt.needed <- TRUE
start.frs.type <- "hclust"
plusF <- hclust.set <- FALSE
C <- 5 # penalty parameter

while(iarg>=6){
  not.an.option <- TRUE
  opt <- args[iarg]
  if(substring(opt,1,1)!='-' && iarg==6) stop(paste(opt,"is not an option"))
  val <- NULL
  if(substring(opt,1,1)!='-'){
    val <- opt
    iarg <- iarg - 1
    opt <- args[iarg]
  }
  if(opt=="-s"){
    if(is.null(val)) stop("seqfile not specified in -s seqfile")
    seqfile <- FnameE(val); not.an.option <- FALSE
  }
  if(opt=="-t"){
    if(is.null(val)) stop("treefile not specified in -t treefile")
    treefile <- FnameE(val); not.an.option <- FALSE
  }
  if(opt=="-c"){
    if(is.null(val)) stop("number of classes not specified in -c nclasses")
    nclass <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-h"){
    start.frs.type <- "hclust"; hclust.set <- TRUE; not.an.option <- FALSE
  }
  if(opt=="-m"){
    iq.freqs <- FALSE; not.an.option <- FALSE
  }
  if(opt=="-l"){
    lwt.needed <- FALSE; not.an.option <- FALSE
  }
  if(opt=="-q"){
    if(is.null(val)) stop("quantile not specified in -q quantile")
    q <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-C"){
    if(is.null(val)) stop("penalty not specified in -C penalty")
    C <- as.numeric(val); not.an.option <- FALSE
  }
  if(opt=="-d"){
    plusF <- TRUE; not.an.option <- FALSE
  }
  if(not.an.option) stop(paste(opt,"is not an option"))
  iarg <- iarg - 1 
}

if(is.null(seqfile)) stop("sequence file needed: -s seqfile")
if(is.null(treefile)) stop("tree file needed: -t treefile")
if(is.null(nclass)) stop("must specify number of classes: -c nclasses")

if(nclass %in% seq(10,60,10) & !hclust.set) start.frs.type <- "c-series"

## Extract high rate sites
system(paste(dgpe,"-i",seqfile,"-t",treefile,"-o tmp.high-rate-seqfile","-q",q,
             " > tmp.out"))

## Likelihood weights
ntaxa <- as.numeric(scan("tmp.high-rate-seqfile",n=1,quiet=TRUE))
if(lwt.needed){
  lwt <- CalculateWeightsE(seqfile="tmp.high-rate-seqfile",bindir=bindir,
                           clean=FALSE)$w
}else{
  lwt <- rep(1,ntaxa)
}
if(sum(lwt<0)>0) lwt[lwt<0] <- 0
lwt <- ntaxa*lwt/sum(lwt)
write(format(lwt,sci=TRUE,digits=16),ncol=1,file="tmp.lwt")

## Starting frequencies
if(start.frs.type=="hclust"){
  system(paste(mult.data," -s -i",seqfile,"> tmp.fr"))
  fr <- matrix(scan("tmp.fr",quiet=TRUE),ncol=22,byrow=TRUE)
  fr <- fr[,1:20]/apply(fr[,1:20],1,sum)
  frs <- HclustCenters(DataFrequencies(seqfile,clean=TRUE,bindir=bindir),
                       hclust.type=start.frs.type,
                       nclass=nclass,dmethod="manhattan")
  ## print(frs)
}
if(start.frs.type=="c-series"){
  frs <- scan(paste(cseries.dir,"C",nclass,".aafreq.dat",sep=""),quiet=TRUE)
  frs <- matrix(frs,ncol=20,byrow=TRUE)
}
write(t(frs),file="tmp.frs",ncol=20)
if(plusF){
  system(paste(bindir,"charfreq 20 < ",seqfile," >> tmp.frs",sep=""))
  nclass <- nclass+1
}

## mult-mix-lwt
cmdline <- paste(bindir,"mult-mix-lwt -i tmp.high-rate-seqfile",
                 " -l tmp.lwt",
                 " -c ",nclass,
                 " -f tmp.frs",
                 ifelse(plusF," -d ",""),
                 " -p ",
                 ifelse(C>0,paste(" -C",C),""),
                 " > tmp.err",sep="")
## cat(cmdline,"\n")
system(cmdline)
lnl <-  scan("tmp.err",quiet=TRUE)

## Output
if(iq.freqs){
  fr <- matrix(scan("tmp.freq",quiet=TRUE),ncol=20,byrow=TRUE)
  CreateIQFreqfile(fr,"esmodel.nex")
}
did.rename <- file.rename("tmp.freq","estimated-frequencies")

## temporary files created:
tmp.files <- c("tmp.ctl","tmp.out","rate_est.dat","tmp.sitefile",
               "tmp.high-rate-seqfile","tmp.Sigma","tmp.frs","tmp.err",
               "tmp.wt","tmp.lwt")
did.remove <- file.remove(tmp.files)
