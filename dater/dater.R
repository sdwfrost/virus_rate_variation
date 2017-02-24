#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=T)
if(length(args) < 5)
{
  stop('Usage: rtt.R <stub> <ncpu> <nreps> <numStart> <searchRoot>')
}

suppressPackageStartupMessages(require(ape))
suppressPackageStartupMessages(require(magrittr))
suppressPackageStartupMessages(require(treedater))
suppressPackageStartupMessages(require(pander))
suppressPackageStartupMessages(require(lubridate))
panderOptions('digits', 4)
panderOptions('round', 4)
panderOptions('keep.trailing.zeros', TRUE)

stub <- args[1]
ncpu <- as.integer(args[2])
nreps <- as.integer(args[3])
numStart <- as.integer(args[4])
searchRoot <- as.integer(args[5])

trfn <- paste("../ml_trees/",stub,".phy_phyml_tree",sep="")
sfn <- paste("../alignments/",stub,".fasta",sep="")

if(grepl("root",stub)){
  tr <- read.nexus(trfn)
  s <- read.dna(gsub("rooted_","",sfn),format="fasta",as.matrix = FALSE)
}else{
  tr <- read.tree(trfn)
  s <- read.dna(sfn,format="fasta",as.matrix = FALSE)
}

sl <- length(s[[1]])

tn <- tr$tip.label
td <- tn %>% strsplit(.,"_") %>% lapply(.,tail,1) %>% unlist %>% as.double
sts <- td
names(sts) <- tn

strict <- dater(tr, sts, s=sl,
                minblen=1./365,
                temporalConstraints=TRUE,
                strictClock=TRUE,
                numStart=numStart,
                searchRoot=searchRoot,
                ncpu=ncpu)
#strict.ci <- parboot.treedater(strict,nreps=nreps,overrideTempConstraint=F,ncpu=ncpu)

relaxed <- dater(tr, sts, s=sl,
                minblen=1./365,
                temporalConstraints=TRUE,
                strictClock=FALSE,
                numStart=numStart,
                searchRoot=searchRoot,
                ncpu=ncpu)

results <- data.frame(Dataset=c(stub,stub),Clock=c("Strict","Relaxed"),
                      TMRCA=as_date(date_decimal(c(strict$timeOfMRCA,relaxed$timeOfMRCA))),
                      Rate=c(strict$meanRate,relaxed$meanRate),
                      CV=c(0,relaxed$coef_of_variation))
pander(results)

write.table(results,paste(stub,".txt",sep=""),col.names=TRUE,row.names=FALSE,sep="\t",quote=F)
save.image(file=paste(stub,".RData",sep=""))
