ROHbasedgrm <-function(inpgenofile,mapinfo,inputformat='beagle',ROHsizeNSNP=100,Nummismatch=0,outputformat='asreml',outROHcount=F,outname){
  ####### importing phased haplotypes
  if(inputformat=='linkage'){
    del <-  read.table(paste(inpgenofile,sep=''),header=F,stringsAsFactors=F,nrow=50)
    classes <- sapply(del,class)
    haptypes <- read.table(paste(inpgenofile,sep=''),header=F,stringsAsFactors=F,colClasses=classes)
    rm(classes,del); gc()
    animid <- as.vector(haptypes[,1])
    haptypes <- haptypes[,-1]
    cat('... phased haplotypes imported and processed ...\n')
    ### edits and infor gathering on data ###
    #animtest <- 
    animtest <- nrow(haptypes)
    pathaplo <- as.matrix(haptypes[1:animtest,seq(1,ncol(haptypes),2)]-1)
    mathaplo <- as.matrix(haptypes[1:animtest,seq(2,ncol(haptypes),2)]-1)
    nanim <- nrow(mathaplo)
    nsnps <- ncol(mathaplo)
    ids <- 1:nrow(mathaplo)
  } else if(inputformat=='beagle'){
    del <-  read.table(paste(inpgenofile,sep=''),skip=1,header=F,stringsAsFactors=F,nrow=50)
    classes <- sapply(del,class)
    haptypes <- read.table(paste(inpgenofile,sep=''),skip=1,header=F,stringsAsFactors=F,colClasses=classes)
    rm(classes,del); gc()
    animid <- read.table(paste(inpgenofile,sep=''),nrows=1,header=F,stringsAsFactors=F)[,-1:-2]
    animid <- as.vector(t(animid[,seq(1,ncol(animid),2)]))
    haptypes <- haptypes[,-1:-2]
    cat('... phased haplotypes imported and processed ...\n')
    haptypes <- apply(haptypes,2,function(x)as.numeric(as.factor(x)))
    dim(haptypes)
    ### edits and infor gathering on data ###
    animtest <- ncol(haptypes)
    pathaplo <- t(as.matrix(haptypes[,seq(1,ncol(haptypes),2)])-1)
    mathaplo <- t(as.matrix(haptypes[,seq(2,ncol(haptypes),2)])-1)
    dim(pathaplo)
    nanim <- nrow(mathaplo)
    nsnps <- ncol(mathaplo)
    ids <- 1:nrow(mathaplo)
  }
  
  ####### importing marker position for splitting into chromosomes ####
  mapfile <-  read.table(paste(mapinfo,sep=''),header=F,stringsAsFactors=F)
  colnames(mapfile) <- c('CHR','SNP','GenPOS','POS')
  mapfile <- mapfile[order(mapfile$CHR,mapfile$POS,decreasing = F),]
  mapfile$NR <- 1:nrow(mapfile)
  mapfilestrmin <- aggregate.data.frame(mapfile$NR,by=list(mapfile$CHR),min)
  mapfilestrmax <- aggregate.data.frame(mapfile$NR,by=list(mapfile$CHR),max)
  mapfilestr <- merge(mapfilestrmin,mapfilestrmax,by=1)
  colnames(mapfilestr) <- c('CHR','MIN','MAX')
  rm(mapfilestrmin,mapfilestrmax)
  chroms <- sort(unique(mapfile$CHR)); chroms <- chroms[which(chroms!=0)]
  cat('... marker cordinates/Positions imported and processed ...\n') 
  
  #########################################################
  #######        Functions for ROH detection    ###########
  #########################################################
  ##### ROH counts in the focal animal
  patmatROHs <- function(parhaplo,ROHsizeNSNP){
    ####### dirty tricks for fast row computation #######
    (a <- as.numeric(parhaplo)); (b=diff(a))
    (c=c(which(b!=0)+1,length(a)+1))
    (d=diff(c)); (d <- c(c[1]-1,d))
    (c <- c(1,c[-length(c)])); (e=(c+d))
    #####################################################
    Rsize <- data.frame(matrix(0,nrow=length(d),4))
    colnames(Rsize) <- c('start','end','ROHL','alleletype')
    Rsize$start <- c; Rsize$end <- e-1; Rsize$ROHL <- d
    Rsize$alleletype <-a[Rsize$start]
    Rsize <- Rsize[which(Rsize$ROHL>=ROHsizeNSNP),]
    return(Rsize)
  }
  
  ##### ROH counts in animal j 
  coansROH <- function(parhaploj,ROHsize,nHET){
    ROHshare <- matrix(0,nrow=nrow(ROHsize))
    for(s in 1:nrow(ROHsize)){
      rohsizehet <- ROHsize$ROHL[s]-nHET
      animjs <- parhaploj[ROHsize$start[s]:ROHsize$end[s]]
      nseqshare <- length(which(animjs==ROHsize$alleletype[s]))
      if(rohsizehet<=nseqshare){ROHshare[s] <- nseqshare}
    }
    return(sum(ROHshare))
  }
  #########################################################################################   
  ROHbasedcoans <- matrix(0,nrow=nanim,ncol=nanim)
  for (o in chroms){
    cat('\n... chromosome ...',o,' ....\n')
    strendchr <- as.numeric(mapfilestr[which(mapfilestr$CHR==o),c('MIN','MAX')])
    mathaplochr <- mathaplo[,c(strendchr[1]:strendchr[2])]
    pathaplochr <- pathaplo[,c(strendchr[1]:strendchr[2])]
    
    #### ROH GRM matrix
    ROHbasedcoanschr <- matrix(0,nrow=nanim,ncol=nanim)
    #### The number of combinations
    ncomb <- (nanim^2 + nanim)/2
    #### Iteration check printout
    iterchecks.ncomb <- round(ncomb/10,digits=0)
    
    for(i in ids){
      ##### determine ROHs for animal i (both paternal and maternal haplotype)  ####
      Rsize.pat <- patmatROHs(parhaplo=pathaplochr[i,],ROHsizeNSNP=ROHsizeNSNP)
      Rsize.mat <- patmatROHs(parhaplo=mathaplochr[i,],ROHsizeNSNP=ROHsizeNSNP)
      for(j in ids){
        ### dealing with diagonals (This is equivalent to ROHs from PLINK) ###
        ### Note that, some toher plink or ROH determinations thresholds are not used here ####
        if(i==j){
          if(nrow(Rsize.pat)==0 | nrow(Rsize.mat)==0){ ### if no ROH 
            coansROHij <- 0
          } else if(nrow(Rsize.pat)>0 & nrow(Rsize.mat)>0) { ### if there are ROHs
            ##### compare paternal and maternal haplotypes if they are IBD (same) ####
            ROHcountij <- coansROH(parhaploj=pathaplochr[j,],ROHsize=Rsize.mat,nHET=Nummismatch)
            ##### 
            coansROHij <- ROHcountij
          }
          ### sum up all ROHs for animal i with animal i
          ROHbasedcoanschr[j,i] <- coansROHij
          
          ############# dealing with off diagonals or ROH based kiship  ####################
        } else if(j>i){ 
          ##### compare paternal haplotype of animal i with paternal and maternal haplotype of animal j if they are IBD (same) ####
          if(nrow(Rsize.pat)==0){### if no ROH
            ROHcountjpat.pati <- 0
            ROHcountjmat.pati <- 0
          } else { ### if there are ROHs
            ROHcountjpat.pati <- coansROH(parhaploj=pathaplochr[j,],ROHsize=Rsize.pat,nHET=Nummismatch)
            ROHcountjmat.pati <- coansROH(parhaploj=mathaplochr[j,],ROHsize=Rsize.pat,nHET=Nummismatch)   
          }
          ##### compare maternal haplotype of animal i with with paternal and maternal haplotype of animal j if they are IBD (same) ####
          if(nrow(Rsize.mat)==0){ ### if no ROH
            ROHcountjpat.mati <- 0
            ROHcountjmat.mati <- 0
          } else { ### if there are ROHs
            ROHcountjpat.mati <- coansROH(parhaploj=pathaplochr[j,],ROHsize=Rsize.mat,nHET=Nummismatch)
            ROHcountjmat.mati <- coansROH(parhaploj=mathaplochr[j,],ROHsize=Rsize.mat,nHET=Nummismatch) 
          }
          ### sum up all ROHs from the 4-way comparison of animal i with animal j
          coansROHij <- sum(ROHcountjpat.pati,ROHcountjmat.pati,ROHcountjpat.mati,ROHcountjmat.mati)
          ROHbasedcoanschr[j,i] <- coansROHij
          if(ncomb %% iterchecks.ncomb==0){cat(ncomb,' combinations left ...\n')}
          ncomb <- ncomb-1
        }
      }
    }
    ROHbasedcoans <- ROHbasedcoans+ROHbasedcoanschr
  }
  if(outROHcount==T){
    write.table(ROHbasedcoans,paste(outname,'_ROHcount.txt',sep=''),quote=F,row.names=F,col.names=F)
    return(ROHbasedcoans)
  } else {
    diag(ROHbasedcoans) <- 1+(diag(ROHbasedcoans)/nsnps)
    ROHbasedcoans[lower.tri(ROHbasedcoans)] <- ROHbasedcoans[lower.tri(ROHbasedcoans)]/(4*nsnps)
    if(outputformat=='matrix'){
      write.table(ROHbasedcoans,paste(outname,'_ROH.grm',sep=''),quote=F,row.names=F,col.names=F)
    } else if(outputformat=='asreml'){
      rowcolwise <- ROHbasedcoans[lower.tri(ROHbasedcoans,diag=T)]
      ROHbasedcoansFF <- as.data.frame(which(row(ROHbasedcoans)>=col(ROHbasedcoans),arr.ind=TRUE))
      ROHbasedcoansFF$ROHG <- rowcolwise
      ROHbasedcoansFF <- ROHbasedcoansFF[,c(2,1,3)]
      ROHbasedcoansFF <- ROHbasedcoansFF[order(ROHbasedcoansFF[,2],ROHbasedcoansFF[,1]),]
      ROHbasedcoansFF <- ROHbasedcoansFF[,c(2,1,3)]
      write.table(ROHbasedcoansFF,paste(outname,'_ROHasreml.grm',sep=''),quote=F,row.names=F,col.names=F)
    }
    return(ROHbasedcoans)
  }
}
