system('QMSim.exe TestROH.prm')

#map <- read.table('TestROH/lm_mrk_001.txt',header=T,stringsAsFactors = F)
#map$pos <- format(map$Position*1000000,scientific = F)
#write.table(map[,c(2,1,3,4)],'ROH.map',quote = F,row.names = F,col.names = F)

inpgenofile = 'TestROH/ROHpop_mrk_001.txt'
mapinfo='ROH.map'
ROHsizeNSNP=30
source('ROHbased_grm.R')
ROHgrm <- ROHbasedgrm(inpgenofile,mapinfo,ROHsizeNSNP,Nummismatch=0,outputformat='asreml',outname='test')


Z <- mathaplo + pathaplo
afreq <- colMeans(Z)/2
K <- sum(2*afreq*(1-afreq))
Z <- scale(Z,center=T,scale=F)
G <- tcrossprod(data.matrix(Z))/K
