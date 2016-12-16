########################################################################
################## Example file generated with QMSim   #################
########################################################################
system('QMSim.exe TestROH.prm')


#### Testing a file phased simulated data  ####
inpgenofile='example.linkgae'
mapinfo='example.map'
ROHsizeNSNP=100
Nummismatch=0
inputformat='beagle'
outputformat='asreml'
outROHcount=F
outname='example'
source('ROHbased_grm.R')
exampleROHgrm <- ROHbasedgrm(inpgenofile,mapinfo,inputformat,ROHsizeNSNP,
                          Nummismatch,outputformat,outROHcount,outname)



#### Testing a file phased and imputed with BEAGLE  ####
inpgenofile = 'Rawchr29.bgl.gz'
mapinfo='Rawchr29.bim'
ROHsizeNSNP=50
Nummismatch=1
inputformat='beagle'
outputformat='asreml'
outROHcount=F
outname='Raw'
source('ROHbased_grm.R')
RawROHgrm <- ROHbasedgrm(inpgenofile,mapinfo,inputformat,ROHsizeNSNP,
                      Nummismatch,outputformat,outROHcount,outname)

