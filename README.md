## ROH based kinship

The script is to compute ROH based kinship for all animals that have been haplotyped.  
The idea is based on the approach by [F. Gomez-Romano et al. (2016)](http://onlinelibrary.wiley.com/doi/10.1111/jbg.12213/epdf).  See the paper for the formulae.

### input data
       - Phased data with ID and Paternal and Maternal haplotypes for all SNPs either 'linkage' or 'beagle' file format
       - The marker map information

### Arguments
       - inpgenofile            :: The input phased 'linkage' or beagle format file 
       - mapinfo                :: MAP file (similar format and columns with PLINK map file. i.e. CHR, SNPNAME, GenPOS(cM), POS(BP)) 
       - inputformat            :: The input file format. i.e. is it 'linkage' or 'beagle' phased file format
       - ROHsizeNSNP            :: ROHsize in number of SNPs
       - Nummismatch            :: Number of heterozygotes allowed 
       - outputformat           :: output file format - 'asreml' or 'matrix' 
       - outROHcount            :: should the output be in ROH counts OR kiship estimate 
       - outname                :: output file name
       
##### outROHcount -- This argument is very important when you want to parallelize the computation in case the marker density and number of samples is large (>1000 animals > 50000 markers). In this case you can send job per chromosome and late combine by summing number of ROH across all the chromosomes. 
       You can compute this in R using the following code
       ################################################################################################################
       nanim=1000           ##### how many animals
       nsnps = 50000        ##### number of snps
       
       ##### empty ROH segment counts (useful later)
       ROHbasedcoans <- matrix(0,nrow=nanim,ncol=nanim)

       ##### list of files
       chromsdata <- c('','','')

       for (o in 1:length(chromsdata)){
         cat('\n... chromosome ...',o,' ....\n')
         chrROH <-  read.table(paste(chromsdata[o],sep=''),header=F,stringsAsFactors=F)
         ROHbasedcoans <- ROHbasedcoans+chrROH
       }
       diag(ROHbasedcoans) <- 1+(diag(ROHbasedcoans)/nsnps)
       ROHbasedcoans[lower.tri(ROHbasedcoans)] <- ROHbasedcoans[lower.tri(ROHbasedcoans)]/(4*nsnps)
       ####################################################################################################################
       
       
       
