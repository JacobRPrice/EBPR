#----------------------------------------------------------
# Notes or References:
#
#----------------------------------------------------------
########
# Preliminary Items
########

# load utilities and functions
source("/PATH/TO/DIR/r-EBPR/R/utilities.R")
source("/PATH/TO/DIR/r-EBPR/R/functions.R")

# set temporary directories, if needed
tempfig_path<-file.path(figs_path,"Build_Phyloseq")
temprds_path<-file.path(rds_path,"Build_Phyloseq")

# prep environment and get data
.cran_packages<-c("ggplot2", "gridExtra","magrittr")
.bioc_packages<-c("dada2", "phyloseq", "DECIPHER", "phangorn")
sapply(c(.cran_packages,.bioc_packages),require,character.only=TRUE)

seqtab123<-readRDS(file.path(temprds_path,"seqtab123.RDS"))

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# Construct a phylogenetic tree
# make tree using updated protocol (v2)
########

# extract unique sequences from sequence table
seqs <- getSequences(seqtab123)
names(seqs) <- seqs 
alignment<-AlignSeqs(DNAStringSet(seqs),anchor=NA)
saveRDS(alignment,file.path(temprds_path,"alignment.RDS"))
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
saveRDS(phang.align,file.path(temprds_path,"phang.align.RDS"))
dm <- dist.ml(phang.align)
saveRDS(dm,file.path(temprds_path,"dm.RDS"))
treeNJ <- NJ(dm) 
saveRDS(treeNJ,file.path(temprds_path,"treeNJ.RDS"))
fit <- pml(treeNJ, data=phang.align) 
saveRDS(fit,file.path(temprds_path,"fit.RDS"))
fitGTR <- update(fit, k=4, inv=0.2)
saveRDS(fitGTR,file.path(temprds_path,"fitGTR.RDS"))
fitGTR <- optim.pml(
	fitGTR, 
	model="GTR", 
	optInv=TRUE, 
	optGamma=TRUE,
	rearrangement = "stochastic", 
	control = pml.control(trace = 0)
)
saveRDS(fitGTR,file.path(temprds_path,"fitGTR.RDS"))

#----------------------------------------------------------
########
# clear temporary paths
########
rm(temprds_path)
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########