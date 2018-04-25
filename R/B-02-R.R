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
tempfig_path<-file.path(figs_path,"Build_Phyloseq2")
temprds_path<-file.path(rds_path,"Build_Phyloseq2")

# prep environment and get data
.cran_packages<-c("ggplot2", "gridExtra")
.bioc_packages<-c("dada2", "phyloseq")
sapply(c(.cran_packages,.bioc_packages),require,character.only=TRUE)

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# import derep
########
derepRs <- readRDS(file.path(temprds_path,"derepRs.RDS"))

#----------------------------------------------------------
########
# use dada2 to model substitution errors and distinguish sequencing errors from real biological variation. 
########
errR<-learnErrors(derepRs, multithread=TRUE)
plot.err.r<-plotErrors(errR,nominalQ=TRUE)
ggsave(file.path(tempfig_path,"dada.err.r.pdf"),plot.err.r,width=10,height=10,limitsize=FALSE)

dadaRs<-dada(derepRs, err=errR, multithread=TRUE)
dada2:::checkConvergence(dadaRs[[1]])
saveRDS(dadaRs,file.path(temprds_path,"dadaRs.RDS"))
for (i in seq(1,length(dadaRs))) {
	print(names(dadaRs)[[i]])
	print(dadaRs[[i]])
}

#----------------------------------------------------------
########
# clear temporary paths
########
rm(tempfig_path)
rm(temprds_path)
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########