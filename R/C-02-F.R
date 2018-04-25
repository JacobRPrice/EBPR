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
tempfig_path<-file.path(figs_path,"Build_Phyloseq3")
temprds_path<-file.path(rds_path,"Build_Phyloseq3")

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
derepFs <- readRDS(file.path(temprds_path,"derepFs.RDS"))

#----------------------------------------------------------
########
# use dada2 to model substitution errors and distinguish sequencing errors from real biological variation. 
########
errF<-learnErrors(derepFs, multithread=TRUE)
plot.err.f<-plotErrors(errF,nominalQ=TRUE)
ggsave(file.path(tempfig_path,"dada.err.f.pdf"),plot.err.f,width=10,height=10,limitsize=FALSE)

dadaFs<-dada(derepFs, err=errF, multithread=TRUE)
dada2:::checkConvergence(dadaFs[[1]])
saveRDS(dadaFs,file.path(temprds_path,"dadaFs.RDS"))
for (i in seq(1,length(dadaFs))) {
	print(names(dadaFs)[[i]])
	print(dadaFs[[i]])
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