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
if(!file_test("-d",tempfig_path)) dir.create(tempfig_path)
temprds_path<-file.path(rds_path,"Build_Phyloseq")
if(!file_test("-d",temprds_path)) dir.create(temprds_path)

# prep environment and get data
.cran_packages<-c("ggplot2", "gridExtra")
.bioc_packages<-c("dada2", "phyloseq")
sapply(c(.cran_packages,.bioc_packages),require,character.only=TRUE)

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# import files
########
seqtab1<-readRDS(file.path(file.path(rds_path,"Build_Phyloseq1"),"seqtab.RDS"))
seqtab2<-readRDS(file.path(file.path(rds_path,"Build_Phyloseq2"),"seqtab.RDS"))
seqtab3<-readRDS(file.path(file.path(rds_path,"Build_Phyloseq3"),"seqtab.RDS"))

#----------------------------------------------------------
########
# merge seqtabs
########
seqtab123<-mergeSequenceTables(seqtab1,seqtab2,seqtab3)
saveRDS(seqtab123,file.path(temprds_path,"seqtab123.RDS"))

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