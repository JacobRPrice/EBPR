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
temprds_path<-file.path(rds_path,"Build_Phyloseq")

# prep environment and get data
.bioc_packages<-c("dada2", "phyloseq", "DECIPHER", "phangorn")
sapply(c(.bioc_packages),require,character.only=TRUE)

taxtab<-readRDS(file.path(temprds_path,"taxtab.RDS"))
samdf<-readRDS(file.path(temprds_path,"samdf.RDS"))
seqtab123<-readRDS(file.path(temprds_path,"seqtab123.RDS"))
fitGTR<-readRDS(file.path(temprds_path,"fitGTR.RDS"))

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# combine all produced objects into one PhyloSeq object
########

ps_orig <- phyloseq(
	tax_table(taxtab),
	sample_data(samdf),
	otu_table(seqtab123, taxa_are_rows = FALSE),
	phy_tree(fitGTR$tree)
)
ps_orig
sample_names(ps_orig)
# save phyloseq object in temprds_path as a backup
saveRDS(ps_orig,file.path(temprds_path,"ps_orig_BAK.RDS"))
# save phyloseq object in data_path for use going forward
saveRDS(ps_orig,file.path(data_path,"ps_orig.RDS"))

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