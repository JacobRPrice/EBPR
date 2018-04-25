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
.cran_packages<-c("ggplot2", "gridExtra")
.bioc_packages<-c("dada2", "phyloseq")
sapply(c(.cran_packages,.bioc_packages),require,character.only=TRUE)

seqtab<-readRDS(file.path(temprds_path,"seqtab123.RDS"))

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# assign taxonomy and create taxonomy table
########
########
# assign taxonomy - only goes to genus
taxtab.1to6<-assignTaxonomy(
	seqtab,
	refFasta=file.path(data_path,"assignTaxonomy/silva_nr_v123_train_set.fa.gz"),
	minBoot=80, # default=50
	verbose=TRUE,
	multithread=TRUE
)
head(unname(taxtab.1to6))
saveRDS(taxtab.1to6,
	file.path(temprds_path,"taxtab.1to6.RDS"))

########
# merge taxonomy and genus-species classification information
taxtab<-addSpecies(
	taxtab.1to6,
	refFasta=file.path(data_path,"assignTaxonomy/silva_species_assignment_v123.fa.gz"),
	verbose=TRUE,
	allowMultiple=TRUE # DO WE WANT TO ALLOW MULTIPLE? YES FOR NOW
)
head(unname(taxtab))
colnames(taxtab)
head(rownames(taxtab))
table(unname(taxtab)[,1],useNA="always")
table(unname(taxtab)[,2],useNA="always")
table(unname(taxtab)[,3],useNA="always")
table(unname(taxtab)[,4],useNA="always")
table(unname(taxtab)[,5],useNA="always")
table(unname(taxtab)[,6],useNA="always")
table(unname(taxtab)[,7],useNA="always")
saveRDS(taxtab,
	file.path(temprds_path,"taxtab.RDS"))

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