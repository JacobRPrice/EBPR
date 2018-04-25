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

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# create sample data object
########
dat = read.csv(file.path(data_path,"MolecularSampleData.txt"),sep="\t",header=TRUE)
dat
names(dat)
dim(dat)

######## 
# perform calculations/conversions
######## 
dat$Day.cat<-as.factor(dat$Day)


######## 
# name rows so Sample ID's match
######## 
rownames(dat)<-dat$SampleName
dat$SampleID<-rownames(dat)

######## 
# inspect dat file before saving
######## 
dat
dim(dat)
names(dat)
attributes(dat)
str(dat,list.len=ncol(dat))

######## 
# Save dat object
######## 
saveRDS(dat,file.path(temprds_path,"samdf.RDS"))
write.table(dat,file.path(output_path,"samdf.txt"),sep="\t")

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