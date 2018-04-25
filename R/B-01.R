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

# set paths for raw data
# raw read directory
rawread_path<-file.path(data_path,"rawreads2")
# filtered read directory
filt_path<-file.path(data_path,"filtreads2")
if(!file_test("-d",filt_path)) dir.create(filt_path)

# set temporary directories, if needed
tempfig_path<-file.path(figs_path,"Build_Phyloseq2")
if(!file_test("-d",tempfig_path)) dir.create(tempfig_path)
temprds_path<-file.path(rds_path,"Build_Phyloseq2")
if(!file_test("-d",temprds_path)) dir.create(temprds_path)

# prep environment and get data
.cran_packages<-c("ggplot2", "gridExtra")
.bioc_packages<-c("dada2", "phyloseq")
sapply(c(.cran_packages,.bioc_packages),require,character.only=TRUE)

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# import sequences
########
fns <- sort(list.files(rawread_path, full.names = TRUE))
fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]
head(fns)
head(fnFs)
head(fnRs)

#----------------------------------------------------------
########
# trim and filter 
########

########
# visualize quality of sequences
########
# plots of forward reads
plot.qual.f<-NULL
plot.qual.f<-list()
for(i in 1:21) { 
	plot.qual.f[[i]]<-plotQualityProfile(fnFs[i]) + 
	ggtitle("Fwd") + 
	geom_hline(yintercept=30, size=0.5) + 
	geom_hline(yintercept=25, size=0.5) + 
	geom_hline(yintercept=20, size=0.5) +
	geom_vline(xintercept=230, size=0.5, color="red") + 
	geom_vline(xintercept=70, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=80, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=90, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=100, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=110, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=120, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=130, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=140, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=150, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=160, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=170, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=180, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=190, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=200, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=210, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=220, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=230, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=240, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=250, size=0.5, color="gray",linetype="dashed") +
	xlim(0,250) + ylim(0,40)
}
plot.read.qual.f<-grid.arrange(
	plot.qual.f[[1]],plot.qual.f[[2]],plot.qual.f[[3]],plot.qual.f[[4]],
	plot.qual.f[[5]],plot.qual.f[[6]],plot.qual.f[[7]],plot.qual.f[[8]],
	plot.qual.f[[9]],
	plot.qual.f[[10]],plot.qual.f[[11]],plot.qual.f[[12]],plot.qual.f[[13]],
	plot.qual.f[[14]],plot.qual.f[[15]],plot.qual.f[[16]],plot.qual.f[[17]],
	plot.qual.f[[18]],plot.qual.f[[19]],
	plot.qual.f[[20]],plot.qual.f[[21]], nrow=3,ncol=7)
ggsave(file.path(tempfig_path,"read.qual.f.pdf"),plot.read.qual.f,width=2*7,height=2*3,limitsize=FALSE)

# plots of reverse reads
plot.qual.r<-NULL
plot.qual.r<-list()
for(i in 1:21) { 
	plot.qual.r[[i]]<-plotQualityProfile(fnRs[i]) + 
	ggtitle("Rev") + 
	geom_hline(yintercept=30, size=0.5) + 
	geom_hline(yintercept=25, size=0.5) + 
	geom_hline(yintercept=20, size=0.5) +
	geom_vline(xintercept=220, size=0.5, color="red") + 
	geom_vline(xintercept=70, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=80, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=90, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=100, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=110, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=120, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=130, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=140, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=150, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=160, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=170, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=180, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=190, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=200, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=210, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=220, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=230, size=0.5, color="gray",linetype="dashed") + 
	geom_vline(xintercept=240, size=0.5, color="gray",linetype="dashed") +
	geom_vline(xintercept=250, size=0.5, color="gray",linetype="dashed") +
	xlim(0,250) + ylim(0,40)
}
plot.read.qual.r<-grid.arrange(
	plot.qual.r[[1]],plot.qual.r[[2]],plot.qual.r[[3]],plot.qual.r[[4]],
	plot.qual.r[[5]],plot.qual.r[[6]],plot.qual.r[[7]],plot.qual.r[[8]],
	plot.qual.r[[9]],
	plot.qual.r[[10]],plot.qual.r[[11]],plot.qual.r[[12]],plot.qual.r[[13]],
	plot.qual.r[[14]],plot.qual.r[[15]],plot.qual.r[[16]],plot.qual.r[[17]],
	plot.qual.r[[18]],plot.qual.r[[19]],
	plot.qual.r[[20]],plot.qual.r[[21]], nrow=3,ncol=7)
ggsave(file.path(tempfig_path,"read.qual.r.pdf"),plot.read.qual.r,width=2*7,height=2*3,limitsize=FALSE)

########
# carry out trimming and filtering of reads
########
# make list of names with paths for filtered read files
filtFs<-file.path(filt_path,basename(fnFs))
filtRs<-file.path(filt_path,basename(fnRs))
cbind(fnFs,filtFs)
cbind(fnRs,filtRs)

# # carry out trimming and filtering of reads
for(i in seq_along(fnFs)) {
  fastqPairedFilter(
  	c(fnFs[[i]], fnRs[[i]]),
  	c(filtFs[[i]], filtRs[[i]]),
  	trimLeft=10, 
  	truncLen=c(230,220),
  	maxN=0, # none allowed for dada2
  	maxEE=2, 
  	truncQ=2,
  	compress=TRUE,
  	verbose=TRUE)
}

#----------------------------------------------------------
########
# dereplicate forward and reverse reads
########
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
sam.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
cbind(basename(fnFs),basename(fnRs))
cbind(basename(filtFs),basename(filtRs))
cbind(
	basename(filtFs),
	sam.names,
	c("O1.1", "O1.2", "O1.3", "O2.1", "O2.2_A", "O2.3", 
	"O3.1_A", "O3.2", "O3.3", "O4.1", "O4.2", "O4.3", 
	"O5.1", "O5.2_A", "O5.3_A", "O6.1_A", "O6.2_A", "O6.3_A",
	"O7.1_A", "O7.2", "O7.3")
)
sam.names<-c("O1.1", "O1.2", "O1.3", "O2.1", "O2.2_A", "O2.3", 
	"O3.1_A", "O3.2", "O3.3", "O4.1", "O4.2", "O4.3", 
	"O5.1", "O5.2_A", "O5.3_A", "O6.1_A", "O6.2_A", "O6.3_A",
	"O7.1_A", "O7.2", "O7.3")
names(derepFs) <- sam.names
names(derepRs) <- sam.names
derepFs[[1]]
derepRs[[1]]
saveRDS(derepFs,file.path(temprds_path,"derepFs.RDS"))
saveRDS(derepRs,file.path(temprds_path,"derepRs.RDS"))

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