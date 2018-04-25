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
# import files
########
derepFs <- readRDS(file.path(temprds_path,"derepFs.RDS"))
dadaFs <- readRDS(file.path(temprds_path,"dadaFs.RDS"))
derepRs <- readRDS(file.path(temprds_path,"derepRs.RDS"))
dadaRs <- readRDS(file.path(temprds_path,"dadaRs.RDS"))

plot.err.f<-plotErrors(dadaFs,nominalQ=TRUE)
plot.err.r<-plotErrors(dadaRs,nominalQ=TRUE)
plot.err.model<-grid.arrange(plot.err.f,plot.err.r,nrow=1,ncol=2)
ggsave(file.path(tempfig_path,"dada.err.FULL.pdf"), plot.err.model,width=20,height=10,limitsize=FALSE)

#----------------------------------------------------------
########
# merge the inferred forward and reverse sequences, while removing paired sequnces that do not perfectly overlap as a final control against residual error
########
mergers <- mergePairs(
	dadaFs, 
	derepFs, 
	dadaRs, 
	derepRs,
	minOverlap=20,
	verbose=TRUE
)
head(mergers[[1]])
saveRDS(mergers,file.path(temprds_path,"mergers.RDS"))

#----------------------------------------------------------
########
# Construct the sequence table and remove chimeras
########
names(mergers)
# any effect of changing min overlap?
table(mergers[[1]]$nmatch)
seqtab.all <- makeSequenceTable(mergers)

########
# inspect the sequence table
########
dim(seqtab.all)
table(nchar(colnames(seqtab.all)))
summary(nchar(colnames(seqtab.all)))
pdf(file.path(tempfig_path,"seqtab.all.hist.pdf"))
hist(nchar(colnames(seqtab.all)),breaks=100)
abline(v=summary(nchar(colnames(seqtab.all)))[2],col="blue")
abline(v=summary(nchar(colnames(seqtab.all)))[3],col="red")
abline(v=summary(nchar(colnames(seqtab.all)))[5],col="blue")
dev.off()
# is there anything funky going on? Any surprises?
saveRDS(seqtab.all,file.path(temprds_path,"seqtab.all.RDS"))

########
# remove chimeras
########
seqtab <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)
# inspect sequence table with chimera's removed
dim(seqtab.all)
dim(seqtab)
table(nchar(colnames(seqtab)))
summary(nchar(colnames(seqtab)))
pdf(file.path(tempfig_path,"seqtab.hist.pdf"))
hist(nchar(colnames(seqtab)),breaks=100)
abline(v=summary(nchar(colnames(seqtab)))[2],col="blue")
abline(v=summary(nchar(colnames(seqtab)))[3],col="red")
abline(v=summary(nchar(colnames(seqtab)))[5],col="blue")
dev.off()
saveRDS(seqtab,file.path(temprds_path,"seqtab.RDS"))

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