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

# prep environment and get data
analysis_prep(c("phyloseq","ggplot2","DESeq2"))

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# Differential Abundance modeling - DESeq2
########
########
# prep phyloseq obj for analysis
########
#save phyloseq object
ps_saved<-ps

names(sample_data(ps))
str(sample_data(ps))
sample_data(ps)$Experiment

stat<-c(
	rep("Operation",3),
	rep("Start",3),
	rep("Start",3),
	rep("Attenuation",3),
	rep("Attenuation",3),
	rep("Operation",3),
	rep("Failing",3),
	rep("Failing",3),
	rep("Crashed",3),
	"Attenuation",
	"Attenuation",
	"Failing",
	"Failing",
	"Failing",
	"Failing",
	"Failing",
	"Crashed")
stat<-factor(stat)

sample_data(ps)$Status<-stat
sample_data(ps)$Status

########
# run DESeq2
########
dds<-phyloseq_to_deseq2(ps, ~ Experiment + Status)
dds<-DESeq(dds, test="Wald", fitType="parametric")
resultsNames(dds)
results(dds)
summary(results(dds, contrast=c("Status","Operation","Start")))
summary(results(dds, contrast=c("Status","Crashed","Operation")))

########
# create results objects
########
# pull results
res.startup<-results(dds, contrast=c("Status","Operation","Start"),alpha=0.1)
res.crash<-results(dds, contrast=c("Status","Crashed","Operation"),alpha=0.1)
summary(res.startup)
summary(res.crash)
# add taxonomy
res.startup<-cbind(as(res.startup,"data.frame"),	as(tax_table(ps),"matrix"))
res.crash<-cbind(as(res.crash,"data.frame"),as(tax_table(ps),"matrix"))
# pull only significant taxa
sig.startup<-res.startup[which(res.startup$padj < 0.05),]
sig.crash<-res.crash[which(res.crash$padj < 0.05),]
dim(res.startup)
dim(sig.startup)
dim(res.crash)
dim(sig.crash)
# save results
saveRDS(res.startup,file=file.path(rds_path,"res.startup.RDS"))
saveRDS(res.crash,file=file.path(rds_path,"res.crash.RDS"))
saveRDS(sig.startup,file=file.path(rds_path,"sig.startup.RDS"))
saveRDS(sig.crash,file=file.path(rds_path,"sig.crash.RDS"))

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# create figures
########
library(theseus)
comp1lab <- c("Decreased 1 to 4", "No Change 1 to 4", "Increased 1 to 4")
comp2lab <- c("Decreased 4 to 7", "No Change 4 to 7", "Increased 4 to 7")
st14<-sig.startup
st47<-sig.crash
psra <- transform_sample_counts(ps,function(x) {x / sum(x)})

#----------------------------------------------------------
########
# Cohort_Phylum
########

p.Cohort_Phylum<-
cohort_relabund(
	PS=ps, comp1=st14, comp2=st47, 
	comp1lab = 
		c("Decreased 1 to 4", "No Change 1 to 4", "Increased 1 to 4"),
	comp2lab = 
		c("Decreased 4 to 7", "No Change 4 to 7", "Increased 4 to 7")) +
	ylab("Relative Abundance (Phylum)") +
	geom_vline(xintercept=35-5-0.5, linetype="solid", alpha=1.0, size=0.5) +
	scale_y_continuous(limits=c(0,0.50)) +
	geom_rect(aes(xmin=0.5,xmax=3,ymin=0.49,ymax=0.5), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=12,xmax=14,ymin=0.49,ymax=0.5), inherit.aes = FALSE, fill="blue") +
	geom_rect(aes(xmin=26,xmax=29,ymin=0.49,ymax=0.5), inherit.aes = FALSE, fill="green") +
	geom_rect(aes(xmin=30,xmax=32,ymin=0.49,ymax=0.5), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=33,xmax=35,ymin=0.49,ymax=0.5), inherit.aes = FALSE, fill="blue") +
	theme_bw() +
	theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0)) +
		theme(legend.position = "bottom") + 
		theme(axis.title.x=element_blank()) +
	theme(axis.text.y=element_text(angle=90, vjust=0, hjust=0.5)) +
	theme(axis.text.y = element_text(size=rel(0.8))) +
	theme(axis.text.x = element_text(size=rel(0.76))) +
	theme(axis.title.y = element_text(size=rel(0.8))) #+

p.Cohort_Phylum

ggsave(
	file.path(figs_path,"Cohort_Phylum.eps"),
	p.Cohort_Phylum,
	width=190, height=160, units="mm")

#----------------------------------------------------------
########
# Cohort_PAO-GAO-Cyto.eps
########
targ.g<-c("Aquicella","Candidatus_Accumulibacter", "Candidatus_Competibacter", "Candidatus_Protochlamydia", "Dechloromonas", "Defluviicoccus")
targ.c<-c("Cytophagia")

pstargra <- subset_taxa(psra, ((Genus %in% targ.g) | (Class %in% targ.c)))
pstargra
tax_table(pstargra)[tax_table(pstargra)[,3]=="Cytophagia",6] <- "Cytophagia (Class)"

ggsave(file.path(figs_path,"Cohort_PAO-GAO-Cyto.eps"),
cohort_relabund(PSisRelAbund=TRUE, taxfill="Genus",
	PS=pstargra, xvar="SampleID",comp1=st14, comp2=st47, 
	comp1lab = 
		c("Decreased 1 to 4", "No Change 1 to 4", "Increased 1 to 4"),
	comp2lab = 
		c("Decreased 4 to 7", "No Change 4 to 7", "Increased 4 to 7")) +
	ylab("Relative Abundance (Genus)") +
	geom_vline(xintercept=35-5-0.5, linetype="solid", alpha=1.0, size=0.5) +
	scale_y_continuous(limits=c(0,0.35)) +
	geom_rect(aes(xmin=0.5,xmax=3,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=12,xmax=14,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="blue") +
	geom_rect(aes(xmin=26,xmax=29,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="green") +
	geom_rect(aes(xmin=30,xmax=32,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=33,xmax=35,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="blue") +
	theme_bw() +
	theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0)) +
	theme(legend.position = "bottom") + 
	theme(axis.title.x=element_blank()) +
	theme(axis.text.y=element_text(angle=90, vjust=0, hjust=0.5)) +
 theme(axis.text.y = element_text(size=rel(0.8))) +
 theme(axis.text.x = element_text(size=rel(0.80))) +
 theme(axis.title.y = element_text(size=rel(0.8))) +
	theme(panel.spacing=unit(0.2,"lines"))
			 ,
			 width=190, height=160, units="mm")


#----------------------------------------------------------
########
# Cohort_PAO-GAO-Cyto_O_slim.eps
########
targ.g<-c("Aquicella","Candidatus_Accumulibacter", "Candidatus_Competibacter", "Dechloromonas")
targ.c<-c("Cytophagia")

pstargra <- subset_taxa(psra, ((Genus %in% targ.g) | (Class %in% targ.c)))
pstargra
tax_table(pstargra)[tax_table(pstargra)[,3]=="Cytophagia",6] <- "Cytophagia (Class)"

ggsave(file.path(figs_path,"Cohort_PAO-GAO-Cyto_O_slim.eps"),
cohort_relabund(PSisRelAbund=TRUE, taxfill="Genus",
	PS=prune_samples(sample_data(pstargra)$Experiment=="O",pstargra), 
	xvar="SampleID",comp1=st14, comp2=st47, 
	comp1lab = 
		c("Decreased 1 to 4", "No Change 1 to 4", "Increased 1 to 4"),
	comp2lab = 
		c("Decreased 4 to 7", "No Change 4 to 7", "Increased 4 to 7")) +
	ylab("Relative Abundance (Genus)") +
	scale_y_continuous(limits=c(0,0.35)) +
	geom_rect(aes(xmin=0.5,xmax=3,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=12,xmax=14,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="blue") +
	geom_rect(aes(xmin=26,xmax=29,ymin=0.34,ymax=0.35), inherit.aes = FALSE, fill="green") +
	theme_bw() +
	theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0)) +
	theme(legend.position = "bottom") + 
	theme(axis.title.x=element_blank()) +
	theme(axis.text.y=element_text(angle=90, vjust=0, hjust=0.5)) +
 theme(axis.text.y = element_text(size=rel(0.8))) +
 theme(axis.text.x = element_text(size=rel(0.80))) +
 theme(axis.title.y = element_text(size=rel(0.8))) +
	theme(panel.spacing=unit(0.2,"lines")) +
	guides(fill = guide_legend(nrow=2))
			 ,
			 width=190, height=160, units="mm")

#----------------------------------------------------------
########
# Cohort_PredBact.eps
########
targ.g<-c("Bdellovibrio")
targ.o<-c("Myxococcales")

pstargra <- subset_taxa(psra, ((Genus %in% targ.g) | (Order %in% targ.o)))
pstargra
tax_table(pstargra)[((tax_table(pstargra)[,4]=="Myxococcales") & (is.na(tax_table(pstargra)[,6]))),6] <- "Myxococcales (Order, Unannotated Genus)"

ggsave(file.path(figs_path,"Cohort_PredBact.eps"),
cohort_relabund(PSisRelAbund=TRUE, taxfill="Genus",
	PS=pstargra, xvar="SampleID",comp1=st14, comp2=st47, 
	comp1lab = 
		c("Decreased 1 to 4", "No Change 1 to 4", "Increased 1 to 4"),
	comp2lab = 
		c("Decreased 4 to 7", "No Change 4 to 7", "Increased 4 to 7")) +
	ylab("Relative Abundance (Genus)") +
	geom_vline(xintercept=35-5-0.5, linetype="solid", alpha=1.0, size=0.5) +
	scale_y_continuous(limits=c(0,0.02)) +
	geom_rect(aes(xmin=0.5,xmax=3,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=12,xmax=14,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="blue") +
	geom_rect(aes(xmin=26,xmax=29,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="green") +
	geom_rect(aes(xmin=30,xmax=32,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=33,xmax=35,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="blue") +
	theme_bw() +
	theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0)) +
	theme(legend.position = "bottom") + 
	theme(axis.title.x=element_blank()) +
	theme(axis.text.y=element_text(angle=90, vjust=0, hjust=0.5)) +
 theme(axis.text.y = element_text(size=rel(0.8))) +
 theme(axis.text.x = element_text(size=rel(0.80))) +
 theme(axis.title.y = element_text(size=rel(0.8))) +
	theme(panel.spacing=unit(0.2,"lines"))
			 ,
			 width=190, height=160, units="mm")

#----------------------------------------------------------
########
# Cohort_PredBact_O.eps
########

ggsave(file.path(figs_path,"Cohort_PredBact_O.eps"),
cohort_relabund(PSisRelAbund=TRUE, taxfill="Genus",
	PS=prune_samples(sample_data(pstargra)$Experiment=="O",pstargra), 
	xvar="SampleID",comp1=st14, comp2=st47, 
	comp1lab = 
		c("Decreased 1 to 4", "No Change 1 to 4", "Increased 1 to 4"),
	comp2lab = 
		c("Decreased 4 to 7", "No Change 4 to 7", "Increased 4 to 7")) +
	ylab("Relative Abundance (Genus)") +
	scale_y_continuous(limits=c(0,0.02)) +
	geom_rect(aes(xmin=0.5,xmax=3,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="red") +
	geom_rect(aes(xmin=12,xmax=14,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="blue") +
	geom_rect(aes(xmin=26,xmax=29,ymin=0.0195,ymax=0.02), inherit.aes = FALSE, fill="green") +
	theme_bw() +
	theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0)) +
	theme(legend.position = "bottom") + 
	theme(axis.title.x=element_blank()) +
	theme(axis.text.y=element_text(angle=90, vjust=0, hjust=0.5)) +
 theme(axis.text.y = element_text(size=rel(0.8))) +
 theme(axis.text.x = element_text(size=rel(0.80))) +
 theme(axis.title.y = element_text(size=rel(0.8))) +
	theme(panel.spacing=unit(0.2,"lines"))
			 ,
			 width=190, height=160, units="mm")


#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# Gather Data for tables
########
########
# Target species relative abundance
########
names(sample_data(ps))
cbind(sample_data(ps)$SampleGroup, sample_names(ps))
psmerge <- merge_samples(ps,"SampleGroup")
psmerge
sample_names(psmerge)

# merge samples
psmergera <- transform_sample_counts(psmerge,function(x) {x / sum(x)})

#subset taxa of interest
targ.g<-c("Aquicella","Candidatus_Accumulibacter", "Candidatus_Competibacter", "Candidatus_Protochlamydia", "Dechloromonas", "Defluviicoccus", "Bdellovibrio")
targ.c<-c("Cytophagia")
targ.o<-c("Myxococcales")

pstargra <- subset_taxa(psmergera, ((Genus %in% targ.g) | (Class %in% targ.c) | (Order %in% targ.o)))
pstargra
tax_table(pstargra)[tax_table(pstargra)[,3]=="Cytophagia",6] <- "Cytophagia (Class)"
tax_table(pstargra)[((tax_table(pstargra)[,4]=="Myxococcales") & (is.na(tax_table(pstargra)[,6]))),6] <- "Myxococcales (Order, Unannotated Genus)"
pstargra
# verify that everything is annotated
table(tax_table(pstargra)[,6],useNA="always")

dat<-cbind(
	tax_table(pstargra),
	t(otu_table(pstargra)))
dim(dat)
dat[1,]

write.csv(as.data.frame(dat), file=file.path(output_path,"TargSpec_RelAbund.txt"))

########
# cohort membership data
########
comp1<-st14
comp2<-st47
ls.comp1.up <- rownames(comp1[comp1$log2FoldChange>0,])
ls.comp1.down <- rownames(comp1[comp1$log2FoldChange<0,])
ls.comp2.up <- rownames(comp2[comp2$log2FoldChange>0,])
ls.comp2.down <- rownames(comp2[comp2$log2FoldChange<0,])
# make cohort data table
codat <- as.data.frame(tax_table(ps))
str(codat)
codat$comp1to4 <-NULL
codat$comp1to4 <- comp1lab[2]
codat$comp1to4[rownames(codat) %in% ls.comp1.up] <-comp1lab[3]
codat$comp1to4[rownames(codat) %in% ls.comp1.down] <-comp1lab[1]
table(codat$comp1to4)
codat$comp4to7 <-NULL
codat$comp4to7 <- comp2lab[2]
codat$comp4to7[rownames(codat) %in% ls.comp2.up] <-comp2lab[3]
codat$comp4to7[rownames(codat) %in% ls.comp2.down] <-comp2lab[1]
table(codat$comp4to7)
str(codat)
table(tax_table(ps)[,6])

dim(otu_table(psra))
str(otu_table(psra))
colnames(otu_table(psra)) == rownames(codat)
write.csv(as.data.frame(codat),file=file.path(output_path,"DiffAbund_CohortData.txt"))

