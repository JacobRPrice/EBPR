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
analysis_prep(c("phyloseq","ggplot2","DESeq2","vegan"))

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# prep data
########
# save unaltered ps 
psnorm<-ps
ps_orig
ps

# prep dds object
dds<-phyloseq_to_deseq2(ps,~1)
dds
dds<-estimateSizeFactors(dds)
dds<-estimateDispersions(dds)
vst<-getVarianceStabilizedData(dds)
dim(vst)
dim(otu_table(ps))
vst[vst < 0.0] <- 0.0
otu_table(ps)<-otu_table(vst, taxa_are_rows=TRUE)

#----------------------------------------------------------
########
# initial clustering and ordinations 
# getting an idea of what the data looks like
########

dim(vst)
dis<-vegdist(t(vst), method="bray")
clus<-hclust(dis, method="single")
plot(clus)
cluc<-hclust(dis, method="complete")
plot(cluc)
clua<-hclust(dis, method="average")
plot(clua)
cor(dis,cophenetic(clus))
cor(dis,cophenetic(cluc))
cor(dis,cophenetic(clua)) # average looks best
plot(clua, hang=-1)
rect.hclust(clua,5) # 5 clusters joins seqences nicely
grp<-cutree(clua,5)
grp

#----------------------------------------------------------
########
# ordinations
########
sample_data(ps)$Cluster<-as.factor(grp)

########
# Principal Components Analysis (PCA)
########
out<-ordinate(ps,method="MDS",distance="bray")
evals<-out$values$Eigenvalues
ptmat<-vegan::scores(out.dpcoa, display="sites")
ptdf<-data.frame(labels=rownames(ptmat),ptmat)
ptmap<-aes(xend=Axis1,yend=Axis2,x=0,y=0,shape=NULL,col=NULL,label=labels)
labelmap<-aes(x=4*Axis1,y=4*Axis2,shape=NULL,col=NULL,label=labels)
plot.pca<-
	  plot_ordination(ps, out, color="Cluster", shape="Experiment") +  
  coord_fixed(sqrt(evals[2]/evals[1])) +
  theme(legend.position="bottom") + 
  geom_text_repel(aes(Axis.1, Axis.2, label = SampleID))
ggsave(
	file.path(figs_path,"PCoA.bray.eps"), 
	plot.pca, 
	width=190, height=190/2, units="mm")

########
# DPCoA 
########
out.dpcoa<-ordinate(ps,method="DPCoA")
evals<-out.dpcoa$eig
ptmat<-vegan::scores(out.dpcoa, display="sites")
ptdf<-data.frame(labels=rownames(ptmat),ptmat)
ptmap<-aes(xend=Axis1,yend=Axis2,x=0,y=0,shape=NULL,col=NULL,label=labels)
labelmap<-aes(x=4*Axis1,y=4*Axis2,shape=NULL,col=NULL,label=labels)
library(ggrepel)
p.dpcoa.overlay<-
  plot_ordination(ps,
                  out.dpcoa,
                  type="taxa",
                  color="Order") +

	geom_text_repel(labelmap, data=ptdf) +
	geom_point(mapping=labelmap, data=ptdf) +
  theme(legend.position="bottom") + 
  theme(legend.title=element_blank())
ggsave(
	file.path(figs_path,"DPCoA_overlay.eps"),p.dpcoa.overlay,
	width=190, height=220, units="mm")

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########