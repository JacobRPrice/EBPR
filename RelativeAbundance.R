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
# pkgs.to.load<-c("phyloseq","ggplot2","gridExtra","DESeq2","ggrepel")
# analysis_prep(pkgs.to.load)
analysis_prep(c("phyloseq","ggplot2","gridExtra"))
#theme_set(theme_bw())

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# Phylum level
########
psra<-transform_sample_counts(ps,function(x) {x / sum(x)})
ps_orig
psra
psra.mergeP<-tax_glom(psra,taxrank="Phylum",NArm=FALSE)
psra.mergeP
names(sample_data(psra.mergeP))

plot.p<-
	plot_bar(
			prune_taxa(
				(taxa_sums(psra.mergeP)/nsamples(psra.mergeP))>0.01,
				psra.mergeP),
			x="SampleID",
			fill="Phylum"#,
#			facet_grid=Experiment~.
			) +
	theme_bw() +
	theme(axis.title.x=element_blank()) +
  geom_vline(xintercept=35-5-0.5, linetype="solid", alpha=1.0, size=1.5) +
  theme(legend.position="bottom") +
	ylab("Relative Abundance (Phylum)") +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	#xlab("Condition") +
	theme(axis.text.x=element_text(angle=90,vjust=0.5)) #+
		# theme_bw()

#----------------------------------------------------------
########
# Order level
########
ps_orig
psra
psra.mergeO<-tax_glom(psra,taxrank="Order",NArm=FALSE)
psra.mergeO

plot.o<-
	plot_bar(
			prune_taxa(
				(taxa_sums(psra.mergeO)/nsamples(psra.mergeO))>0.01,
				psra.mergeO),
			x="SampleID",
			fill="Order"#,
#			facet_grid=Experiment~.
			) +
	theme_bw() +
	theme(axis.title.x=element_blank()) +
  scale_y_continuous(limits = c(0, 1.0)) +
	geom_vline(xintercept=35-5-0.5, linetype="solid", alpha=1.0, size=1.5) +
  theme(legend.position="bottom") +
	ylab("Relative Abundance (Order)") +
	theme(axis.text.y=element_text(angle=90,hjust=0.5)) +
	#xlab("Condition") +
	theme(axis.text.x=element_text(angle=90,vjust=0.5)) #+
		# theme_bw()

#----------------------------------------------------------
########
# plot phylum and order together
########

ggsave(
	file.path(figs_path, "RelAbund_Phylum&Order.eps"),
	cowplot::plot_grid(
		plot.p +
			#theme(legend.title=element_blank()) + 
			theme(legend.text=element_text(size=rel(0.8))) + 
			theme(axis.text.x=element_text(size=rel(0.8))) + 
			theme(axis.text.y=element_text(size=rel(0.8))) +
			theme(axis.title.y=element_text(size=rel(0.9))), 
		plot.o+
			#theme(legend.title=element_blank()) + 
			theme(legend.text=element_text(size=rel(0.8))) + 
			theme(axis.text.x=element_text(size=rel(0.8))) + 
			theme(axis.text.y=element_text(size=rel(0.8))) +
			theme(axis.title.y=element_text(size=rel(0.9))),
		align="v", nrow=2, rel_heights=c(9/20, 11/20)
	)
	,
	width=190, height=210, units="mm")


#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########