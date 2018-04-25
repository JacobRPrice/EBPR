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
analysis_prep(c("phyloseq","ggplot2","gridExtra"))

# set seed
set.seed(19570530)

#----------------------------------------------------------
########
# Alpha Diversity
########
AlphaDiversity<-
  plot_richness(ps_orig, measures=c("Chao1","Shannon"), 
                color="SeqRun", shape="Experiment", x="Day", nrow=2) + 
  #  theme_bw() +
  theme(legend.position = "right") +
	geom_smooth(method=lm) +
	theme(axis.text.x=element_text(angle=90, hjust=0.5)) +
	theme(axis.text.y=element_text(angle=90, hjust=0.5)) 
  

ggsave(
	file.path(figs_path,"AlphaDiversity.pdf"), 
	AlphaDiversity, 
	width=190, height=120, units="mm")

AlphaDiversity$data
AlphaDiversity$data$SampleID
AlphaDiversity$data$variable
AlphaDiversity$data$value

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########