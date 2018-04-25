#----------------------------------------------------------
#----------------------------------------------------------
# Notes or References:
# 
#----------------------------------------------------------
########
# FUNCTION: read raw data
########
raw_data<-function() {
	readRDS(file.path(data_path,"ps_orig.RDS"))
}

#----------------------------------------------------------
########
# FUNCTION: preprocess data (phyloseq object)
########
preprocessed_data<-function() {
	ps_orig<-raw_data()
	ps<-prune_taxa(taxa_sums(ps_orig)>0,ps_orig)
	ps<-subset_taxa(ps,!is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))
	ps<-subset_taxa(ps,!Phylum %in% c("Cyanobacteria"))
	ps<-filter_taxa(ps, function(x) sum(x > 3) >= 10, TRUE)
	list(ps_orig = ps_orig, ps = ps)
}

#----------------------------------------------------------
########
# FUNCTION: prepare the R environment, load packages, and import phyloseq objects
########
analysis_prep<-function(pkgs) {
	all_obj<-ls(envir=.GlobalEnv)
	all_obj<-setdiff(
			all_obj,
			c(
				"proj_path", "data_path", "figs_path", "R_path", "output_path", "rds_path",
				"raw_data", "preprocessed_data", "analysis_prep", "veganotu", "vegansd"
			))
	rm(list=all_obj,envir=.GlobalEnv)
	sapply(pkgs, require, character = TRUE)
	attach(preprocessed_data())
}

#----------------------------------------------------------
########
# vegan related functions
########

veganotu <- function(physeq) {
	require("vegan")
	OTU <- otu_table(physeq)
	if (taxa_are_rows(OTU)) {
		OTU <- t(OTU)
	}
	return(as(OTU, "matrix"))
}

vegansd <- function(physeq) {
	require("vegan")
	sd <- sample_data(physeq)
	return(as(sd,"data.frame"))
}

#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
#----------------------------------------------------------
########
# 
########