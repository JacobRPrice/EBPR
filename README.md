# Evidence of predation and parasitism affecting EBPR performance through microbial community instability.

## Introduction
This repository contains all of the scripts needed to replicate the results and output presented in:   

Yaping, N., Price, J.R., Wang, Y., Cheng, M., Keshani Langroodi, S., Woloszynek, S., Rosen, G., Yuan, L., Sales, C.M. **2018 (IN REVIEW)**. Evidence of predation and parasitism affecting EBPR performance through microbial community instability.

Please read this markdown file in it's entirety before proceeding. 

All scripts have been written for and tested on Mac machines; some modification may be required for use with PC or Linux/Unix systems.

___  

## Replicate bioinformatics pipeline  
Follow the instructions in this section to carry out all of the bioinformatics processes up to and including the creation of the phyloseq object.
1) Download or clone this repository to your local machine
2) Raw reads - verify they are located in the correct directory    
* Raw sequencing files are currently hosted on github for reviewers/users convenience. Github may request the removal of this data.   
* Raw sequence reads have also been uploaded to the NCBI Sequence Read Archive (SRA ID: SRP131717, Bioproject ID: PRJNA431315)   
* Raw sequencing files should be placed within the "./data/rawreads/" directory
3) SILVA files   
* Users must obtain a SILVA liscense and copies of the SILVA training and species assignment datasets.   
* These datasets should be placed in a directory called "assignTaxonomy" within the data/ directory.    
* The specific name/version of training set used for this paper: silva_nr_v123_train_set.fa.gz    
* The specific name/version of the species assignment set used for this paper: silva_species_assignment_v123.fa.gz   
* Refer to dada2's taxonomic training page: http://benjjneb.github.io/dada2/training.html   
4) **Change directory paths before starting**   
Prior to carrying out any replication of results, users should change the file paths found in the scripts, i.e. change "PATH/TO/DIR/" so that R is able to locate data and scripts as they are sourced.   
5) Execute BioinformaticsPathway.bash 
* The file "BioinformaticsPathway.bash" is a bash script that will sequentially execute the required scripts in the "/R/" directory. This step is computationally heavy. It is suggested that users open up a second terminal window and excecute the "caffeinate" command to prevent the computer from going to sleep. 

___  

## Perform analysis peice-meal
A copy of the phyloseq object is included in the ./data/ directory. Once the phyloseq object has been successfully created (or you decide to use the one we've uploaded), use the provided scripts to replicate the results from this manuscript.   

The scripts for the results section generally follow this order:   
* Kinetics.R  
* AlphaDiversity.R    
* RelativeAbundance.R    
* Ordination_exploratory.R  
* DiffAbund.R  
