!/bin/bash

date
Rscript ./R/A-01.R 
Rscript ./R/B-01.R 
Rscript ./R/C-01.R 


date
Rscript ./R/A-02-F.R 
Rscript ./R/A-02-R.R 


date
Rscript ./R/B-02-F.R 
Rscript ./R/B-02-R.R 


date
Rscript ./R/C-02-F.R 
Rscript ./R/C-02-R.R 


date
Rscript ./R/A-03.R 
Rscript ./R/B-03.R 
Rscript ./R/C-03.R 

# merge seqtab together
date
Rscript ./R/04.R 

# build phy tree and tax table
date
Rscript ./R/05-PhyTree.R 
Rscript ./R/05-Tax.R 


# prep sample dataframe
date
Rscript ./R/05-SampDF.R 

# build phyloseq object
date
Rscript ./R/06.R
