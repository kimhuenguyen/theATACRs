###################################################
#######This is a script by theATACRs to analyse whether chromosome accessibility of a promotor region correlates with RNA expression



###########install GEOquery
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GEOquery")

