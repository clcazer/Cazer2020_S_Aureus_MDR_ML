# ---
# title: "Master Do"
# author: "Casey Cazer"
# Last updated: May 4, 2020
# ---

#this script runs several analysis for creating a multidrug resistance antibiogram using association mining
#specifically, this analysis accompanies the manuscript: 
##"Analysis of multidrug resistance with a machine learning generated antibiogram for Staphylococcus aureus at one New York hospital between 2007 and 2018"

#install packages required for analysis
##Rtools and Java may be required
install.packages(c("Amelia", "visdat", "naniar", "arules", "plyr", "dplyr", "tidyr", "purrr", "tibble", "stringr",
 "ggplot2", "xlsx", "scales", "ggpubr", "igraph", "splitstackshape",
 "DataCombine", "RColorBrewer", "Cairo", "mgsub", "knitr", "png", "gridExtra"))

#library path
normalizePath(.libPaths(), winslash="/")

#check the packages installed
installed.packages(.libPaths()[1])[,"Package"]

#libraries and functions required
library(plyr)
library(tidyr)
library(purrr)
library(reshape2)
library(arules)
library(dplyr)
library(tibble)
library(stringr)
library(ggplot2)
library(scales)
library(ggpubr)
library(igraph)
library(splitstackshape)
library(mgsub)
library(RColorBrewer)
library(Cairo)
library(png)
library(grid)
library(gridExtra)
library(xlsx)
library(knitr)


#load custom functions
knit('scripts/Mining Functions.Rmd')
knit('scripts/PvalueSets.Rmd')
source("scripts/AMR Data Cleaning Functions.R",local=FALSE, echo=TRUE, spaced=TRUE)
save.image("Rdata/mining functions.RData")
#load("Rdata/mining functions.RData")


#prepare data for analysis
source("scripts/data prep.R",local=FALSE, echo=TRUE, spaced=TRUE)

#create antibiogram and analyze prevalence of multidrug resistance
source("scripts/antibiogram and mdr profiles.R",local=FALSE, echo=TRUE, spaced=TRUE)
    ##will warn that mean.default(., na.rm=TRUE) argument is not logical for antibiogram by Infection Type.
    ##this is ok, it is trying to calculate prevalence of Infection.Type, this column is later dropped

#mine association sets    
source("scripts/mining sets.R",local=FALSE, echo=TRUE, spaced=TRUE)

#calculate P-values for each association set
source("scripts/PValues.R",local=FALSE, echo=TRUE, spaced=TRUE)

#filter sets based on P-values
source("scripts/filtering and summarizing sets.R",local=FALSE, echo=TRUE, spaced=TRUE)
    ##warning of QM_for_all_itemsets about NA's is normal
  
#calculate bootstrap confidence intervals on set quality measure
source("scripts/bootstrap.R",local=FALSE, echo=TRUE, spaced=TRUE)

#visualize the association sets by antimicrobial drugs
source("scripts/AM circle plots.R",local=FALSE, echo=TRUE, spaced=TRUE)

#visualize the association sets by antimicrobial class
source("scripts/class circle plots.R",local=FALSE, echo=TRUE, spaced=TRUE)
#will warn NAs introduced by coercion--ok, this occurs for node classes but these are the same as node names

#visualize the association sets of SST by MRSA v MSSA
source("scripts/SST MRSA v MSSA.R",local=FALSE, echo=TRUE, spaced=TRUE)

#descriptive statistics for manuscript
source("scripts/descriptive analysis.R",local=FALSE, echo=TRUE, spaced=TRUE)

#for methods appendix, check example tables
source("scripts/methods example.R",local=FALSE, echo=TRUE, spaced=TRUE)

#save environment
save.image("RData/SA year and MRSA MSSA and Type Analysis.RData")
