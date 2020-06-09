# title: "Descriptive Analysis"
# author: "Casey Cazer"
# Last updated: May 4, 2020

sink("results/descriptive.rtf", append=FALSE)
#this script provides descriptive output for the manuscript

###in abstract/intro/methods
#number of isolates
nrow(SA.db)

#tet-clinda resistance in SST (also cited in results)
filter(all.sets_combined_boot, Category=="SA.SST", items=="Clindamycin,Tetracycline")
#were these isolates MRSA or MSSA? (cited in results)
filter(SA.SST_db, Clindamycin==TRUE, Tetracycline==TRUE) %>% select(Oxacillin)

#number of antimicrobials tested
length(AM_include)

#years covered
min(SA.db$Study.Year)
max(SA.db$Study.Year)

#number of antimicrobials with clincal breakpoints
length(AM_include_with_bp)


###in Results
#overall prevalence of methicillin-resistance (from antibiogram, SupTable 1)
filter(as.data.frame(SupTable1), Category=="Overall") %>% select("Oxacillin")

#overall prevalence of MDR (>= 3 class resistances)
filter(mdr.classtable, !(NumClassRes %in% c("0", "1", "2"))) %>% summarize(sum=sum(NumIsolates)/sum(mdr.classtable$NumIsolates))

#overall prevalence of pan-susceptible
filter(mdr.classtable, NumClassRes=="0") %>% summarize(sum=sum(NumIsolates)/sum(mdr.classtable$NumIsolates))

#five most common class resistance pattern, with cLift interval excluding 1
filter(all.sets_combined_boot, LiftCrosses1==FALSE) %>% group_by(ClassCodes) %>% dplyr::summarise(count=n(), avgcLift=mean(cLift)) %>% arrange(desc(count)) %>% head(5)

#BL*FQ patterns by category (e.g. pneumonia vs MRSA)
filter(all.sets_combined_boot, LiftCrosses1==FALSE, ClassCodes=="BL,FQ") %>% group_by(Category) %>% dplyr::summarise(count=n(), avgcLift=mean(cLift)) %>% arrange(desc(count))

#top 5 class resistance patterns with greatest average cLift
filter(all.sets_combined_boot, LiftCrosses1==FALSE) %>% group_by(ClassCodes) %>% dplyr::summarise(count=n(), avgcLift=mean(cLift)) %>% arrange(desc(avgcLift)) %>% head(5)

#specific drugs in the high cLift patterns
filter(all.sets_combined_boot, LiftCrosses1==FALSE) %>% group_by(ClassCodes) %>% select(Category,items, support, count, cLift, Classes)%>% filter(ClassCodes=="BL,TE,SF")

#demonstration of cLift calculation using BL-TE-SF in MRSA (ceftaroline, omadacycline, tri-sul)
#need to calculate the conditional support--the support (%R) among isolates tested against all 3 drugs
filter(MRSA_db, !is.na(Ceftaroline), !is.na(Omadacycline), !is.na(Trimethoprim.sulfamethoxazole)) %>%
  summarise(Ceft=mean(Ceftaroline), 
            Omad=mean(Omadacycline), 
            TMZ=mean(Trimethoprim.sulfamethoxazole), 
            expected.prev=signif(mean(Ceftaroline)*mean(Omadacycline)*mean(Trimethoprim.sulfamethoxazole)*100,1),
            num.isol.req=1/(mean(Ceftaroline)*mean(Omadacycline)*mean(Trimethoprim.sulfamethoxazole)))

#number of resistance patterns at drug level for MRSA and MSSA after eCSR and cLift filters
nrow(filter(all.sets_combined_boot, LiftCrosses1==FALSE, Category=="MRSA"))
nrow(filter(all.sets_combined_boot, LiftCrosses1==FALSE, Category=="MSSA"))

#subnetwork analysis of BL-FQ-MC-LC-LG by infection site
#all possible combinations of the class codes from 2 to 5 way
sub1 <- c("BL", "FQ", "MC", "LC", "LG")
sub1.combn <- c(apply(combn(sub1,2), MARGIN =2, FUN=paste, collapse=","),
                apply(combn(sub1,3), MARGIN =2, FUN=paste, collapse=","),
                apply(combn(sub1,4), MARGIN =2, FUN=paste, collapse=","),
                apply(combn(sub1,5), MARGIN =2, FUN=paste, collapse=","))
filter(all.sets_combined_boot, LiftCrosses1==FALSE, ClassCodes %in% sub1.combn, Category %in% c("SA.SST", "SA.blood", "SA.pneum")) %>% group_by(Category) %>% dplyr::summarise(n.pattern = n_distinct(items))


#oritavancin sets did not meet eCSR P-value <=0.05
summary(subset(SA2016_sets, items %in% "Oritavancin")) #minimum Pval.eCSR = 0.78


###Discussion
#number of non-MDR MRSA and MSSA phenotypes
sum(m.mdr.classtable[which(m.mdr.classtable$NumClassRes %in% c("0", "1", "2") & m.mdr.classtable$Category=="MRSA"),"NumUniqPatterns"])
sum(m.mdr.classtable[which(m.mdr.classtable$NumClassRes %in% c("0", "1", "2") & m.mdr.classtable$Category=="MSSA"),"NumUniqPatterns"])

#diversity of class phenotypes and 3-class resistance patterns in MDR MRSA and MSSA
filter(m.mdr.classtable, !(NumClassRes %in% c("0", "1", "2"))) %>% group_by(Category) %>% dplyr::summarise(MDR.uniq.pheno=sum(NumUniqPatterns))
filter(all.sets_combined_boot, LiftCrosses1==FALSE, NumClasses==3, Category %in% c("MRSA", "MSSA")) %>% group_by(Category) %>% dplyr::summarise(Num3ClassPatterns=n())

#six-class and seven-class phenotypes in MRSA
filter(m.mdr.classtable, (NumClassRes %in% c("6", "7")) & Category=="MRSA")

#the specific patterns for six-class and seven-class MRSA phenotypes
mrsa.index <- SA.db[,"Oxacillin"]==TRUE #MRSA isolate index in whole dataset
mrsa.6.7 <- filter(MDRprofiles[mrsa.index,], NumClass %in% c("6", "7"))
mrsa.6.7
mrsa.index.6.7 <- mrsa.index & MDRprofiles[, "NumClass"] %in% c("6", "7") #index of the six-class and seven-class MRSA isolates
all_AM_panels[mrsa.index.6.7] #panels tested for those isolates. comparing mdrprofiles and panels shows that it is not necessarily different panels that is driving the discovered phenotype diversity
Reduce(intersect, lapply(mrsa.6.7$Classes.sort, unlist)) #the 6-class and 7-class MRSA isolates...
#have a background of B-lactam (a), FQ (b), macrolide (c) resistance (common to all)

#false negatives explain why not all expected cross-resistances show up in a pattern 
#(e.g. beta-lactam in a pattern just represented by ertapenem or ceftaroline)
inspect(subset(MRSA_sets, items %ain% c("Ertapenem", "Clindamycin", "Moxifloxacin", "Levofloxacin", "Ciprofloxacin")))
##the ertapenem sets generally include other B-lactams (penicillin, cefatroline; note oxacillin is not included as item in MRSA_db)
filter(MRSA_db, Ertapenem==TRUE & Clindamycin==TRUE & Moxifloxacin==TRUE & Levofloxacin==TRUE & Ciprofloxacin==TRUE)
filter(SA.db, Ertapenem==TRUE) %>% select(Penicillin, Ceftaroline)
##all the isolates with ertapenem R also have penicillin R. Only some ertapenem R also have ceftaroline R
filter(SA.db, Ceftaroline==TRUE) %>% select(Penicillin, Ertapenem)
##all ceftaroline R are resistant to Penicillin and Ertapenem (excluding missing data)
summary(subset(MRSA_sets, items %ain% c("Ertapenem", "Clindamycin", "Moxifloxacin", "Levofloxacin", "Ciprofloxacin", "Penicillin")))
##eCSR P-val all > 0.05

#ceftaroline*TMS association in skin/soft tissue: count(n), cLift
filter(all.sets_combined_boot, Category=="SA.SST", items=="Ceftaroline,Trimethoprim.sulfamethoxazole")

#association sets recovered from most common resistance phenotypes
common.pheno <- names(sort(table(MDRprofiles$Profile), decreasing=TRUE)[1:20])
common.pheno <- gsub("\\s+", ",", gsub("^\\s+|\\s+$", "", common.pheno))
#common resistance phenotypes recovered exactly
recovered1 <- common.pheno %in% filter(all.sets_combined_boot, LiftCrosses1==FALSE)$items
#common resistance phenotypes recovered from MRSA, excluding oxacillin because that was excluded from MRSA sets but all MRSA isolates are oxacillin-resistant
recovered2 <- str_remove(common.pheno, ",Oxacillin") %in% filter(all.sets_combined_boot, Category=="MRSA", LiftCrosses1==FALSE)$items
#recovered subsets from common resistance phenotypes
full <- str_split(common.pheno, ",") #need pattern as vector with each element = 1 AM
sets <- str_split(filter(all.sets_combined_boot, LiftCrosses1==FALSE)$items, ",") #need sets as vector with each element = 1 AM
recovered3 <- NULL
for (i in 1:length(common.pheno)){ #for each full phenotype
compare <- sapply(sets, function(x) x %in% full[[i]]) #compare elements of full phenotype to each set
  recovered3[i] <- any(sapply(compare, sum)>=2) #for each set, calculate number of AM from the full phenotype. do any sets contain >=2 AM (a partial of the full phenotype)
}
#recovered3 becomes the third column of Table 2

#table of top MDR profiles (typical MDR analysis) = table 2
table2 <- MDRprofiles$Profile
table2[nchar(table2)==0] <- "Pan-Susceptible"
table2 <- cbind(as.data.frame(sort(table(table2), decreasing=TRUE)[1:20]),recovered3) #bind to recovered 3
colnames(table2) <- c("Resistance_Phenotype", "n", "Recovered_Resistance_Pattern") #name columns
table2$Recovered_Resistance_Pattern <- dplyr::recode(as.character(table2$Recovered_Resistance_Pattern), "TRUE"="Yes", "FALSE"="No") #relabel recovered3 as Yes/No
table2$Recovered_Resistance_Pattern[!str_detect(table2$Resistance_Phenotype, " ")] <- "*" #if pattern is <2 antimicrobials, pattern cannot be recovered by mining
write.xlsx(table2, "results/Table 2_most common MDR phenotypes.xlsx", row.names=FALSE)

#number of common resistance phenotypes
length(recovered3)
#number of common resistance phenotypes with at least a partial set recovered
sum(recovered3)

#erythromycin*penicillin sets
#number of common phenotypes with erythromycin*penicillin
sum(colSums(sapply(full, function(x) c("Erythromycin", "Penicillin") %in% x))==2) #compare ery*pen to each full phenotype, find full phenotypes with both ery and pen, add them up
#cLift in MRSA, bloodstream, SST
filter(MRSA_MSSA_edges, edgeName=="Erythromycin Penicillin") %>% select(cat, cLift)
filter(type_edges, edgeName=="Erythromycin Penicillin") %>% select(cat, cLift)

#erythromycin*oxacillin trends
filter(SA_edges, edgeName=="Erythromycin Oxacillin") %>% select(cat, cLift)


######Sup Table 2, values copied into Word
#effect of eCSR and cLift filters on number of rules and itemset order
#pre-filtering
prefilter.sets_combined <- rbind(all_sets(setNames, varNames)[[1]], 
                                 all_sets(msetNames, mvarNames)[[1]], 
                                 all_sets(type.setNames, type.varNames)[[1]])
prefilter.sets_combined <- itemset_to_class_set(prefilter.sets_combined, AM_class)

prefilter.sets_combined %>% mutate(Analysis = case_when(
  Category %in% varNames ==TRUE ~ "Year",
  Category %in% mvarNames ==TRUE ~ "MRSA.v.MSSA",
  Category %in% type.varNames ==TRUE ~ "InfectionType")) %>%
  group_by(Analysis) %>% dplyr::summarise(n=n(), avg.n=n()/n_distinct(Category), 
                                   avg.order=mean(order), min.order=min(order), max.order=max(order), 
                                   avg.class=mean(as.numeric(NumClasses)), min.class=min(as.numeric(NumClasses)), max.class=max(as.numeric(NumClasses)))



#post-filtering with eCSR
all.sets_combined_boot %>% mutate(Analysis = case_when(
  Category %in% varNames ==TRUE ~ "Year",
  Category %in% mvarNames ==TRUE ~ "MRSA.v.MSSA",
  Category %in% type.varNames ==TRUE ~ "InfectionType")) %>%
  group_by(Analysis) %>% dplyr::summarise(n=n(), avg.n=n()/n_distinct(Category), 
                                   avg.order=mean(order), min.order=min(order), max.order=max(order), 
                                   avg.class=mean(as.numeric(NumClasses)), min.class=min(as.numeric(NumClasses)), max.class=max(as.numeric(NumClasses)))

#post-filtering with eCSR and cLift
all.sets_combined_boot %>% filter(LiftCrosses1==FALSE) %>% mutate(Analysis = case_when(
  Category %in% varNames ==TRUE ~ "Year",
  Category %in% mvarNames ==TRUE ~ "MRSA.v.MSSA",
  Category %in% type.varNames ==TRUE ~ "InfectionType")) %>%
  group_by(Analysis) %>% dplyr::summarise(n=n(), avg.n=n()/n_distinct(Category), 
                                   avg.order=mean(order), min.order=min(order), max.order=max(order), 
                                   avg.class=mean(as.numeric(NumClasses)), min.class=min(as.numeric(NumClasses)), max.class=max(as.numeric(NumClasses)))

### 1-10-20
filter(SA.db, Ceftaroline==TRUE, Trimethoprim.sulfamethoxazole==TRUE)
filter(SA.db, Telavancin==TRUE) %>% write.xlsx("results/Telavancin NS.xlsx")

#1-29-20
#TMS-tetracycline resistance in MRSA
filter(all.sets_combined_boot, Category=="MRSA", items=="Omadacycline,Tetracycline,Trimethoprim.sulfamethoxazole") #none with all three
filter(all.sets_combined_boot, Category=="MRSA", items=="Tetracycline,Trimethoprim.sulfamethoxazole")
filter(all.sets_combined_boot, Category=="MRSA", items=="Omadacycline,Trimethoprim.sulfamethoxazole")

sink()