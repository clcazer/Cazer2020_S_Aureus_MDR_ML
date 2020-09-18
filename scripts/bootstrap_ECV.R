# ---
# title: "bootstrap QM"
# author: "Casey Cazer"
# Last updated: Sep 16, 2020
# ---

#create bootstrap intervals on cLift
QM <- "cLift"

bootNames <- c("SA2008_bestsets_boot", "SA2009_bestsets_boot", "SA2010_bestsets_boot", "SA2011_bestsets_boot", "SA2012_bestsets_boot", "SA2013_bestsets_boot", "SA2014_bestsets_boot", "SA2015_bestsets_boot", "SA2016_bestsets_boot", "SA2017_bestsets_boot", "SA2018_bestsets_boot")

for (i in seq_along(bootNames)){ #for each db
  bestsets <- get(best_setNames[i]) #get the bestsets
  setList <- LIST(items(bestsets)) #list all the sets in the bestsets
  dbName <- dbNames[i] #get the relevant binary db
  transName <- varNames[i] #get the relevant transaction db
  boot <- QM_bootstrap_CI(dbName, transName, setList, 1000) #bootstrap procedure, rep=1000. returns a list (one entry per set in setList); each list entry is a datagrame of bootstrap values for all enabled QM (sup, eSup, CSR, eCSR, lift, eLift, cLift)
  bootstrapCI <- ldply(boot, function(x) quantile(pull(x,QM), c(0.025,0.975), names=FALSE)) #for each set, returns the 0.025 and 0.975 quantiles of the QM. saves as dataframe
  names(bootstrapCI) <- c("items", paste(QM, "Boot0.025", sep=""), paste(QM,"Boot0.975", sep=""))
  label <- paste('SA',as.character(paste(2007+i)), "_bestsets_boot", sep="")
  assign(label,bootstrapCI) #save
  rm(boot, bootstrapCI, label, bestsets, setList, dbName, transName)
  print(i)
}



#merge sets into one dataframe
bestsets_boot_lift <- data.frame()
for (i in seq_along(bootNames)){
  sets_boot <- get(bootNames[i])
  sets_boot$Category <- rep(varNames[i], nrow(sets_boot)) #label with category
  bestsets_boot_lift <- rbind(bestsets_boot_lift, sets_boot) #combine
  rm(sets_boot)
}



###########################
#MRSA v MSSA
mbootNames <- c("MRSA_bestsets_boot", "MSSA_bestsets_boot")

for (i in seq_along(mbootNames)){
  bestsets <- get(mbest_setNames[i])
  setList <- LIST(items(bestsets))
  dbName <- mdbNames[i]
  transName <- mvarNames[i]
  boot <- QM_bootstrap_CI(dbName, transName, setList, 1000)
  bootstrapCI <- ldply(boot, function(x) quantile(pull(x,QM), c(0.025,0.975), names=FALSE)) 
  names(bootstrapCI) <- c("items", paste(QM, "Boot0.025", sep=""), paste(QM,"Boot0.975", sep=""))
  label <- mbootNames[i]
  assign(label,bootstrapCI)
  rm(boot, bootstrapCI, label, bestsets, setList, dbName, transName)
}

#merge sets into one dataframe
mbestsets_boot_lift <- data.frame()
for (i in seq_along(mbootNames)){
  sets_boot <- get(mbootNames[i])
  sets_boot$Category <- rep(mvarNames[i], nrow(sets_boot))
  mbestsets_boot_lift <- rbind(mbestsets_boot_lift, sets_boot)
  rm(sets_boot)
}


#Infection type
type.bootNames <- c("SA.blood_bestsets_boot", "SA.inabd_bestsets_boot", "SA.pneum_bestsets_boot", "SA.SST_bestsets_boot")

for (i in seq_along(type.bootNames)){
  bestsets <- get(type.best_setNames[i])
  setList <- LIST(items(bestsets))
  dbName <- type.dbNames[i]
  transName <- type.varNames[i]
  boot <- QM_bootstrap_CI(dbName, transName, setList, 1000)
  bootstrapCI <- ldply(boot, function(x) quantile(pull(x,QM), c(0.025,0.975), names=FALSE)) 
  names(bootstrapCI) <- c("items", paste(QM, "Boot0.025", sep=""), paste(QM,"Boot0.975", sep=""))
  label <- paste('SA.', type.abbrev[i] , "_bestsets_boot", sep="")
  assign(label,bootstrapCI)
  rm(boot, bootstrapCI, label, bestsets, setList, dbName, transName)
}

#merge sets into one dataframe
type.bestsets_boot_lift <- data.frame()
for (i in seq_along(type.bootNames)){
  sets_boot <- get(type.bootNames[i])
  sets_boot$Category <- rep(type.varNames[i], nrow(sets_boot))
  type.bestsets_boot_lift <- rbind(type.bestsets_boot_lift, sets_boot)
  rm(sets_boot)
}

#save all in one
bestsets_boot_lift_combined <- rbind(bestsets_boot_lift, mbestsets_boot_lift, type.bestsets_boot_lift)

#append boot columns to all_sets_combined
#first strip {} from all_sets_combined in order to match by 'items'
all.sets_combined$items <- str_replace_all(all.sets_combined$items,c("\\{|\\}"),"")

all.sets_combined_boot <- merge(all.sets_combined, bestsets_boot_lift_combined,
                                by=c("Category", "items"), all=T)

#mark those where lift interval crosses 1, the null value of no association
all.sets_combined_boot$LiftCrosses1 <- all.sets_combined_boot$cLiftBoot0.025 <=1 & all.sets_combined_boot$cLiftBoot0.975 >=1

#for the sup table, drop those with Lift Crossing 1 and drop that column, Classes column, and eLift (single letter codes)
SupTable5 <- filter(all.sets_combined_boot, LiftCrosses1==FALSE) %>% select(-LiftCrosses1, -Classes, -eLift)
#change category abbrevitions to match manuscript
##strip SA from beginning
SupTable5$Category <- as.character(SupTable5$Category)
SupTable5$Category <- str_replace(SupTable5$Category, "^SA", "")
SupTable5$Category <- dplyr::recode(SupTable5$Category, .pneum="PIHP",
                                                  .blood="BSI",
                                                  .inabd="IAI",
                                                  .SST="SSSI")
#reorder to put class columns earlier
SupTable5 <- select(SupTable5, "Category", "items", "order", "ClassCodes", "NumClasses", everything())

#name QM columns consistent with manuscript
SupTable5 <- dplyr::rename(SupTable5, Resistance_Pattern=items,
                    Number_of_Resistance_Traits_in_Pattern=order,
                    Antimicrobial_Classes_in_Pattern=ClassCodes,
                    Support=support,
                    Number_of_Isolates=count,
                    CSR=csr,
                    Lift=lift,
                    eSupport=eSup,
                    Number_of_Antimicrobial_Classes=NumClasses,
                    cLift_Bootstrap_Lower=cLiftBoot0.025,
                    cLift_Bootstrap_Upper=cLiftBoot0.975)


#save to table
filename <- "results/Sup Table 5_all resistance patterns_ECV.xlsx"
write.xlsx(SupTable5, filename, row.names = FALSE)
