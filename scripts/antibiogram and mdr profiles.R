# ---
# title: "Antibiogram and MDR Profiles"
# author: "Casey Cazer"
# Last updated: "Sep 11, 2020"
# output: html_document
# ---


#overall AM panels tested, regardless of breakpoint availability
#not published
all_AM_panels <- abx_profile(SA, AM_col)
write.xlsx(table(all_AM_panels, SA$Study.Year), "results/SA panels all AM.xlsx", sheetName="By Year", row.names=FALSE)
write.xlsx(table(all_AM_panels, SA.db$Oxacillin), "results/SA panels all AM.xlsx", sheetName="By Oxacillin", append=TRUE, row.names=FALSE)
write.xlsx(table(all_AM_panels, SA$Infection.Type), "results/SA panels all AM.xlsx", sheetName="By Infection Type", append=TRUE, row.names=FALSE)


#create antibiogram of %R (percent resistant) by year
#note that R and NS are used interchangeably, but truly we examine NS (non-susceptible)
group.yr <- match("Study.Year", names(SA.db)) #group index
abgm <- antibiogram(SA.db, AM_col_db, group.yr)

#remove extra Study.Year column, label grouping variable as Category
abgm <- select(abgm,-Study.Year_Prevalence)
abgm$Category <- abgm$`data[, group]`

#order columns
abgm <- select(abgm, Category, Study.Year_Number_Isolates_Tested, everything(), -`data[, group]`) %>% 
          rename(`Total Isolates`=Study.Year_Number_Isolates_Tested)


#MRSA v MSSA antibiogram
group.m <- match("Oxacillin", names(MSSA.bin)) #group index
MRSA_MSSA_antibiogram <- antibiogram(SA.db, AM_col_db, group.m)

#remove extra grouping columns, label grouping variable as Category
MRSA_MSSA_antibiogram$Category <- c("MSSA", "MRSA", "Overall")
MRSA_MSSA_antibiogram <- select(MRSA_MSSA_antibiogram, Category, Oxacillin.1_Number_Isolates_Tested, everything(), -`data[, group]`, -`Oxacillin.1_Prevalence`) %>%
                        rename(`Total Isolates`=Oxacillin.1_Number_Isolates_Tested)

#Infection Type antibiogram
group.type <- match("Infection.Type", names(SA.db)) #group index
type.abgm <- antibiogram(SA.db, AM_col_db, group.type)

#remove extra grouping column, label grouping variable as Category
type.abgm <- select(type.abgm, -Infection.Type_Prevalence)
type.abgm$Category <- type.abgm$`data[, group]`
type.abgm <- select(type.abgm, Category, Infection.Type_Number_Isolates_Tested, everything(), -`data[, group]`) %>%
  rename(`Total Isolates`=Infection.Type_Number_Isolates_Tested)


#combine into one antibiogram file
overall_abgm <- merge(abgm, MRSA_MSSA_antibiogram, by=intersect(names(abgm), names(MRSA_MSSA_antibiogram)), all=TRUE)
overall_abgm <- merge(overall_abgm, type.abgm, by=intersect(names(overall_abgm), names(type.abgm)), all=TRUE)

#sort to put by year, then MRSA/MSSA, then infection type
order <- c(levels(as.factor(SA.db$Study.Year)), "MRSA", "MSSA", levels(SA.db$Infection.Type), "Overall")
overall_abgm <- slice(overall_abgm, match(order, Category))

#round resisance prevalences
overall_abgm_round <- overall_abgm
overall_abgm_round[,seq(4,ncol(overall_abgm_round),2)] <- signif(overall_abgm_round[,seq(4,ncol(overall_abgm_round),2)], digits=2) #rounds prevalence columns to 2 significant digits


#write nice text table for manuscript with: percent (Num tested)
SupTable1 <- overall_abgm_round[,1:2] #start with category and total isolates
#collapse columns of 'number tested' and 'percent resistant' into one column for each AM
for (i in seq(3,(ncol(overall_abgm_round)-1),2)){ #for every other column in the overall_abgm (each AM)
  tab <- apply(overall_abgm_round[,c(i+1,i)],1, paste, collapse=" (") #paste the column of Prevalence with column of Number isolates tested; put N tested in parentheses
  tab <- as.data.frame(paste(tab, ")", sep="")) #close parentheses
  colnames(tab) <- str_split(colnames(overall_abgm_round)[i], "_")[[1]][1] #make column name the AM
  SupTable1 <- cbind(SupTable1, name=tab) #add to Table1
  rm(tab)
}

#NaN to '-' for %R when no isolates tested
SupTable1 <- sapply(SupTable1, function(x) {gsub("NaN", "-", x)})

filename <- "results/SupTable1_Antibiogram.xlsx"
#save to excel spreadsheet
write.xlsx(SupTable1, filename, row.names=FALSE)
#in results, the overall row of Sup Table 1 is used to report susceptibility and n isolates tested

#Overall antibiogram plus class and breakpoint for manuscript
step1 <- mutate(AM_class, AM=as.character(AM), Abbreviation=as.character(Abbreviation)) %>% left_join(SA.bp, by=c("AM"="Antimicrobial", "Abbreviation"="Abbreviation")) %>% select(AM, Abbreviation, NSbp, Resource, Code) #get AM, Abbreviation, NSbp, Resource of NSBP and Class Code
step2 <- filter(as.data.frame(SupTable1), Category=="Overall") %>% select(-Category, -"Total Isolates") %>% t() %>% as.data.frame() %>% tibble::rownames_to_column(var="Antimicrobial") #get overall antibiogram (prevalence, num tested)
Table1 <- full_join(step1, as.data.frame(step2), by=c("AM"="Antimicrobial")) %>% dplyr::rename(PropNS=V1, Class=Code) %>% arrange(Class) #combine

#save
filename <- "results/Table1_Antibiogram.xlsx"
write.xlsx(Table1, filename, row.names=FALSE)
rm(step1, step2)

########################################################################
#examine MDR profiles by number of drug resistances and class resistances
#abx profile (drugs tested) and mdr profile (drugs resistant) use interpretated MICs (binary R/S)--hence drugs without bp will be excluded


#first for all isolates
MDRprofiles <- mdr_profile(SA.db, AM_col_db, AM_class) #generates resistance profile for each isolate

#for each number of drug resistances, calculate number of isolates and number of unique drug phenotypes
mdr.drugtable <- as.data.frame(table(MDRprofiles$NumRes)) #resistance pattern lengths and number of isolates
names(mdr.drugtable) <- c("NumDrugRes", "NumIsolates")
#number of unique res patterns by num of resistances involved
for (i in 0:max(MDRprofiles$NumRes)){ #for each resistance pattern length starting at 0 (row 1)
  mdr.drugtable[i+1,3] <- nrow(unique(MDRprofiles[which(MDRprofiles$NumRes==i),])) #number of unique mdr profiles in each resistance length
}
names(mdr.drugtable)[3]="NumUniqPatterns"

#same thing but by number of class resistances
mdr.classtable <- as.data.frame(table(MDRprofiles$NumClass)) #number of isolates in each mdr profile by number of classes

for (i in 0:max(MDRprofiles$NumClass)){ #for class resistance pattern length
  mdr.classtable[i+1,3] <- length(unique(MDRprofiles$Classes.sort[which(MDRprofiles$NumClass==i)])) #number of unique mdr profiles in each class pattern length
}
names(mdr.classtable) <- c("NumClassRes", "NumIsolates", "NumUniqPatterns")
#note, there is an isolate (row 988) that has a one-class pattern of Oxacillin--it was not tested against penicillin


#MRSA v MSSA mdr profiles
MRSA_mdrprofile <- mdr_profile(MRSA.bin, AM_col_db, AM_class)
MSSA_mdrprofile <- mdr_profile(MSSA.bin, AM_col_db, AM_class)

#combine mdr profiles
MRSA_mdrprofile$Category <- rep("MRSA", nrow(MRSA_mdrprofile))
MSSA_mdrprofile$Category <- rep("MSSA", nrow(MSSA_mdrprofile))
MRSA_MSSA_mdrprofile <- rbind(MRSA_mdrprofile, MSSA_mdrprofile)

#analyze mdr profiles by MRSA/MSSA status
#for each number of drug resistances, calculate number of isolates and number of unique drug phenotypes
m.mdr.drugtable <- as.data.frame(table(MRSA_MSSA_mdrprofile$NumRes, MRSA_MSSA_mdrprofile$Category))
names(m.mdr.drugtable) <- c("NumDrugRes", "Category", "NumIsolates")

#also want to know number of unique res patterns by num of resistances involved
for (i in as.numeric(row.names(m.mdr.drugtable[which(m.mdr.drugtable$Category=="MRSA"),]))){ #for MRSA rows
  m.mdr.drugtable[i,4] <-  nrow(unique(MRSA_MSSA_mdrprofile[which(MRSA_MSSA_mdrprofile$NumRes==m.mdr.drugtable$NumDrugRes[i] & MRSA_MSSA_mdrprofile$Category=="MRSA"),])) #number of unique mdr profiles for given resistance length
}
for (i in as.numeric(row.names(m.mdr.drugtable[which(m.mdr.drugtable$Category=="MSSA"),]))){ #for MSSA rows
  m.mdr.drugtable[i,4] <- nrow(unique(MRSA_MSSA_mdrprofile[which(MRSA_MSSA_mdrprofile$NumRes==m.mdr.drugtable$NumDrugRes[i] & MRSA_MSSA_mdrprofile$Category=="MSSA"),])) #number of unique mdr profiles in each resistance length
}
names(m.mdr.drugtable)[4] <- "NumUniqPatterns"


#by classes
m.mdr.classtable <- as.data.frame(table(MRSA_MSSA_mdrprofile$NumClass, MRSA_MSSA_mdrprofile$Category)) #number of isolates in each mdr profile by number of classes
names(m.mdr.classtable) <- c("NumClassRes", "Category", "NumIsolates")

for (i in as.numeric(row.names(m.mdr.classtable[which(m.mdr.classtable$Category=="MRSA"),]))){ #for MRSA rows
  m.mdr.classtable[i,4] <- length(unique(MRSA_MSSA_mdrprofile$Classes.sort[which(MRSA_MSSA_mdrprofile$NumClass==m.mdr.classtable$NumClassRes[i] & MRSA_MSSA_mdrprofile$Category=="MRSA")])) #number of unique mdr profiles in each class pattern length
}

for (i in as.numeric(row.names(m.mdr.classtable[which(m.mdr.classtable$Category=="MSSA"),]))){ #for MSSA rows
  m.mdr.classtable[i,4] <- length(unique(MRSA_MSSA_mdrprofile$Classes.sort[which(MRSA_MSSA_mdrprofile$NumClass==m.mdr.classtable$NumClassRes[i] & MRSA_MSSA_mdrprofile$Category=="MSSA")])) #number of unique mdr profiles in each class pattern length
}
names(m.mdr.classtable)[4] <- "NumUniqPatterns"



#expected distribution of mdr classes and patterns under null hypothesis of all resistances being independent of each other
#simulate data and account for different AM panels (e.g. missingness). 
#for the entire SA_db, generate 100 random databases with the same item frequencies as the database of interest and maintaining AM panels (missingness)

 set.seed(500) #set seed so that analysis is repeatable
 
 #100 databases
 rand_SA.db <- rep(list(SA.db[AM_col_db]),100) #copy actual AM data (exclude metadata) to maintain NA positions (maintains AM susceptibility panels)
 
 n_obs <- nrow(SA.db) #number of observations
 n_nonNA <- n_obs - colSums(is.na(SA.db[AM_col_db])) #vector of number of non-NA positions (e.g. not missing) for each antimicrobial
 items <- sapply(SA.db[AM_col_db], sum, na.rm=TRUE)/n_nonNA #prevalence of each resistance


  for (i in 1:100){ #for each of the 100 databases
   for (k in 1:length(items)){ #for each AM
     bin_dat <- rbinom(n=n_nonNA[k], size=1, prob=items[k]) #generate random binary data, length=number of non-NA observations, probability=%R
     rand_dat <- rand_SA.db[[i]][,k] #get column for AM k from database (has appropriate NA positions)
     rand_dat[!is.na(rand_dat)] <- bin_dat #replace non-NA values with random binary data
     rand_SA.db[[i]][,k] <- as.logical(rand_dat) #save back to database
     
   }
  }

 #get the MDR profile for each isolate
 rand_MDRProfiles <- data.frame()
for (i in 1:100){ #for each of the 100 databases
 x <- mdr_profile(rand_SA.db[[i]], seq(1:ncol(rand_SA.db[[i]])), AM_class)
 x <- mutate(x, db=rep(i, nrow(x))) #add identifier of database #
 rand_MDRProfiles <- rbind(rand_MDRProfiles, x) #bind all MDR profiles from all 100 databases together
}


#create mdr.drugtable
 #calculate number of resistances for each isolate
 rand_MDRProfiles$NumRes <- as.factor(rand_MDRProfiles$NumRes)
 
 #calculate number of unique MDR drug profiles
rand_MDRProfiles %>% 
  group_by(db, NumRes) %>% 
  dplyr::summarise(profiles=n_distinct(Profile)) %>% 
  complete(db, NumRes, fill=list(profiles=0)) -> unique_prof

#calculate number of isolates with each unique MDR drug profile
rand_MDRProfiles %>% 
  group_by(db, NumRes) %>% 
  dplyr::summarise(isolates = n()) %>% 
  complete(db, NumRes, fill=list(isolates=0)) -> count

rand_mdr.drugtable <- inner_join(count, unique_prof, by=c("db", "NumRes")) 
rm(unique_prof, count)
names(rand_mdr.drugtable) <- c("db", "NumDrugRes", "NumIsolates", "NumUniqPatterns")


#create mdr.classtable
#repeat above but for AM classes rather than individual drugs
rand_MDRProfiles$NumClass <- as.factor(rand_MDRProfiles$NumClass)
rand_MDRProfiles %>% 
  group_by(db, NumClass) %>% 
  dplyr::summarise(profiles=length(unique(Classes.sort))) %>% 
  complete(db, NumClass, fill=list(profiles=0)) -> unique_prof

rand_MDRProfiles %>% 
  group_by(db, NumClass) %>% 
  dplyr::summarise(isolates = n()) %>% 
  complete(db, NumClass, fill=list(isolates=0)) -> count

rand_mdr.classtable <- inner_join(count, unique_prof, by=c("db", "NumClass")) 
rm(unique_prof, count)
names(rand_mdr.classtable) <- c("db", "NumClassRes", "NumIsolates", "NumUniqPatterns")


#overlapping histograms of MDR profiles in simulated data vs actual observed data
#need actual and simulated data in same dataframe. actual data stored as db=101
#plot red line (actual data) last by making actual data db=101
plot_mdr.classtable <- bind_rows(rand_mdr.classtable, as.data.frame(cbind("db" = rep(101, nrow(mdr.classtable)), mdr.classtable)))

#for each db: standardize number of unique patterns in each NumClass category by the total number of unique patterns in the db
plot_mdr.classtable <- group_by(plot_mdr.classtable, db) %>% mutate(StdNumUniqPatterns = NumUniqPatterns/sum(NumUniqPatterns))

#setup for histograms
colors <- c(rep("black", 100), "red") #actual data is red; simulated data is black
names(colors) <- seq(1,101,1)
size <- c(rep(0.1,100), 1) #actual data is thicker line than simulated data
names(size) <- seq(1,101,1)
alpha <- c(rep(0.25,100), 1) #simulated data is partially transparent
names(alpha) <- seq(1,101,1)

#plot the percent of isolates at each MDR classification (by number of resistant AM classes)
SA_mdr_PercIsolates_class_plot <- ggplot(plot_mdr.classtable, aes(x=NumClassRes, group = db, color=as.character(db), size=as.character(db), alpha=as.character(db))) +
  geom_line(aes(y=NumIsolates/nrow(SA.db)))+ #all db have same number of isolates (nrow(SA.db))
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Percent of Isolates",
       x="Number of Class Resistances")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


#plot the relative frequency of unique MDR profiles at each MDR classification
SA_mdr_RelFreqPatterns_class_plot <- ggplot(plot_mdr.classtable, aes(x=NumClassRes, group = db, color=as.character(db), size=as.character(db), alpha=as.character(db))) +
  geom_line(aes(y=StdNumUniqPatterns))+
    scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(expand=c(0,0), limits=c(0,0.45))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Relative Frequency\n Number of Unique Phenotypes",
       x="Number of Class Resistances")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))



#for plots by number of drugs rather than number of classes. Not published.
#need to make NumDrugRes from factor to character, otherwise unequal factor levels
rand_mdr.drugtable$NumDrugRes <- as.character(rand_mdr.drugtable$NumDrugRes)
mdr.drugtable$NumDrugRes <- as.character(mdr.drugtable$NumDrugRes)
plot_mdr.drugtable <- bind_rows(rand_mdr.drugtable, as.data.frame(cbind("db" = rep(101, nrow(mdr.drugtable)), mdr.drugtable)))

#standardize by the total number of unique patterns in each db
plot_mdr.drugtable <- group_by(plot_mdr.drugtable, db) %>% mutate(StdNumUniqPatterns = NumUniqPatterns/sum(NumUniqPatterns))

colors2 <- c(rep("red", 13), rep("black", 1100))
size2 <- c(rep(2,13), rep(1,1100))
alpha2 <- c(rep(1,13), rep(0.25,1100))



#MRSA v MSSA; same analysis as above
#for MRSA.bin and MSSA.bin, generate 100 random databases with the same AM resistance frequencies as the database of interest and maintaining AM panels
 set.seed(500) #set seed so that function is repeatable
 rand_MRSA.bin <- rep(list(MRSA.bin[AM_col_db]),100) #copy actual data to maintain NA positions
 rand_MSSA.bin <- rep(list(MSSA.bin[AM_col_db]),100)
 n_obs_MRSA <- nrow(MRSA.bin) #number of observations
 n_obs_MSSA <- nrow(MSSA.bin)
 n_nonNA_MRSA <- n_obs_MRSA - colSums(is.na(MRSA.bin[AM_col_db])) #vector of number of non-NA positions for each antimicrobial
 n_nonNA_MSSA <- n_obs_MSSA - colSums(is.na(MSSA.bin[AM_col_db]))
 items_MRSA <- sapply(MRSA.bin[AM_col_db], sum, na.rm=TRUE)/n_nonNA_MRSA #frequency of each resistance
 items_MSSA <- sapply(MSSA.bin[AM_col_db], sum, na.rm=TRUE)/n_nonNA_MSSA 


  for (i in 1:100){ #100 databases
   for (k in 1:length(items_MRSA)){ #for each AM
     bin_dat_MRSA <- rbinom(n=n_nonNA_MRSA[k], size=1, prob=items_MRSA[k]) #generate data
     bin_dat_MSSA <- rbinom(n=n_nonNA_MSSA[k], size=1, prob=items_MSSA[k])
     rand_dat_MRSA <- rand_MRSA.bin[[i]][,k] #get column for AM k
     rand_dat_MSSA <- rand_MSSA.bin[[i]][,k] 
     rand_dat_MRSA[!is.na(rand_dat_MRSA)] <- bin_dat_MRSA #replace non-NA values with random binomial data
     rand_dat_MSSA[!is.na(rand_dat_MSSA)] <- bin_dat_MSSA 
     rand_MRSA.bin[[i]][,k] <- as.logical(rand_dat_MRSA) #save
     rand_MSSA.bin[[i]][,k] <- as.logical(rand_dat_MSSA)
   }
  }

 #identify MDR profiles in simulated data
 rand_MDRProfiles_MRSA <- data.frame()
 rand_MDRProfiles_MSSA <- data.frame()
for (i in 1:100){ #for each random database
 x <- mdr_profile(rand_MRSA.bin[[i]], seq(1:ncol(rand_MRSA.bin[[i]])), AM_class)
 y <- mdr_profile(rand_MSSA.bin[[i]], seq(1:ncol(rand_MSSA.bin[[i]])), AM_class)
 x <- mutate(x, db=rep(i, nrow(x)))
 y <- mutate(y, db=rep(i, nrow(y)))
 rand_MDRProfiles_MRSA <- rbind(rand_MDRProfiles_MRSA, x)
 rand_MDRProfiles_MSSA <- rbind(rand_MDRProfiles_MSSA, y)
}

 
 #combine mdr profiles
rand_MDRProfiles_MRSA$Category <- rep("MRSA", nrow(rand_MDRProfiles_MRSA))
rand_MDRProfiles_MSSA$Category <- rep("MSSA", nrow(rand_MDRProfiles_MSSA))
rand_MDRProfiles_MRSA_MSSA <- rbind(rand_MDRProfiles_MRSA, rand_MDRProfiles_MSSA)


#create mdr.drugtable
rand_m.mdr.drugtable <- as.data.frame(table(rand_MDRProfiles_MRSA_MSSA$db, rand_MDRProfiles_MRSA_MSSA$NumRes, rand_MDRProfiles_MRSA_MSSA$Category))
names(rand_m.mdr.drugtable) <- c("db", "NumDrugRes", "Category", "NumIsolates")

for (i in 1:nrow(rand_m.mdr.drugtable)){ #for each row
  rand_m.mdr.drugtable[i,5] <- nrow(unique(rand_MDRProfiles_MRSA_MSSA[which(
    rand_MDRProfiles_MRSA_MSSA$NumRes==rand_m.mdr.drugtable$NumDrugRes[i] & 
    rand_MDRProfiles_MRSA_MSSA$Category==rand_m.mdr.drugtable$Category[i] & 
    rand_MDRProfiles_MRSA_MSSA$db==rand_m.mdr.drugtable$db[i]),])) #number of unique mdr profiles for given resistance length, db and category
}
names(rand_m.mdr.drugtable)[5] <- "NumUniqPatterns"


#by classes
rand_m.mdr.classtable <- as.data.frame(table(rand_MDRProfiles_MRSA_MSSA$db, rand_MDRProfiles_MRSA_MSSA$NumClass, rand_MDRProfiles_MRSA_MSSA$Category)) #number of isolates in each mdr profile by number of classes
names(rand_m.mdr.classtable) <- c("db", "NumClassRes", "Category", "NumIsolates")

for (i in 1:nrow(rand_m.mdr.classtable)){ #for each row
  rand_m.mdr.classtable[i,5] <- length(unique(rand_MDRProfiles_MRSA_MSSA[which(
    rand_MDRProfiles_MRSA_MSSA$NumClass==rand_m.mdr.classtable$NumClassRes[i] & 
      rand_MDRProfiles_MRSA_MSSA$Category==rand_m.mdr.classtable$Category[i] & 
      rand_MDRProfiles_MRSA_MSSA$db==rand_m.mdr.classtable$db[i]), "Classes.sort"])) #number of unique mdr profiles in each class pattern length
  #unique returns list here, hence use length not nrow
  }

names(rand_m.mdr.classtable)[5] <- "NumUniqPatterns"

rand_m.mdr.classtable$db <- as.numeric(as.character(rand_m.mdr.classtable$db)) #factor to numeric
 

#overlapping histograms of simulated and observed data
#need actual and rand data in same dataframe. actual data stored as db=101
plot_m.mdr.classtable <- bind_rows(rand_m.mdr.classtable, as.data.frame(cbind("db" = rep(101, nrow(m.mdr.classtable)), m.mdr.classtable)))

#for each db and category (MRSA/MSSA): standardize num unique patterns per NumClassRes by the total number of unique patterns in each db
plot_m.mdr.classtable <- group_by(plot_m.mdr.classtable, db, Category) %>% mutate(StdNumUniqPatterns = NumUniqPatterns/sum(NumUniqPatterns))

#histogram set-up
colors <- c(rep("black", 100), "red")
names(colors) <- seq(1,101,1)
size <- c(rep(0.1,100), 1)
names(size) <- seq(1,101,1)
alpha <- c(rep(0.25,100), 1)
names(alpha) <- seq(1,101,1)

#percent of isolates at each MDR classification
MRSA_mdr_PercIsolates_class_plot <- ggplot(plot_m.mdr.classtable[which(plot_m.mdr.classtable$Category=="MRSA"),], aes(x=NumClassRes, group = db, color=as.character(db), size=as.character(db), alpha=as.character(db))) +
  geom_line(aes(y=NumIsolates/nrow(MRSA_db)))+ #all db have same number of isolates (nrow(MRSA_db))
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Percent of Isolates",
       x="Number of Class Resistances")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))

MSSA_mdr_PercIsolates_class_plot <- ggplot(plot_m.mdr.classtable[which(plot_m.mdr.classtable$Category=="MSSA"),], aes(x=NumClassRes, group = db, color=as.character(db), size=as.character(db), alpha=as.character(db))) +
  geom_line(aes(y=NumIsolates/nrow(MSSA_db)))+ #all db have same number of isolate
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(labels=scales::percent_format(accuracy=1), expand=c(0,0), limits=c(0,0.55))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Percent of Isolates",
       x="Number of Class Resistances")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


#plot relative frequency of MDR profiles at each MDR classification
MRSA_mdr_RelFreqIsolates_class_plot <- ggplot(plot_m.mdr.classtable[which(plot_m.mdr.classtable$Category=="MRSA"),], aes(x=NumClassRes, group = db, color=as.character(db), size=as.character(db), alpha=as.character(db))) +
  geom_line(aes(y=StdNumUniqPatterns))+
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(expand=c(0,0), limits=c(0,0.45))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Relative Frequency \n Number of Unique Phenotypes",
       x="Number of Class Resistances")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))

MSSA_mdr_RelFreqIsolates_class_plot <- ggplot(plot_m.mdr.classtable[which(plot_m.mdr.classtable$Category=="MSSA"),], aes(x=NumClassRes, group = db, color=as.character(db), size=as.character(db), alpha=as.character(db))) +
  geom_line(aes(y=StdNumUniqPatterns))+
  scale_color_manual(values=colors)+
  scale_size_manual(values=size)+
  scale_alpha_manual(values=alpha)+
  scale_y_continuous(expand=c(0,0), limits=c(0,0.45))+
  scale_x_discrete(expand=c(0,0))+
  labs(y="Relative Frequency \n Number of Unique Phenotypes",
       x="Number of Class Resistances")+
  theme_light()+
  theme(legend.position = "none", plot.margin=unit(c(1,1,1,1), "cm"),
        axis.title=element_text(size=22), axis.text=element_text(size=20))


#by number of drugs, not classes. Not published.
#need to convert NumDresRes to character and db to numeric to avoid warning about unequal factor levels in bind_rows
rand_m.mdr.drugtable$NumDrugRes <- as.character(rand_m.mdr.drugtable$NumDrugRes)
m.mdr.drugtable$NumDrugRes <- as.character(m.mdr.drugtable$NumDrugRes)
rand_m.mdr.drugtable$db <- as.numeric(as.character(rand_m.mdr.drugtable$db))
plot_m.mdr.drugtable <- bind_rows(rand_m.mdr.drugtable, as.data.frame(cbind("db" = rep(101, nrow(m.mdr.drugtable)), m.mdr.drugtable)))

#standardize by the total number of unique patterns in each db
plot_m.mdr.drugtable <- group_by(plot_m.mdr.drugtable, db, Category) %>% mutate(StdNumUniqPatterns = NumUniqPatterns/sum(NumUniqPatterns))

colors2 <- c(rep("red", 13), rep("black", 1400))
size2 <- c(rep(2,13), rep(1,1400))
alpha2 <- c(rep(1,13), rep(0.25,1400))


#plots
fig_SA <- ggarrange(SA_mdr_PercIsolates_class_plot, SA_mdr_RelFreqPatterns_class_plot, nrow=1, ncol=2, labels=c("A","B"), font.label=list(size=24, face="bold")) %>% annotate_figure(left=text_grob(expression(paste("All  ", italic("S. aureus"))), rot=90, size=24, face="bold"))

fig_MRSA <- ggarrange(MRSA_mdr_PercIsolates_class_plot, MRSA_mdr_RelFreqIsolates_class_plot, nrow=1, ncol=2, labels=c("E","F"), font.label=list(size=24, face="bold")) %>% annotate_figure(left=text_grob(expression(paste("Methicillin-Resistant  ", italic("S. aureus"))), rot=90, size=24, face="bold"))

fig_MSSA<- ggarrange(MSSA_mdr_PercIsolates_class_plot, MSSA_mdr_RelFreqIsolates_class_plot, nrow=1, ncol=2, labels=c("C","D"), font.label=list(size=24, face="bold")) %>% annotate_figure(left=text_grob(expression(paste("Methicillin-Susceptible  ", italic("S. aureus"))), rot=90, size=24, face="bold"))

filename <- "figures/Fig1_Class Resistance and Pattern Distribution.png"

png(filename, width=14, height=21, units='in', res=600)
plot.new()

ggarrange(fig_SA, fig_MSSA, fig_MRSA, nrow=3)

dev.off()