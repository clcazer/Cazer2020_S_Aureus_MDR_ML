# ---
# title: "Antibiogram and MDR Profiles"
# author: "Casey Cazer"
# Last updated: "Sep 16, 2020"
# output: html_document
# ---


#overall AM panels tested, regardless of breakpoint availability
#not published
all_AM_panels <- abx_profile(SA, AM_col)
write.xlsx(table(all_AM_panels, SA$Study.Year), "results/SA panels all AM_ECV.xlsx", sheetName="By Year", row.names=FALSE)
write.xlsx(table(all_AM_panels, SA.db$Oxacillin), "results/SA panels all AM_ECV.xlsx", sheetName="By Oxacillin", append=TRUE, row.names=FALSE)
write.xlsx(table(all_AM_panels, SA$Infection.Type), "results/SA panels all AM_ECV.xlsx", sheetName="By Infection Type", append=TRUE, row.names=FALSE)


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
SupTable4 <- overall_abgm_round[,1:2] #start with category and total isolates
#collapse columns of 'number tested' and 'percent resistant' into one column for each AM
for (i in seq(3,(ncol(overall_abgm_round)-1),2)){ #for every other column in the overall_abgm (each AM)
  tab <- apply(overall_abgm_round[,c(i+1,i)],1, paste, collapse=" (") #paste the column of Prevalence with column of Number isolates tested; put N tested in parentheses
  tab <- as.data.frame(paste(tab, ")", sep="")) #close parentheses
  colnames(tab) <- str_split(colnames(overall_abgm_round)[i], "_")[[1]][1] #make column name the AM
  SupTable4 <- cbind(SupTable4, name=tab) #add to Table1
  rm(tab)
}

#NaN to '-' for %R when no isolates tested
SupTable4 <- sapply(SupTable4, function(x) {gsub("NaN", "-", x)})

#add breakpoint to table
add_bp <- mutate(AM_class, AM=as.character(AM)) %>% left_join(SA.bp, by=c("AM"="Antimicrobial")) %>% select(AM, NSbp) %>% spread(AM, NSbp) %>% mutate(Category = "ECV", 'Total Isolates' = "-") #get AM, NSbp
SupTable4 <- merge(add_bp, SupTable4, by=colnames(SupTable4), all=T, sort=F)

filename <- "results/SupTable4_Antibiogram_ECV.xlsx"
#save to excel spreadsheet
write.xlsx(SupTable4, filename, row.names=FALSE)
#in results, the overall row of Sup Table 1 is used to report susceptibility and n isolates tested