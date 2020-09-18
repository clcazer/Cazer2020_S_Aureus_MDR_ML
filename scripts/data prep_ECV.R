
# title: "Data Prep"
# author: "Casey Cazer"
# Last Updated: "September 10, 2020"

#this script prepares MIC data for analysis

#load data
#contains blank cells which are interpreted as a factor level instead of missing (NA)--use na.strings
jmi <- read.table("data/Linelist SENTRY site_136.csv", header=TRUE, skipNul =TRUE, sep=",", na.strings = c(""))

#only Staph aureus (SA)
SA <- jmi[which(jmi$Organism.Code=="SA"),]
rm(jmi)

#antimicrobial (AM) columns, identified through examining dataset
AM_col <- c(18:85)

#summarise percent of isolates missing MIC values (pmissing)
pmissing <- function(x) round(sum(is.na(x))/length(x),2)

#tabulate the percent of missing values for each AM (antimicrobials) by Study Year
SA.tab <- ddply(SA[,AM_col], .(SA$Study.Year), colwise(pmissing))

#if all missing (pmissing=1), drop from consideration. keep only those without all missing
not.all.miss <- colSums(SA.tab)!=nrow(SA.tab)
SA.tab2 <- SA.tab[,not.all.miss, drop=FALSE]

#permit all antimicrobials that were tested in at least 2 isolates--excludes nitrofurantoin (only 1 isolate tested)
AM_include <- colnames(SA.tab2)[2:ncol(SA.tab2)] #Study year is first column
AM_include <- AM_include[!AM_include %in% "Nitrofurantoin"]


#AM resistance breakpoints
SA.bp=read.table("data/Staph aureus breakpoints.csv", header=TRUE, skipNul = TRUE, sep=",", colClasses=c("character", rep("NULL", 5), "character", "numeric", "character", "NULL" ), encoding = "UTF-8")
names(SA.bp)=col.names = c("Antimicrobial", "Abbreviation","NSbp", "Resource")
#non-susceptible breakpoint column is > #. 
#include only AM with ECV (EUCAST)
SA.bp <- filter(SA.bp, Resource=="EUCAST")

#create SA database with only metadata and AM_include
#relevant metadata: Collection.Number, Isolate.Number, Objective, Study.Year, Bank.Number, Infection.Source, Infection.Type, Nosocomial, Admit.Date, Age, Gender, Culture.Date, Specimen.Type, VAP, ICU
names(SA)[1] <- "Collection.Number"
metadata_include <- c("Collection.Number", "Isolate.Number", "Objective", "Study.Year", "Bank.Number", "Infection.Source", "Infection.Type", "Nosocomial", "Admit.Date", "Age", "Gender", "Culture.Date", "Specimen.Type", "VAP", "ICU")
SA.bin <- select(SA, metadata_include, AM_include)

#convert MIC to numeric for AM with breakpoints (bp)
AM_include_with_bp <- match(SA.bp$Antimicrobial, names(SA.bin)) #MIC_to_interpretation function requires a numerical index: AM_include_with_bp
AM_include_with_bp <- AM_include_with_bp[!is.na(AM_include_with_bp)]

#include only AM that were used in clinical bp analysis
AM_include_with_bp <- AM_include_with_bp[names(SA.bin)[AM_include_with_bp] %in% AM_with_clinical_bp]

#since the susceptibility testing was done with clinical breakpoints in mind, some dilution ranges may not cover the ECV
for (i in AM_include_with_bp){ #for each MIC column
  bp.index<-match(names(SA.bin)[i],SA.bp$Antimicrobial) #index of which AM bp should be applied to the i MIC column
  
  #Check that the dilutions tested cover the breakpoint appropriately (bp must be > smallest MIC value and < largest MIC value tested). Print warning if the bp does not fall within the tested dilutions
  min_dilution_too_big = as.numeric(str_match(as.character(SA.bin[,i]), "<=(.*)")[,2])>SA.bp$NSbp[bp.index] #TRUE if minimum dilution is greater than bp; FALSE if it is less than bp; NA if row does not contain a minimum dilution ("<=")
  #str_match splits the MIC into two columns at "<=". The original string is in the first column returned. If there is a "<=", the MIC is in the second column. if there is not a "<=", it generates an NA in both columns
  
  #report affected AM
  if (any(min_dilution_too_big==TRUE, na.rm=TRUE)){
    warning("minimum dilution greater than breakpoint for ", colnames(SA.bin)[i], ". ", sum(min_dilution_too_big==TRUE, na.rm=TRUE), " values replaced by NA \n")
  }
  
  max_dilution_too_small=as.numeric(str_match(as.character(SA.bin[,i]), ">(.*)")[,2])<SA.bp$NSbp[bp.index] #TRUE if max dilution is less than bp; FALSE if max dilution is greter than bp; NA if does not contain a max dilution (">")
  #str_match splits the MIC into two columns at ">". The original string is in the first column. If there is a ">", the MIC is in the second column. if there is not a ">", it generates an NA in both columns
  
  #report affected AM
  if (any(max_dilution_too_small==TRUE, na.rm=TRUE)){
    warning("maximum dilution less than breakpoint for ", colnames(data)[i], ". ", sum(max_dilution_too_small==TRUE, na.rm=TRUE)," values replaced by NA \n")
  }
  
  #turn NA into FALSE for indexing
  min_dilution_too_big <- replace_na(min_dilution_too_big, FALSE)
  max_dilution_too_small <- replace_na(max_dilution_too_small, FALSE)
  
  #replace MIC values with NA when it is a min dilution that is too big or a max dilution that is too small
  SA.bin[min_dilution_too_big,i] <- as.character(NA)
  SA.bin[max_dilution_too_small,i] <- as.character(NA)
}

SA.bin <- MIC_to_interpretation(SA.bin, AM_include_with_bp, SA.bp) #if no bp, column will not be logical T/F

#drop AM columns that are not logical (no breakpoints)
SA.db <- SA.bin %>% select(metadata_include,names(SA.bin[sapply(SA.bin,is.logical)]))

#AM column indices in SA.db for future applications
AM_col_db <- match(names(SA.db)[sapply(SA.db, is.logical)], names(SA.db))

#split by year
dbNames <- c("SA2008_db", "SA2009_db", "SA2010_db", "SA2011_db", "SA2012_db", "SA2013_db", "SA2014_db", "SA2015_db", "SA2016_db", "SA2017_db", "SA2018_db")
for (i in 1:length(dbNames)){
  x <- SA.db[which(SA.db$Study.Year==as.numeric(paste(2007+i))), AM_col_db]
  label <- paste('SA',as.character(paste(2007+i)), "_db", sep="")
  assign(label, x)
  rm(x)
  rm(label)
}



#MRSA (methicillin-resistant) vs MSSA (methicillin-susceptible); oxacillin used to measure methicillin resistance
MRSA.bin=SA.db[which(SA.db$Oxacillin==TRUE), ]
MSSA.bin=SA.db[which(SA.db$Oxacillin==FALSE),]

#create db with just AM col to include
mdbNames <- c("MRSA_db", "MSSA_db")
m.AM_col_db <- AM_col_db
MRSA_db <- MRSA.bin[,m.AM_col_db]
MSSA_db <- MSSA.bin[,m.AM_col_db]


#Infection type, only those with >30 isolates (concordant with CLSI antibiogram guidelines), excluding "other"
type <- c("bloodstream infection", "intra-abdominal infection",  "pneumonia in hospitalized patients", "skin/soft tissue infection")
type.abbrev <- c("blood", "inabd", "pneum", "SST")
type.dbNames <- c("SA.blood_db", "SA.inabd_db",  "SA.pneum_db", "SA.SST_db")
for (i in 1:length(type.dbNames)){
  x <- SA.db[which(SA.db$Infection.Type==type[i]), AM_col_db]
  label <- paste('SA.',type.abbrev[i], "_db", sep="")
  assign(label, x)
  rm(x)
  rm(label)
}



#for looking at multidrug resistance, need to classify drugs into classes
###AM classes
#a: beta-lactams: ceftaroline, ceftriaxone, doripenem, imipenem, meropenem, oxacillin, penicillin, piperacillin-tazobactam
#b: FQ: ciprofloxacin, levofloxacin, moxifloxacin
#c: Macrolide: azithromycin, clarithromycin, erythromycin, telithromycin
#d: lincosamide: clindamycin
#e: lipoglycopeptide: telavancin, oritavancin
#f: aminoglycoside: gentamicin
#g: tetracycline
#h: Sulfa
#i: lipopeptide: daptomycin
#j: fusidic acid
#k: oxazolidinone: linezolid, tedizolid
#l: streptogramin: quindalfo
#m: glycopeptide: teicoplanin, vancomycin
#n: mupirocin
#o: pleuromutilins: retapamulin

#code each AM as a one-letter class (for sorting) and two-letter code (for easier interpretation)

AM_class <- data.frame("AM" = c(colnames(SA.db)[AM_col_db]))
AM_class$Class <- dplyr::recode(AM_class$AM, 
                         "Cefazolin" = "a",
                         "Cefepime" = "a",
                         "Ceftaroline" = "a",
                         "Ceftazidime" = "a",
                         "Ceftazidime.avibactam" = "a",
                         "Ceftriaxone" = "a",
                         "Cefuroxime" = "a",
                         "Doripenem" = "a",
                         "Ertapenem" = "a",
                         "Imipenem" = "a",
                         "Meropenem" = "a",
                         "Oxacillin" = "a",
                         "Penicillin" = "a",
                         "Piperacillin.tazobactam" = "a",
                         "Ciprofloxacin" = "b",
                         "Levofloxacin" = "b",
                         "Moxifloxacin" = "b",
                         "Azithromycin" = "c",
                         "Clarithromycin" = "c",
                         "Erythromycin" = "c",
                         "Telithromycin" = "c",
                         "Clindamycin" = "d",
                         "Oritavancin" = "e",
                         "Telavancin" = "e",
                         "Gentamicin" = "f",
                         "Doxycycline" = "g",
                         "Minocycline" = "g", 
                         "Omadacycline" = "g",
                         "Tetracycline" = "g",
                         "Tigecycline" = "g",
                         "Trimethoprim.sulfamethoxazole" = "h",
                         "Daptomycin" = "i",
                         "FusidicAcid" = "j",
                         "Linezolid" = "k",
                         "Tedizolid" = "k",
                         "QuinDalfo" = "l",
                         "Teicoplanin" = "m",
                         "Vancomycin" = "m",
                         "Mupirocin" = "n",
                         "Retapamulin" = "o"
)

#eventually want to use two-letter codes for AM rather than one letter codes (e.g. labels, etc)
AM_class$Code <- dplyr::recode_factor(AM_class$Class, 
                               a = "BL",
                               b = "FQ",
                               c = "MC",
                               d = "LC",
                               e = "LG",
                               f = "AG",
                               g = "TE",
                               h = "SF",
                               i = "LP",
                               j = "FA",
                               k = "OZ",
                               l = "SG",
                               m = "GP",
                               n = "MP",
                               o = "PM"
)

#get abbreviation for antimicrobial from SA.bp
AM_class$Abbreviation <- with(SA.bp, Abbreviation[match(AM_class$AM, Antimicrobial)])