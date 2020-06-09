
# title: "Data Prep"
# author: "Casey Cazer"
# Last Updated: "Feb 6, 2020"

#this script prepares MIC data for analysis


#load data
#contains blank cells which are interpreted as a factor level instead of missing (NA)--use na.strings
jmi <- read.table("data/Linelist SENTRY site_136.csv", header=TRUE, skipNul =TRUE, sep=",", na.strings = c(""))

#only Staph aureus (SA)
SA <- jmi[which(jmi$Organism.Code=="SA"),]
rm(jmi)


#use only AM that were tested in >1 isolate
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
SA.bp=read.table("data/Staph aureus breakpoints.csv", header=TRUE, skipNul = TRUE, sep=",", colClasses=c(rep("character", 4), "numeric", "character", "NULL" ), encoding = "UTF-8")
names(SA.bp)=col.names = c("Antimicrobial","S", "I", "R", "NSbp", "Resource")
#S breakpoint column is <= #; R breakpoint column is >= #. non-susceptible breakpoint column is > #. 


#ignore epidemiologic cut-off values (ECOFF) and non-CLSI/FDA/EUCAST cutoffs (e.g. Hetem13 for mupirocin)
SA.bp <- dplyr::filter (SA.bp, Resource!="ECOFF", Resource!="Hetem13")

#create SA database with only metadata and AM_include
#relevant metadata: Collection.Number, Isolate.Number, Objective, Study.Year, Bank.Number, Infection.Source, Infection.Type, Nosocomial, Admit.Date, Age, Gender, Culture.Date, Specimen.Type, VAP, ICU
names(SA)[1] <- "Collection.Number"
metadata_include <- c("Collection.Number", "Isolate.Number", "Objective", "Study.Year", "Bank.Number", "Infection.Source", "Infection.Type", "Nosocomial", "Admit.Date", "Age", "Gender", "Culture.Date", "Specimen.Type", "VAP", "ICU")
SA.bin <- select(SA, metadata_include, AM_include)

#convert MIC to numeric for AM with breakpoints (bp)
AM_include_with_bp <- match(SA.bp$Antimicrobial, names(SA.bin)) #MIC_to_interpretation function requires a numerical index: AM_include_with_bp
AM_include_with_bp <- AM_include_with_bp[!is.na(AM_include_with_bp)]
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
#exclude Oxacillin since it is used to define the groups
mdbNames <- c("MRSA_db", "MSSA_db")
oxacillin.index <- match("Oxacillin", colnames(MRSA.bin))
m.AM_col_db <- AM_col_db[!AM_col_db %in% oxacillin.index] #exclude oxacillin index
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

#code each AM as a one-letter class (for sorting) and two-letter code (for easier interpretation)

AM_class <- data.frame("AM" = c(colnames(SA.db)[AM_col_db]))
AM_class$Class <- recode(AM_class$AM, 
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
                         "Vancomycin" = "m"
)

#eventually want to use two-letter codes for AM rather than one letter codes (e.g. labels, etc)
AM_class$Code <- recode_factor(AM_class$Class, 
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
                               m = "GP"
)