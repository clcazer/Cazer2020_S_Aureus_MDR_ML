# title: "Methods Example"
# author: "Casey Cazer"
# Last updated: May 4, 2020

#for Supplementary
sink("results/supplementary methods.rtf", append=FALSE)
#example data
example <- data.frame("Ampicillin" = c("TRUE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE", "FALSE", "FALSE", "TRUE", "TRUE"),
                      "Ciprofloxacin" = c("TRUE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE", "TRUE", "FALSE", "FALSE", "FALSE"),
                      "Azithromycin" = c("TRUE", "FALSE", "FALSE", "TRUE", "TRUE", "TRUE", "FALSE", "TRUE", "FALSE", "FALSE"),
                      "Tetracycline" = c("FALSE", "FALSE", "FALSE", "TRUE", "FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE"))

example <- as.data.frame(sapply(example, as.logical))
example.trans <- as(example, "transactions")
summary(example.trans)
example #Supplementary Table6
summary(example) #Supplementary Table6 Support

#itemsets
example.sets <- apriori(example.trans, parameter=list(support=1/length(example.trans),  maxlen=4, minlen=2, target="frequent itemsets"))
example.sets@quality$CrossSupRatio <- interestMeasure(example.sets, "crossSupportRatio", example.trans, reuse=TRUE)
example.sets@quality$lift <- interestMeasure(example.sets, "lift", example.trans, reuse=TRUE)
inspect(example.sets) #Supplementary Table7


#introduce missing data
example2 <- data.frame("Ampicillin" = c("TRUE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE", "FALSE", "FALSE", "TRUE", "TRUE"),
                       "Ciprofloxacin" = c("TRUE", "NA", "FALSE", "FALSE", "TRUE", "TRUE", "NA", "FALSE", "FALSE", "FALSE"),
                       "Azithromycin" = c("TRUE", "FALSE", "NA", "NA", "TRUE", "TRUE", "FALSE", "TRUE", "NA", "FALSE"),
                       "Tetracycline" = c("FALSE", "FALSE", "FALSE", "TRUE", "FALSE", "FALSE", "FALSE", "FALSE", "TRUE", "TRUE"))

example2 <- as.data.frame(sapply(example2, as.logical))
example2.trans <- as(example2, "transactions")
summary(example2.trans)
example2 #Supplementary Table8
summary(example2) #Supplementary Table8 eSupport row; Support row is same as Supplementary Table6

#optimistic data
example2.o <- example2
example2.o[is.na(example2.o)] <- as.logical("TRUE")
example2.o.trans <- as(example2.o, "transactions")

#itemsets
example2.sets <- apriori(example2.trans, parameter=list(support=1/length(example2.trans),  maxlen=4, minlen=2, target="frequent itemsets"))
example2.sets@quality$CrossSupRatio <- interestMeasure(example2.sets, "crossSupportRatio", example2.trans, reuse=TRUE)
example2.sets@quality$lift <- interestMeasure(example2.sets, "lift", example2.trans, reuse=TRUE)
itemset_list <- LIST(items(example2.sets), decode = FALSE)
eQM <- expectedQM(example2.sets, example2.trans, example2.o.trans, example2, itemset_list)
example2.sets@quality$eSup <- eQM$eSup
example2.sets@quality$eCSR <- eQM$eCSR
example2.sets@quality$cLift <- eQM$cLift

inspect(example2.sets) #Supplementary Table9
sink()