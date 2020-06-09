# ---
# title: "Circle Plots"
# author: "Casey Cazer"
# Last updated: May 4, 2020
# ---


###########SET UP ##############
#For plots using specific AM
#get coordinates for vetices
  good.coords <- as.matrix(read.xlsx("data/circle graph coordinates_reduced vertex.xlsx", sheetName = "Sheet1", colIndex=c(1:2)))
  colnames(good.coords) <- NULL

#need vertex dataframe of vertex attributes
  #include only AM that were not all susceptible
  not_all_S <- sapply(SA.db[AM_col_db], sum, na.rm=TRUE)!=0 #all susceptible wouuld have column sum of 0 in binary db
  vertices <- AM_class[not_all_S,]
  vertices$Code <- droplevels(vertices$Code) #drop unused codes

#nice labels
vertices$label <- vertices$AM
vertices$label <- dplyr::recode(vertices$label, 
                         "Trimethoprim.sulfamethoxazole"="Trimethoprim \n Sulfamethoxazole", 
                         "FusidicAcid"="Fusidic \n Acid", 
                         "QuinDalfo"="Quinupristin \n Dalfopristin")

#vertex colors
#nice colors for colorblindness
#http://mkweb.bcgsc.ca/colorblind/
#https://shared-assets.adobe.com/link/21f9419a-625f-4863-7a02-f9eb22784f45
colors <- c("#FCAE91", "#CD913C", "#46AA96",  "#F569D7", "#EFF3FF", "#08519C", "#A03764", "#5F820A", "#A05FBE", "#A50F15", "#193C37")
names(colors) <- levels(vertices$Code)
vertices$color <- colors[match(vertices$Code, names(colors))] #color by AM class

#colors in long list, needed for vertex pie charts
longcolors <- list(c(vertices$color, "#DCDCDC", "#FFFFFF")) #gray DCDCDC=susceptible, white FFFFFF=not tested, color=resistant

#edge colors
edge_colors <- brewer.pal(11,"PuOr")[c(2:5,7:10)] #cLift <1 is orange, cLift >1 is purple. all edges should have cLift>1 based on large eCSR representing positive association


#plot only edges/sets that meet eCSR P-val (those are already in bestsets and seelcted in all.sets_combined) and cLift bootstrap interval does not cross 1
keep_edges <- subset(all.sets_combined_boot, all.sets_combined_boot$LiftCrosses1==FALSE)

#average all quality measures, need index of QM columns and QM names for plot_sets_weighted function
QM.names <- c("support", "count", "csr", "lift", "eSup", "eLift", "cLift", "order", "cLiftBoot0.025", "cLiftBoot0.975")
QM.index <- match(QM.names, colnames(keep_edges))



###################### By Year
SA_keep_sets <- filter(keep_edges, Category %in% varNames) #edges in year databases

#create edge dataframe for igraph by decomposing itemsets into edges connecting pairwise combinations of items (AM)
#QM are averaged among itemsets containing the same edge within each category
SA_edges <- plot_sets_weighted(SA_keep_sets, "items", "Category", QM.index, QM.names, AM_class)

#edges colored by cLift. not continuous coloring--binned by cLift
SA_edges$cLiftBin <- cut(SA_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf)) #bins of cLift: 0-0.25; 0.25-0.5; 0.5-0.75; 0.75-1; 1-2; 2-5; 5-10; 10-Inf
names(edge_colors) <- levels(SA_edges$cLiftBin)


#make a graph for each year
for (i in seq_along(varNames)){ #for each year category
  #create igraph; columns 3:4 are the Node columns, which is all that is needed to define the edge locations
graph <- graph_from_data_frame(d=as.matrix(SA_edges[which(SA_edges$cat==varNames[i]),3:4]), vertices=vertices$AM, directed=FALSE)

#assign edge properties
E(graph)$width <- SA_edges[which(SA_edges$cat==varNames[i]),"n"]+1 #edges weighted by number of contributing sets, start weight at 2 (min n=1)
E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10. note that weight = 10 corresponds to n = 9 because of above line
E(graph)$QM <- as.character(SA_edges[which(SA_edges$cat==varNames[i]),"cLiftBin"]) #cLift of the edge, binned
E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))] #edge color based on binned cLift

#communicate the percentage of isolates tested against each AM by percent fill of the node 
#vertex is pie chart of R (colored), S (gray), NA (white)
dat <- get(dbNames[i]) #get data of R/S/NA from binary db
dat <- dat[,colnames(dat) %in% vertices$AM] #drop cols not in vertices (because all susceptible)
vertices$pieNA <- sapply(dat, function(x) sum(is.na(x))/length(x)) #percent NA
vertices$pieR <- sapply(dat, function(x) sum(x, na.rm=TRUE)/length(x)) #percent R
vertices$pieS <- 1-vertices$pieNA - vertices$pieR #percent S

#for piechart, pie values must be in one column as a list with length = length of vertex colors.
#make vector = length of number of vertex colors (longcolors)--equal to length of number of vertices plus gray and white
#contains pie values for NA and S
vertices$pieVal <- mapply(append, lapply(1:nrow(vertices), function(x) matrix(0, nrow=1, ncol=nrow(vertices))), 
                          mapply(append, as.list(vertices$pieS), as.list(vertices$pieNA), SIMPLIFY=FALSE), SIMPLIFY=FALSE)

#then fill in appropriate R value for corresponding color values
for (j in 1:nrow(vertices)){ #for each node
  vertices$pieVal[[j]][j] <- vertices$pieR[j] #fill in pieR value into the pieVal list
}
names(vertices$pieVal) <- vertices$AM

#assign vertex attributes to graph
V(graph)$label <- as.character(vertices$label)
V(graph)$label.color <- "black"
V(graph)$label.cex <- 1.25
V(graph)$frame.color <- vertices$color
V(graph)$shape <- "pie"
V(graph)$pie <- vertices$pieVal
V(graph)$pie.color <- longcolors

#title
graph_attr(graph, "main") <- str_split(varNames[i], "SA")[[1]][2] #title is category

#layout
graph$layout <- good.coords

#save plots to environment
  label <- paste("g_AM_", varNames[i], sep="")
  assign(label, graph)

#save as png to maintain proportions across graphs
#individual circle graphs are 7x8 inches
  filename <- paste('./figures/', 'FigS3_', LETTERS[i], '_', varNames[i], '_AMsets_cLift.png', sep="")
  Cairo(file=filename, type="png", units="in", width=7, height=8, pointsize=8, res=300)
  plot.new()
  par(cex.main=4)
  plot(graph)
  mtext(LETTERS[i], side=3, adj=0, cex=4) #panel labels
  dev.off()
}

#plot together
#two figures required to accomodate all years
filename <- "figures/FigS3_Year AM sets_1.png"

#arrange on grid, each circle graph is 3x3 in
Cairo(file=filename, type="png", units="in", width=9, height=6, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/FigS3_A_SA2008_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_B_SA2009_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_C_SA2010_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_D_SA2011_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_E_SA2012_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_F_SA2013_AMsets_cLift.png")),
             ncol=3, nrow=2)

dev.off()



filename <- "figures/FigS3_Year AM sets_2.png"

Cairo(file=filename, type="png", units="in", width=9, height=6, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/FigS3_G_SA2014_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_H_SA2015_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_I_SA2016_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_J_SA2017_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS3_K_SA2018_AMsets_cLift.png")),
             ncol=3, nrow=2)
dev.off()




################### MRSA v MSSA
MRSA_MSSA_keep_sets <- filter(keep_edges, Category %in% mvarNames)

MRSA_MSSA_edges <- plot_sets_weighted(MRSA_MSSA_keep_sets, "items", "Category", QM.index, QM.names, AM_class)

MRSA_MSSA_edges$cLiftBin <- cut(MRSA_MSSA_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf))


#plot
for (i in rev(seq_along(mvarNames))){ #reverse order to plot MSSA first
graph <- graph_from_data_frame(d=as.matrix(MRSA_MSSA_edges[which(MRSA_MSSA_edges$cat==mvarNames[i]),3:4]), vertices=vertices$AM, directed=FALSE)

E(graph)$width <- MRSA_MSSA_edges[which(MRSA_MSSA_edges$cat==mvarNames[i]),"n"]+1 #start weight at 2
E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10
E(graph)$QM <- as.character(MRSA_MSSA_edges[which(MRSA_MSSA_edges$cat==mvarNames[i]),"cLiftBin"])
E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))]

dat <- get(mdbNames[i])
dat <- dat[,colnames(dat) %in% vertices$AM] #drop cols not in vertices (because all susceptible)
m.vertices <- filter(vertices, AM!="Oxacillin") #oxacillin is not in binary MRSA/MSSA databases, filter it out from vertices then add back in
m.vertices$pieNA <- sapply(dat, function(x) sum(is.na(x))/length(x))
m.vertices$pieR <- sapply(dat, function(x) sum(x, na.rm=TRUE)/length(x))
m.vertices$pieS <- 1-m.vertices$pieNA - m.vertices$pieR

#add back oxacillin in appropriate column location
m.vertices <- rbind(m.vertices[1:16,], vertices[which(vertices$AM=="Oxacillin"),], m.vertices[17:22,])

#oxacillin pie values
if (i==1){ #for MRSA
  m.vertices[which(m.vertices$AM=="Oxacillin"),]$pieR=1
  m.vertices[which(m.vertices$AM=="Oxacillin"),]$pieS=0
} else{ #for MSSA
  m.vertices[which(m.vertices$AM=="Oxacillin"),]$pieR=0
  m.vertices[which(m.vertices$AM=="Oxacillin"),]$pieS=1
}

m.vertices$pieVal <- mapply(append, lapply(1:nrow(m.vertices), function(x) matrix(0, nrow=1, ncol=nrow(m.vertices))), 
                            mapply(append, as.list(m.vertices$pieS), as.list(m.vertices$pieNA), SIMPLIFY=FALSE), SIMPLIFY=FALSE)

for (j in 1:nrow(m.vertices)){
  m.vertices$pieVal[[j]][j] <- m.vertices$pieR[j]
}

names(m.vertices$pieVal) <- m.vertices$AM

#assign to graph
V(graph)$label <- as.character(m.vertices$label)
V(graph)$label.color <- "black"
V(graph)$label.cex <- 1.25
V(graph)$frame.color <- m.vertices$color
V(graph)$shape <- "pie"
V(graph)$pie <- m.vertices$pieVal
V(graph)$pie.color <- longcolors

#title
graph_attr(graph, "main") <- mvarNames[i]

#layout
graph$layout <- good.coords

#save plots
label <- paste("g_AM_", mvarNames[i], sep="")
assign(label, graph)

#save as png to maintain proportions across graphs
filename <- paste('./figures/', 'FigS1_', if(i==2){LETTERS[1]}else{LETTERS[2]}, '_',  mvarNames[i], '_AMsets_cLift.png', sep="")
Cairo(file=filename, type="png", units="in", width=7, height=8, pointsize=8, res=300)
plot.new()
par(cex.main=4)
plot(graph)
mtext(if(i==2){LETTERS[1]}else{LETTERS[2]}, side=3, adj=0, cex=4) #reverse LETTERS to put MSSA as A
dev.off()
}

#plot together
filename <- "figures/FigS1_MRSA MSSA AM sets.png"

Cairo(file=filename, type="png", units="in", width=6, height=3, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/FigS1_A_MSSA_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS1_B_MRSA_AMsets_cLift.png")),
             ncol=2, nrow=1)
dev.off()


################## Infection Type
type_keep_sets <- filter(keep_edges, Category %in% type.varNames)

type_edges <- plot_sets_weighted(type_keep_sets, "items", "Category", QM.index, QM.names, AM_class)

type_edges$cLiftBin <- cut(type_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf))


for (i in seq_along(type.varNames)){
  graph <- graph_from_data_frame(d=as.matrix(type_edges[which(type_edges$cat==type.varNames[i]),3:4]), vertices=vertices$AM, directed=FALSE)
  
  E(graph)$width <- type_edges[which(type_edges$cat==type.varNames[i]),"n"]+1 #start weight at 2
  E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10
  E(graph)$QM <- as.character(type_edges[which(type_edges$cat==type.varNames[i]),"cLiftBin"])
  E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))]
  
  #vertex is pie chart of R, S, NA
  dat <- get(type.dbNames[i])
  dat <- dat[,colnames(dat) %in% vertices$AM] #drop cols not in vertices (because all susceptible)
  type.vertices <- vertices
  type.vertices$pieNA <- sapply(dat, function(x) sum(is.na(x))/length(x))
  type.vertices$pieR <- sapply(dat, function(x) sum(x, na.rm=TRUE)/length(x))
  type.vertices$pieS <- 1-type.vertices$pieNA - type.vertices$pieR
  
  #make vector = length of number of vertices plus S, NA. 
  type.vertices$pieVal <- mapply(append, lapply(1:nrow(type.vertices), function(x) matrix(0, nrow=1, ncol=nrow(type.vertices))), 
                                 mapply(append, as.list(type.vertices$pieS), as.list(type.vertices$pieNA), SIMPLIFY=FALSE), SIMPLIFY=FALSE)
  
  #then fill in appropriate R value for corresponding color values
  for (j in 1:nrow(type.vertices)){
    type.vertices$pieVal[[j]][j] <- type.vertices$pieR[j]
  }
  names(type.vertices$pieVal) <- type.vertices$AM
  
  #assign to graph
  V(graph)$label <- as.character(type.vertices$label)
  V(graph)$label.color <- "black"
  V(graph)$label.cex <- 1.25
  V(graph)$frame.color <- type.vertices$color
  V(graph)$shape <- "pie"
  V(graph)$pie <- type.vertices$pieVal
  V(graph)$pie.color <- longcolors
  
  #title
  graph_attr(graph, "main") <- dplyr::recode(str_split(type.varNames[i], "SA.")[[1]][2],
                                      "blood" = "BSI",
                                      "inabd" = "IAI",
                                      "pneum" = "PIHP",
                                      "SST" = "SSSI")
  
  #layout
  graph$layout <- good.coords
  
  #save plots
  label <- paste("g_AM_", type.varNames[i], sep="")
  assign(label, graph)
  
  #save as png to maintain proportions across graphs
  filename <- paste('./figures/', 'FigS2_', LETTERS[i], '_', type.varNames[i], '_AMsets_cLift.png', sep="")
  Cairo(file=filename, type="png", units="in", width=7, height=8, pointsize=8, res=300)
  plot.new()
  par(cex.main=4)
  plot(graph)
  mtext(LETTERS[i], side=3, adj=0, cex=4)
  dev.off()
}

#plot together
filename <- "figures/FigS2_Infection Type AM sets.png"

Cairo(file=filename, type="png", units="in", width=6, height=6, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/FigS2_A_SA.blood_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS2_B_SA.inabd_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS2_C_SA.pneum_AMsets_cLift.png")),
             rasterGrob(readPNG("figures/FigS2_D_SA.SST_AMsets_cLift.png")),
             ncol=2, nrow=2)
dev.off()

################## Legends for AM circle plots
#must be saved with same dimensions are circle graphs to maintain proportions.
#insert into documents with figures without resizing the legend or figure images to maintain proportions
filename <- "./figures/legend_color_AMSets.png"
Cairo(file=filename, type="png", units="in", width=7, height=8, pointsize=8, res=300)
plot.new()
legend('center', legend=c("1 - 2", "2 - 5", "5 - 10", "> 10"), text.width = 0.1,
       col=edge_colors[c(5:8)], lty='solid', lwd=4, ncol=4, cex=2,
       title="Line Color: Average cLift")
dev.off()

filename <- "./figures/legend_width_AMSets.png"
Cairo(file=filename, type="png", units="in", width=7, height=8, pointsize=8, res=300)
plot.new()
legend('center', legend=c("1", "3", "6", "\u2265 9 "), text.width = 0.1,
       col=edge_colors[7], lty='solid', lwd=c(2,4,7,10), ncol=4, cex=2,
       title="Line Width: Number of Patterns")
dev.off()

#combine at same proportions as combined circle plots
filename <- "./figures/legend_combined_AmSets.png"
Cairo(file=filename, type="png", units="in", width=6, height=3, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/legend_color_AMSets.png")),
                         rasterGrob(readPNG("figures/legend_width_AMSets.png")),
                         ncol=2, nrow=1)
dev.off()