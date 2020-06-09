# ---
# title: "Circle Plots by Class"
# author: "Casey Cazer"
# Last updated: March 27, 2020
# ---

###########SET UP ##############
#For plots using AM class

#need vertex df 
#include only AM that were not all susceptible
not_all_S <- sapply(SA.db[AM_col_db], sum, na.rm=TRUE)!=0
c.vertices <- AM_class[not_all_S,]
c.vertices$Code <- droplevels(c.vertices$Code) #drop unused codes

c.vertices$AM <- c.vertices$Code
c.vertices <- unique(c.vertices)

#nice labels
c.vertices$label <- c.vertices$Code

#vertex colors
#nice colors for colorblindness
#http://mkweb.bcgsc.ca/colorblind/
#https://shared-assets.adobe.com/link/21f9419a-625f-4863-7a02-f9eb22784f45
colors <- c("#FCAE91", "#CD913C", "#46AA96",  "#F569D7", "#EFF3FF", "#08519C", "#A03764", "#5F820A", "#A05FBE", "#A50F15", "#193C37")
names(colors) <- levels(c.vertices$Code)
c.vertices$color <- colors[match(c.vertices$Code, names(colors))]

#edge colors
edge_colors <- brewer.pal(11,"PuOr")[c(2:5,7:10)] #cLift <1 is orange, cLift >1 is purple


#plot only edges that meet eCSR P-val (in bestsets) and cLift bootstrap interval does not cross 1
c.keep_edges <- subset(all.sets_combined_boot, all.sets_combined_boot$LiftCrosses1==FALSE)

#need to drop 'items' column to avoid conflict inside plot_sets_weighted. ClassCodes will be used instead
c.keep_edges <- select(c.keep_edges, -items)

#average all quality measures, need index of QM columns and QM names for plot_sets_weighted function
QM.names <- c("support", "count", "csr", "lift", "eSup", "eLift", "cLift", "order", "cLiftBoot0.025", "cLiftBoot0.975")
QM.index <- match(QM.names, colnames(c.keep_edges))


################## Legends for class circle plots
#must be saved with same dimensions are circle graphs to maintain proportions.
filename <- "./figures/legend_color_ClassSets.png"
Cairo(file=filename, type="png", units="in", width=9, height=8, pointsize=8, res=300)
plot.new()
legend('center', legend=c("1 - 2", "2 - 5", "5 - 10", "> 10"), text.width = 0.1,
       col=edge_colors[c(5:8)], lty='solid', lwd=4, ncol=4, cex=2,
       title="Line Color: Average cLift")
dev.off()

filename <- "./figures/legend_width_ClassSets.png"
Cairo(file=filename, type="png", units="in", width=9, height=8, pointsize=8, res=300)
plot.new()
legend('center', legend=c("1", "3", "6", "\u2265 9 "), text.width = 0.1,
       col=edge_colors[7], lty='solid', lwd=c(2,4,7,10), ncol=4, cex=2,
       title="Line Width: Number of Patterns")
dev.off()

#combine at same proportions as combined circle plots
filename <- "./figures/legend_combined_ClassSets.png"
Cairo(file=filename, type="png", units="in", width=6, height=3, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/legend_color_ClassSets.png")),
             rasterGrob(readPNG("figures/legend_width_ClassSets.png")),
             ncol=2, nrow=1)
dev.off()


###################### By Year
c.SA_keep_sets <- filter(c.keep_edges, Category %in% varNames)

c.SA_edges <- plot_sets_weighted(c.SA_keep_sets, "ClassCodes", "Category", QM.index, QM.names, AM_class)
#will warn NAs introduced by coercion--ok, this occurs for node classes but these are the same as node names

#if set includes only one class, NA is put as second node. Change to same class code
c.SA_edges[which(c.SA_edges$Node2 == "NA"), "Node2"] <- c.SA_edges[which(c.SA_edges$Node2=="NA"), "Node1"]

#node classes are the same as node names
c.SA_edges$Node1Class <- c.SA_edges$Node1
c.SA_edges$Node2Class <- c.SA_edges$Node2

#edges colored by cLift (same as AM circle plots)
c.SA_edges$cLiftBin <- cut(c.SA_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf))
names(edge_colors) <- levels(c.SA_edges$cLiftBin)



#make a graph for each year
for (i in seq_along(varNames)){ #for each category
graph <- graph_from_data_frame(d=as.matrix(c.SA_edges[which(c.SA_edges$cat==varNames[i]),3:4]), vertices=c.vertices$AM, directed=FALSE)

#edge attributes
E(graph)$width <- c.SA_edges[which(c.SA_edges$cat==varNames[i]),"n"]+1 #start weight at 2
E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10
E(graph)$QM <- as.character(c.SA_edges[which(c.SA_edges$cat==varNames[i]),"cLiftBin"])
E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))]

#vertex attributes
V(graph)$label <- as.character(c.vertices$label)
V(graph)$label.color <- "black"
V(graph)$label.cex <- 2
V(graph)$frame.color <- c.vertices$color
V(graph)$color <- adjustcolor(c.vertices$color, alpha.f=0.5)

#title
graph_attr(graph, "main") <- str_split(varNames[i], "SA")[[1]][2]

#layout
graph$layout <- layout_in_circle(graph, order=order(c.vertices$AM))

#save plots
label <- paste("g_class_", varNames[i], sep="")
assign(label, graph)

#save as png to maintain proportions across graphs
#class circle plots are 9x8. AM circle plots are 7x8. If class circle plots are 7x8 then the one-node loops are cut-off (e.g. edges connecting FQ-FQ)
filename <- paste('./figures/', 'Fig4_', LETTERS[i], "_", varNames[i], '_ClassSets_cLift.png', sep="")
Cairo(file=filename, type="png", units="in", width=9, height=8, pointsize=8, res=300)
plot.new()
par(cex.main=4)
plot(graph)
mtext(LETTERS[i], side=3, adj=0, cex=4)
dev.off()
}

#plot together
filename <- "figures/Fig4_Year Class sets_1.png"

#each plot allocated to 3x3 space
Cairo(file=filename, type="png", units="in", width=9, height=6, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/Fig4_A_SA2008_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_B_SA2009_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_C_SA2010_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_D_SA2011_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_E_SA2012_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_F_SA2013_ClassSets_cLift.png")),
             ncol=3, nrow=2)

dev.off()



filename <- "figures/Fig4_Year Class sets_2.png"

Cairo(file=filename, type="png", units="in", width=9, height=6, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/Fig4_G_SA2014_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_H_SA2015_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_I_SA2016_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_J_SA2017_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig4_K_SA2018_ClassSets_cLift.png")),
             ncol=3, nrow=)
dev.off()


################### MRSA v MSSA
c.MRSA_MSSA_keep_sets <- filter(c.keep_edges, Category %in% mvarNames)

c.MRSA_MSSA_edges <- plot_sets_weighted(c.MRSA_MSSA_keep_sets, "ClassCodes", "Category", QM.index, QM.names, AM_class)

#if set includes only one class, NA is put as second node. Change to same class code
c.MRSA_MSSA_edges[which(c.MRSA_MSSA_edges$Node2 == "NA"), "Node2"] <- c.MRSA_MSSA_edges[which(c.MRSA_MSSA_edges$Node2=="NA"), "Node1"]

#node classes are the same as node names
c.MRSA_MSSA_edges$Node1Class <- c.MRSA_MSSA_edges$Node1
c.MRSA_MSSA_edges$Node2Class <- c.MRSA_MSSA_edges$Node2

c.MRSA_MSSA_edges$cLiftBin <- cut(c.MRSA_MSSA_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf))


#plot
for (i in rev(seq_along(mvarNames))){ #reverse to put MSSA first
graph <- graph_from_data_frame(d=as.matrix(c.MRSA_MSSA_edges[which(c.MRSA_MSSA_edges$cat==mvarNames[i]),3:4]), vertices=c.vertices$AM, directed=FALSE)

E(graph)$width <- c.MRSA_MSSA_edges[which(c.MRSA_MSSA_edges$cat==mvarNames[i]),"n"]+1 #start weight at 2
E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10
E(graph)$QM <- as.character(c.MRSA_MSSA_edges[which(c.MRSA_MSSA_edges$cat==mvarNames[i]),"cLiftBin"])
E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))]

#assign to graph
V(graph)$label <- as.character(c.vertices$label)
V(graph)$label.color <- "black"
V(graph)$label.cex <- 2
V(graph)$frame.color <- c.vertices$color
V(graph)$color <- adjustcolor(c.vertices$color, alpha.f=0.5)

#title
graph_attr(graph, "main") <- mvarNames[i]

#layout
graph$layout <- layout_in_circle(graph, order=order(c.vertices$AM))

#save plots
label <- paste("g_Class_", mvarNames[i], sep="")
assign(label, graph)

#save as png to maintain proportions across graphs
filename <- paste('./figures/', 'Fig2_', if(i==2){LETTERS[1]}else{LETTERS[2]}, '_', mvarNames[i], '_ClassSets_cLift.png', sep="")
Cairo(file=filename, type="png", units="in", width=9, height=8, pointsize=8, res=300)
plot.new()
par(cex.main=4)
plot(graph)
mtext(if(i==2){LETTERS[1]}else{LETTERS[2]}, side=3, adj=0, cex=4) #reverse LETTERS to put MSSA as A
dev.off()
}

#plot together
filename <- "figures/Fig2_MRSA MSSA Class sets.png"

Cairo(file=filename, type="png", units="in", width=6, height=3, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/Fig2_A_MSSA_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig2_B_MRSA_ClassSets_cLift.png")),
             ncol=2, nrow=1)
dev.off()


################## Infection Type
c.type_keep_sets <- filter(c.keep_edges, Category %in% type.varNames)

c.type_edges <- plot_sets_weighted(c.type_keep_sets, "ClassCodes", "Category", QM.index, QM.names, AM_class)

#if set includes only one class, NA is put as second node. Change to same class code
c.type_edges[which(c.type_edges$Node2 == "NA"), "Node2"] <- c.type_edges[which(c.type_edges$Node2=="NA"), "Node1"]

#node classes are the same as node sames
c.type_edges$Node1Class <- c.type_edges$Node1
c.type_edges$Node2Class <- c.type_edges$Node2

c.type_edges$cLiftBin <- cut(c.type_edges$cLift, breaks=c(0,0.25,0.5,0.75,1,2,5,10,Inf))

 
for (i in seq_along(type.varNames)){
  graph <- graph_from_data_frame(d=as.matrix(c.type_edges[which(c.type_edges$cat==type.varNames[i]),3:4]), vertices=c.vertices$AM, directed=FALSE)
  
  E(graph)$width <- c.type_edges[which(c.type_edges$cat==type.varNames[i]),"n"]+1 #start weight at 2
  E(graph)$width <- replace(E(graph)$width, E(graph)$width>10, 10) #cap weight at 10
  E(graph)$QM <- as.character(c.type_edges[which(c.type_edges$cat==type.varNames[i]),"cLiftBin"])
  E(graph)$color <- edge_colors[match(E(graph)$QM, names(edge_colors))]
  
  #assign to graph
  V(graph)$label <- as.character(c.vertices$label)
  V(graph)$label.color <- "black"
  V(graph)$label.cex <- 2
  V(graph)$frame.color <- c.vertices$color
  V(graph)$color <- adjustcolor(c.vertices$color, alpha.f=0.5)
  
  #title
  graph_attr(graph, "main") <- dplyr::recode(str_split(type.varNames[i], "SA.")[[1]][2],
                                             "blood" = "BSI",
                                             "inabd" = "IAI",
                                             "pneum" = "PIHP",
                                             "SST" = "SSSI")
  
  #layout
  graph$layout <- layout_in_circle(graph, order=order(c.vertices$AM))
  
  #save plots
  label <- paste("g_Class_", type.varNames[i], sep="")
  assign(label, graph)
  
  #save as png to maintain proportions across graphs
  filename <- paste('./figures/', 'Fig3_', LETTERS[i], '_', type.varNames[i], '_ClassSets_cLift.png', sep="")
  Cairo(file=filename, type="png", units="in", width=9, height=8, pointsize=8, res=300)
  plot.new()
  par(cex.main=4)
  plot(graph)
  mtext(LETTERS[i], side=3, adj=0, cex=4)
  dev.off()
}


#plot together
filename <- "figures/Fig3_Infection Type Class sets.png"

Cairo(file=filename, type="png", units="in", width=6, height=6, pointsize=8, res=300)
plot.new()
grid.arrange(rasterGrob(readPNG("figures/Fig3_A_SA.blood_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig3_B_SA.inabd_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig3_C_SA.pneum_ClassSets_cLift.png")),
             rasterGrob(readPNG("figures/Fig3_D_SA.SST_ClassSets_cLift.png")),
             ncol=2, nrow=2)
dev.off()