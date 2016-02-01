#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa #
# Description   : Differential Gene Expression Analysis    #
# Rscript DGEA.R --accession GDS5093 --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005 --outputdir ~/Desktop/
# ---------------------------------------------------------#

#############################################################################
#                        Gene Expression  Analysis                          #
#############################################################################

library('argparser')
library('Cairo')
library('dendextend')
library('GEOquery')
library('ggplot2')
library('gplots')
library('rjson')
library('limma')
library('plyr')
library('RColorBrewer')
library('reshape2')

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")
parser <- add_argument(parser, "--dbrdata", help="Full Path to RData containing loaded GEO dataset")
parser <- add_argument(parser, "--geodbpath", help="GEO Dataset full path")
parser <- add_argument(parser, "--accession", help="Accession Number of the GEO Database")
parser <- add_argument(parser, "--outputdir", help="The output directory where graphs get saved")
parser <- add_argument(parser, "--factor", help="Factor type to be classified by")
parser <- add_argument(parser, "--popA", help="Phenotypes in Group A", nargs='+')
parser <- add_argument(parser, "--popB", help="Phenotypes in Group B", nargs='+')
parser <- add_argument(parser, "--popname1", help="name for GroupA")
parser <- add_argument(parser, "--popname2", help="name for GroupB")
parser <- add_argument(parser, "--topgenecount", help="number of top genes to be used")
parser <- add_argument(parser, "--foldchange", help="fold change cut off")
parser <- add_argument(parser, "--thresholdvalue", help="threshold value cut off")


# allow arguments to be run via the command line
argv <- parse_args(parser)


# --------- Geo DataSet Input ------------ #
factor.type <- argv$factor
population1 <- unlist(strsplit(argv$popA, ","))
population2 <- unlist(strsplit(argv$popB, ","))
pop.name1   <- argv$popname1
pop.name2   <- argv$popname2
pop.colour1 <- "#b71c1c"  # Red
pop.colour2 <- "#0d47a1"  # Blue
output.dir  <- argv$outputdir
dbrdata     <- argv$dbrdata

# --------- Volcano Plot ------------ #
no.of.top.genes <- as.numeric(argv$topgenecount)
toptable.sortby <- "p"
fold.change     <- as.numeric(argv$foldchange)
threshold.value <- as.numeric(argv$thresholdvalue)

if (file.exists(dbrdata)){
    load(file = dbrdata)
} else {
  if (is.null(argv$geodbpath)) {
      gse <- getGEO(argv$accession, GSEMatrix = TRUE)            # Automatically Load GEO dataset
  } else {
      gse <- getGEO(filename = argv$geodbpath, GSEMatrix = TRUE) # Load data from downloaded file
  }
  met  <- Meta(gse)                                              # Extract meta data
  eset <- GDS2eSet(gse, do.log2=TRUE)                            # Convert into ExpressionSet Object
  X    <- exprs(eset)                                            # Get Expression Data
}

#############################################################################
#                       Factor Selection                                    #
#############################################################################

gene.names       <- as.character(gse@dataTable@table$IDENTIFIER) # Store gene names
rownames(X)      <- gene.names
pClass           <- pData(eset)[factor.type]
colnames(pClass) <- 'factor.type'
samples          <- rownames(pClass)

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

# Create a data frame with the factors
expression.info <- data.frame(pClass, Sample = samples, row.names = samples)


# Introduce two columns to expression.info -
#   1. population - new two groups/populations
#   2. population.colour - colour for two new two groups/populations
expression.info <- within(expression.info, {
    population        = ifelse (factor.type %in% population1,'Group1', # if true
                                ifelse( factor.type %in% population2, 'Group2', NA) ) # if false
    population.colour = ifelse (factor.type %in% population1,pop.colour1, # if true
                                ifelse( factor.type %in% population2, pop.colour2, '#000000') ) # if false
})

# Convert to a factor
expression.info$population <- as.factor(expression.info$population)

data <- within(melt(X), {
    phenotypes = expression.info[Var2, 'factor.type']
})

# Remove samples that are not belongs to two populations
expression.info <- expression.info[complete.cases(expression.info),]
X <- X[,(colnames(X) %in% rownames(expression.info))]

# Created a new Phenotype class
newPClass           <- expression.info$population
names(newPClass)    <- expression.info$Sample

#############################################################################
#                        Top Table                                          #
#############################################################################

find.toptable <- function(X, newPClass, toptable.sortby, no.of.top.genes, gene.names){

    design  <- model.matrix(~0 + newPClass)

    # plots linear model for each gene and estimate fold changes and standard errors
    fit     <- lmFit(X, design)

    # set contrasts for all classes
    contrasts <- makeContrasts(contrasts="newPClassGroup1-newPClassGroup2",
                               levels = design)

    fit <- contrasts.fit(fit, contrasts)

    # empirical Bayes smoothing to standard errors
    fit <- eBayes(fit)

    # Sort.by shoudl be variable - 'p' or 'LogFC'
    toptable <- topTable(fit, sort.by= toptable.sortby, number=no.of.top.genes, genelist = gene.names)

    return(toptable)
}

#############################################################################
#                        Graphical Representations                          #
#############################################################################


# Heatmap
heatmap <- function(X, sample.colours, cv = TRUE, rv = TRUE){
    # store Heatmap as an .png file in the working directory
    filename <- paste(output.dir,"heatmap.png",sep = "")
    CairoPNG(file = filename, width = 800, height = 800, pointsize = 12)
    color_scale <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))(100)
    heatmap1 <- heatmap.2(X, col=color_scale, scale='row',
                          key=T, keysize=1,
                          dendrogram='column', density.info='none',
                          trace='none', cexCol=0.6, cexRow=0.1,
                          ColSideColors = sample.colours,
                          Colv = cv, Rowv = rv)
    legend("topright", legend = c(pop.name1,pop.name2), horiz = TRUE,
           col = c(pop.colour1,pop.colour2),cex=0.7, lwd = 7)
    dev.off()
}

# Clustering dendogram
clustering <- function(dist.method = "euclidean", clust.method = "average"){
    hc <- hclust(dist(t(X),dist.method), clust.method)
    dend <- as.dendrogram(hc)
    labels_colors(dend) <- expression.info$population.colour[order.dendrogram(dend)]
    filename <- paste(output.dir,"cluster.png",sep = "")
    CairoPNG(file = filename, width = 800, height = 800, pointsize = 12)
    plot(dend, main = "Cluster Dendrogram", xlab = "Samples")
    dev.off()
}

                                                #Bonferroni cut-off
volcanoplot <- function(toptable,fold.change, t = 0.05/length(gene.names)){
    # Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    toptable$threshold = as.factor(abs(toptable$logFC) > fold.change & toptable$P.Value < t)

    # Construct the plot object
    vol = ggplot(data=toptable, aes(x=toptable$logFC, y=-log10(toptable$P.Value), colour=threshold)) +
        geom_point(alpha=0.4, size=1.75)  + xlim(c(-max(toptable$logFC)-0.1, max(toptable$logFC)+0.1)) + ylim(c(0, max(-log10(toptable$P.Value))+0.5)) +
        xlab("log2 fold change") + ylab("-log10 p-value")
    filename <- paste(output.dir,"volcano.png",sep = "")
    ggsave(filename, plot=vol, height = 6, width = 6)
}

get.vol.data <- function(toptable){
    vol.list <- list( genes = toptable$ID,
                      logFC = round(toptable$logFC, 3),
                      pVal  = round(-log10(toptable$P.Value), 3) )
    return(vol.list)
}

# Principal Component Analysis

get.pc.data <- function(X, pc1 = "PC1", pc2 = "PC2"){
    Xpca <- prcomp(t(X), scale= TRUE)
    s <- summary(Xpca)

    # Individual contribution of each princle component
    expVar <- s$importance[2,] * 100   # convert to %

    # Cumulative Variance
    cumVar <- s$importance[3,] * 100

    pcnames <-names(expVar)

    names(expVar) <- NULL
    names(cumVar) <- NULL

    Xscores <- Xpca$x

    filename <- paste(output.dir,"pcscatterplot.png",sep = "")
    CairoPNG(file = filename, width = 800, height = 800, pointsize = 15)
    plot(Xscores[,pc1], Xscores[,pc2], xlab=pc1, ylab=pc2, pch=21, cex=0.9,
          cex.lab=0.9, cex.axis = 0.9, bty='L',bg = expression.info$population.colour)
    par(xpd=TRUE)
    legend("topright", y= -1, c(pop.name1,pop.name2),pch=21,cex=0.9,
            col = c(pop.colour1,pop.colour2), bty = "n", pt.bg =c(pop.colour1,pop.colour2))
    dev.off()

    results <- list(pcnames = pcnames,
                    expVar = expVar,
                    cumVar = cumVar)

    return(results)
}

#############################################################################
#                        Function Calling                                 #
#############################################################################

json.list <- list()
analysis.list <- c("Toptable","Volcano", "PCA","Heatmap", "Clustering")

if ("Toptable" %in% analysis.list){
    toptable <- find.toptable(X, newPClass, toptable.sortby, no.of.top.genes, gene.names)
    #Create Sub Data
    X.toptable <-X[toptable$ID,names(newPClass)]
    toptable.all <- find.toptable(X, newPClass, toptable.sortby, length(gene.names) , gene.names)
    json.list<- append(json.list,list(topgenes = toptable))
}

if ("Volcano" %in% analysis.list){
    volcanoplot(toptable.all,fold.change)
    volcanoplot.data <- get.vol.data(toptable)
    json.list<- append(json.list, list(vol = volcanoplot.data))
}

if ("PCA" %in% analysis.list){
    if(!is.null(X.toptable)){
        pcdata <- get.pc.data(X.toptable)
    }else{
        pcdata <- get.pc.data(X)
    }
   # pcdata <- get.pc.data(newX)
    json.list<- append(json.list, list(pc = pcdata))
}

if ("Heatmap" %in% analysis.list){
    heatmap(X.toptable,expression.info$population.colour)
}

if ("Clustering" %in% analysis.list){
    clustering()
}

if(length(json.list) != 0){
    filename <- paste(output.dir,"data.json",sep = "")
    write(toJSON(json.list), filename)
}
