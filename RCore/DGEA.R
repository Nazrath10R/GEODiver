#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anisa  #
# Description   : Differential Gene Expression Analysis    #
# Rscript DGEA.R --accession GDS5093 --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent,Dengue Fever" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005 --distance "euclidean" --clustering "average" --dbrdata ~/Desktop/GDS5093.rData --outputdir ~/Desktop/ --heatmaprows 100 --dendrow TRUE --dendcol TRUE --analyse "Boxplot,Volcano,PCA,Heatmap" --expsavepath ~/Desktop/topexpr.rData
# ---------------------------------------------------------#

#############################################################################
#                        Import Libraries                                   #
#############################################################################

# silent library loading messages on command line
suppressMessages(library("limma"))
suppressMessages(library("gplots"))
suppressMessages(library("GEOquery"))
suppressMessages(library("pheatmap"))
suppressMessages(library("plyr"))
suppressMessages(library("DMwR"))

# load required libraries
library("argparser")    # Argument passing
library("Cairo")        # Plots saving
library("DMwR")         # Outlier Prediction for clustering
library("GEOquery")     # GEO dataset Retrieval
library("ggplot2")      # Graphs designing
library("gplots")       # Graphs designing
library("jsonlite")     # Convert R object to JSON format
library("pheatmap")     # Heatmap Generating
library("limma")        # Differencial Gene Expression Analysis
library("plyr")         # Splitting, Applying and Combining Data
library("RColorBrewer") # Import Colour Pallete
library("reshape2")     # Prepare dataset for ggplot


#############################################################################
#                        Command Line Arguments                             #
#############################################################################

# set parsers for all input arguments
parser <- arg_parser("This parser contains the input arguments")

# General Parameters
parser <- add_argument(parser, "--outputdir",
    help = "The outout directory where graphs get saved")
parser <- add_argument(parser, "--dbrdata",
    help = "Downloaded GEO dataset full path")
parser <- add_argument(parser, "--analyse",
    help = "List of analysis to be performed", nargs = "+")
parser <- add_argument(parser, "--geodbpath",
    help  =  "GEO Dataset full path")
parser <- add_argument(parser, "--expsavepath",
    help  =  "Toptable Expression data saving full path")

# Sample Parameters
parser <- add_argument(parser, "--accession",
    help = "Accession Number of the GEO Database")
parser <- add_argument(parser, "--factor",
    help = "Factor type to be classified by")
parser <- add_argument(parser, "--popA",
    help = "GroupA - all the selected phenotypes (atleast one)", nargs = "+")
parser <- add_argument(parser, "--popB",
    help = "GroupB - all the selected phenotypes (atleast one)", nargs = "+")
parser <- add_argument(parser, "--popname1",
    help = "name for GroupA")
parser <- add_argument(parser, "--popname2",
    help = "name for GroupB")

# Volcano plot Parameters
parser <- add_argument(parser, "--foldchange",
    help = "fold change cut off")
parser <- add_argument(parser, "--topgenecount",
    help = "number of top genes to be used")
parser <- add_argument(parser, "--thresholdvalue",
    help = "threshold value cut off")

# Heatmap
parser <- add_argument(parser, "--heatmaprows",
    help = "Number of genes show in the heatmap")
parser <- add_argument(parser, "--dendrow",
    help = "Boolean value for display dendogram for Genes")
parser <- add_argument(parser, "--dendcol",
    help = "Boolean value for display dendogram for Samples")
parser <- add_argument(parser, "--distance",
    help = "Distance measurement methods")
parser <- add_argument(parser, "--clustering",
    help = "HCA clustering methods")

# allow arguments to be run via the command line
argv   <- parse_args(parser)

#############################################################################
#                        Command Line Arguments Retrieval                   #
#############################################################################

# General Parameters
output.dir      <- argv$outputdir
dbrdata         <- argv$dbrdata
analysis.list   <- unlist(strsplit(argv$analyse, ","))

# Sample Parameters
factor.type     <- argv$factor
population1     <- unlist(strsplit(argv$popA, ","))
population2     <- unlist(strsplit(argv$popB, ","))
pop.name1       <- argv$popname1
pop.name2       <- argv$popname2
pop.colour1     <- "#b71c1c"  # Red
pop.colour2     <- "#0d47a1"  # Blue

# Volcano plot Parameters
topgene.count   <- as.numeric(argv$topgenecount)
fold.change     <- as.numeric(argv$foldchange)
threshold.value <- as.numeric(argv$thresholdvalue)
toptable.sortby <- "p"

# Clustering
if (argv$distance %in%
    c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")){
    dist.method <- argv$distance
}else{
    dist.method <- "euclidean"
}

if (argv$clustering %in%
    c("ward.D", "ward.D2", "single", "complete",
      "average", "mcquitty", "median", "centroid")){
    clust.method <- argv$clustering
}else{
    clust.method <- "average"
}

# Heatmap
heatmap.rows <- as.numeric(argv$heatmaprows)
dendrow <- as.logical(argv$dendrow)
dendcol <- as.logical(argv$dendcol)


#############################################################################
#                        Load GEO Dataset to Start Analysis                 #
#############################################################################


if (file.exists(dbrdata)){
    load(file = dbrdata)
} else {
    if (is.na(argv$geodbpath)) {
        # Load data from downloaded file
        gse <- getGEO(filename = argv$geodbpath, GSEMatrix = TRUE)
    } else {
        # Automatically Load GEO dataset
        gse <- getGEO(argv$accession, GSEMatrix = TRUE)
    }
    # Convert into ExpressionSet Object
    eset <- GDS2eSet(gse, do.log2 = FALSE)
}

X <- exprs(eset)  # Get Expression Data

# auto-detect if data needs transformation and log2 transform if needed
qx <- as.numeric(quantile(X, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = T))
logc <- (qx[5] > 100) ||
    (qx[6] - qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (logc) {
    X[which(X <= 0)] <- NaN
    exprs(eset) <- log2(X)
}


#############################################################################
#                       Factor Selection                                    #
#############################################################################

# Store gene names
X               <- exprs(eset)
gene.names      <- as.character(gse@dataTable@table$IDENTIFIER)
rownames(X)     <- gene.names

# Phenotype Selection
pclass           <- pData(eset)[factor.type]
colnames(pclass) <- "factor.type"

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

# Create a data frame with the factors
expression.info  <- data.frame(pclass,
                               Sample = rownames(pclass),
                               row.names = rownames(pclass))

# Introduce two columns to expression.info :
#   1. population - new two groups/populations
#   2. population.colour - colour for two new two groups/populations
expression.info <- within(expression.info, {
    population        <- ifelse(factor.type %in% population1, "Group1",
                         ifelse(factor.type %in% population2, "Group2", NA))
    population.colour <- ifelse(factor.type %in% population1, pop.colour1,
                         ifelse(factor.type %in% population2, pop.colour2,
                         "#000000"))
})

# Convert population column to a factor
expression.info$population <- as.factor(expression.info$population)

# Remove samples that are not belongs to two populations
expression.info <- expression.info[complete.cases(expression.info), ]
X <- X[, (colnames(X) %in% rownames(expression.info))]

# Data preparation for ggplot-Boxplot
data <- within(melt(X), {
    phenotypes <- expression.info[Var2, "factor.type"]
    Groups     <- expression.info[Var2, "population.colour"]
})

# Created a new Phenotype class
newpclass           <- expression.info$population
names(newpclass)    <- expression.info$Sample

#############################################################################
#                        Top Table                                          #
#############################################################################

find.toptable <- function(X, newpclass, toptable.sortby, topgene.count){

    design  <- model.matrix(~0 + newpclass)

    # plots linear model for each gene and
    # estimate fold changes and standard errors
    fit     <- lmFit(X, design)

    # set contrasts for two groups
    contrasts <- makeContrasts(contrasts = "newpclassGroup1-newpclassGroup2",
                               levels = design)

    fit <- contrasts.fit(fit, contrasts)

    # empirical Bayes smoothing to standard errors
    fit <- eBayes(fit)

    # Create top Table
    toptable <- topTable(fit,
                         sort.by = toptable.sortby, # 'p' or 'LogFC'
                         number = topgene.count)

    return(toptable)
}

#############################################################################
#                        Graphical Representations                          #
#############################################################################

# Boxplot
samples.boxplot <- function(data, pop.colours, pop.names, path){
    boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = Groups)) + theme(axis.text.x = element_text(angle = 70, hjust = 1), legend.position = "right")+ labs(x = "Samples", y = "Expression Levels") + scale_color_manual(name = "Groups", values = pop.colours, labels = pop.names)
    filename <- paste(path, "boxplot.png", sep = "")
    ggsave(filename, plot = boxplot, width = 8, height = 4)
}

# Heatmap
heatmap <- function(X.matix, exp, heatmap.rows = 100, dendogram.row, dendogram.col,
    dist.method, clust.method, path){

    col.pal <- colorRampPalette(rev(
               RColorBrewer::brewer.pal(11, "RdYlGn")))(100)

    # Column dendogram
    if (dendogram.col == TRUE){
        hc <- hclust(dist(t(X.matix),
                    method = dist.method),
                    method = clust.method)

        # Rank outliers using distance and clustering parameters
        o <- outliers.ranking(t(X.matix),
                              test.data   = NULL,
                              method      = "sizeDiff", # Outlier finding method
                              method.pars = NULL,
                              clus = list(dist = dist.method,
                                          alg  = "hclust",
                                          meth = clust.method))

        ann.col <- data.frame(  Population = exp[, "population"],
                                Factor     = exp[, "factor.type"],
                                Outliers   = o$prob.outliers)
        column.gap <- 0
    }else{
        hc <- FALSE

        ann.col <- data.frame( Population = exp[, "population"],
                               Factor     = exp[, "factor.type"])
        column.gap <- length( (which(ann.col[, "Population"] == "Group1") == T))

    }

    rownames(ann.col) <- exp[, "Sample"]

    filename <- paste(path, "heatmap.svg", sep = "")
    CairoSVG(file = filename)

    pheatmap(X.matix[1:heatmap.rows, ],
             cluster_row    = dendogram.row,
             cluster_cols   = hc,
             annotation_col = ann.col,
             color          = col.pal,
             fontsize       = 6.5,
             fontsize_row   = 3.5,
             fontsize_col   = 3.5,
             gaps_col       = column.gap)
    dev.off()
}

# Apply Bonferroni cut-off as the default thresold value
volcanoplot <- function(toptable, fold.change, t = 0.05 / length(gene.names), path){

    # Highlight genes that have an logFC greater than fold change
    # a p-value less than Bonferroni cut-off
    toptable$threshold <- as.factor(abs(toptable$logFC)  > fold.change &
                                        toptable$P.Value < t)

    # Construct the plot object
    vol <- ggplot(data = toptable, aes(x = toptable$logFC, y = -log10(toptable$P.Value), colour = threshold)) +
        geom_point(alpha = 0.4, size = 1.75)  + xlim(c(-max(toptable$logFC) - 0.1, max(toptable$logFC) + 0.1)) + ylim(c(0, max(-log10(toptable$P.Value)) + 0.5)) + xlab("log2 fold change") + ylab("-log10 p-value")

    # File saving
    filename <- paste(path, "volcano.png", sep = "")
    ggsave(filename, plot = vol, height = 6, width = 6)
}


get.volcanodata <- function(toptable){

    vol.list <- list( genes = toptable$ID,
                      logFC = round(toptable$logFC, 3),
                      pVal  = -log10(toptable$P.Value))
    return(vol.list)
}

# Principal Component Analysis
get.pcdata <- function(Xpca){

    s <- summary(Xpca)

    # Individual contribution of each princle component in percentages
    exp.var <- s$importance[2, ] * 100

    # Cumulative Variance in percentages
    cum.var <- s$importance[3, ] * 100

    # PC names
    pcnames <- names(exp.var)

    names(exp.var) <- NULL
    names(cum.var) <- NULL

    results <- list(pcnames = pcnames,
                    expVar = exp.var,
                    cumVar = cum.var)

    return(results)
}

get.pcplotdata <- function(Xpca, populations){

    Xscores <- Xpca$x

    rownames(Xscores) <- NULL
    cols <- colnames(Xscores)

    Xscores <- lapply(1:nrow(Xscores),
                      function(y) split(Xscores[, y], populations))
    names(Xscores) <- cols


   return(unlist(Xscores, recursive = FALSE))
}

#############################################################################
#                        Function Calling                                 #
#############################################################################

json.list <- list()

# Toptable
toptable <- find.toptable(X, newpclass, toptable.sortby, topgene.count)

# Adding to JSON file
temp.toptable <- toptable
names(temp.toptable) <- NULL
json.list <- append(json.list, list(tops = temp.toptable))

# Filter toptable data from X
X.toptable <- X[as.numeric(rownames(toptable)), ]

# save toptable expression data
if (! is.na(argv$expsavepath)){
    save(X.toptable, expression.info, file = argv$expsavepath )
}

if ("Boxplot" %in% analysis.list){
    samples.boxplot(data, c(pop.colour1, pop.colour2),
        c(pop.name1, pop.name2), path = output.dir)
}

if ("Volcano" %in% analysis.list){

    # Create complete volcano plot
    toptable.all <- find.toptable(X, newpclass,
        toptable.sortby, length(gene.names))

    volcanoplot(toptable.all, fold.change, threshold.value, output.dir)

    # save volcanoplot top data as JSON
    volcanoplot.data <- get.volcanodata(toptable)
    json.list <- append(json.list, list(vol = volcanoplot.data))
}

if ("PCA" %in% analysis.list){

    Xpca <- prcomp(t(X.toptable), scale = TRUE)

    # PC individual and cumulative values
    pcdata <- get.pcdata(Xpca)
    json.list <- append(json.list, list(pc = pcdata))

    # PC scatter plot
    pcplotdata <- get.pcplotdata(Xpca, expression.info[, "population"])

    # adding both data to the json list
    json.list <- append(json.list, list(pcdata = pcplotdata))
}

if ("Heatmap" %in% analysis.list){
    heatmap(X.toptable, expression.info, heatmap.rows = heatmap.rows,
            dendrow, dendcol, dist.method, clust.method, output.dir)
}

if (length(json.list) != 0){
    filename <- paste(output.dir, "data.json", sep = "")
    write(toJSON(json.list, digits=I(4)), filename )
}