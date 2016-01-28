#!/usr/bin/Rsript
# ---------------------------------------------------------#
# Filename      : DGEA.R                   #
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa #
# Description   : Differential Gene Expression Analysis    # 
# ---------------------------------------------------------#

#### TO DO LIST ####
# .Rdata
####

#### TO RUN IN CLI ####
# Rscript DGEA.R --accession GDS5093 --factor "disease.state" --popA "Dengue Hemorrhagic Fever,Convalescent" --popB "healthy control" --popname1 "Dengue" --popname2 "Normal" --topgenecount 250 --foldchange 0.3 --thresholdvalue 0.005
###

#############################################################################
#       Load necessary dependancies, if not previously installed            #
#############################################################################

# source('http://bioconductor.org/biocLite.R')
# biocLite('GEOquery')
# install.packages("argparser")
# install.packages('Cairo')
# install.packages('dendextend')
# install.packages('GEOquery')
# install.packages('ggplot2')
# install.packages('gplots')
# install.packages('rjson')
# install.packages('limma')
# install.packages('plyr')
# install.packages('RColorBrewer')
# install.packages('reshape2')

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
parser <- add_argument(parser, "--accession"    , help="input file")    # GEO Accession ID
parser <- add_argument(parser, "--factor"     , help="input file")    # Factor type to be classified by
parser <- add_argument(parser, "--popA", nargs='+', help="input file")    # GroupA - all the selected phenotypes (atleast one)
parser <- add_argument(parser, "--popB", nargs='+', help="input file")    # GroupB - all the selected phenotypes (atleast one)
parser <- add_argument(parser, "--popname1"     , help="input file")    # name for GroupA
parser <- add_argument(parser, "--popname2"     , help="input file")    # name for GroupB
parser <- add_argument(parser, "--topgenecount"   , help="input file")    # number of top genes to be used
parser <- add_argument(parser, "--foldchange"   , help="input file")    # fold change cut off
parser <- add_argument(parser, "--thresholdvalue" , help="input file")    # threshold value cut off
parser <- add_argument(parser, "--working_dir"    , help="input file")    # GEO Accession ID

# allow arguments to be run via the command line
argv <- parse_args(parser)

# output directory (CHANGE)
output.dir     <- argv$working_dir


# --------- Geo DataSet Input ------------ #
accession.id    <- argv$accession           # 'GDS5093' # GDS5092 GDS5091 GDS5088 GDS5086 GDS3795
factor.type     <- argv$factor            # 'disease.state'
population1     <- unlist(strsplit(argv$popA, ",")) # c('Dengue Hemorrhagic Fever','Convalescent')
population2     <- unlist(strsplit(argv$popB, ",")) # c('healthy control')
pop.name1       <- argv$popname1          # "Dengue"
pop.name2       <- argv$popname2          # "Normal"
pop.colour1     <- "#b71c1c"            # Red  
pop.colour2     <- "#0d47a1"            # Blue 


# --------- Volcano Plot ------------ #
no.of.top.genes <- as.numeric(argv$topgenecount)   # 250
toptable.sortby <- "p"                        # sort by p-value (default)
fold.change   <- as.numeric(argv$foldchange)       # 0.3
threshold.value <- as.numeric(argv$thresholdvalue)   # 0.005 # 0.05/no.of.top.genes -  Bonferroni cut-off


#############################################################################
#                        Testing Variables                              #
#############################################################################

#factor_type  <- 'genotype/variation'
#factor_type  <- 'development.stage'
#factor_type  <- 'infection

# --------- Geo DataSet Input ------------ #
# accession.id    <- 'GDS5093' # GDS5092 GDS5091 GDS5088 GDS5086 GDS3795
# factor.type     <- 'disease.state'
# population1     <- c('Dengue Hemorrhagic Fever','Convalescent')
# population2     <- c('healthy control')
# pop.name1       <- "Dengue"
# pop.name2       <- "Normal"
# pop.colour1     <- "#b71c1c"            # Red  
# pop.colour2     <- "#0d47a1"            # Blue 


# # --------- Volcano Plot ------------ #

# toptable.sortby <- "p"                 # sort by p-value (default)
# fold.change <- 0.3
# threshold.value <- 0.005 # 0.05/no.of.top.genes -  Bonferroni cut-off


#############################################################################
#                        GEO Input                                  #
#############################################################################

# import data sets and process into expression data
gse               <- getGEO(accession.id, GSEMatrix = TRUE)       # Load GEO data
met               <- Meta(gse)                    # Extract meta data
eset              <- GDS2eSet(gse, do.log2=TRUE)                  # Convert into ExpressionSet Object
X                 <- exprs(eset)                                  # Get Expression Data
gene.names        <- as.character(gse@dataTable@table$IDENTIFIER) # Store gene names
names(gene.names) <- rownames(X)
pClass            <- pData(eset)[factor.type]
colnames(pClass)  <- 'factor.type'
samples           <- rownames(pClass)

#############################################################################
#                        Two Population Preparation                         #
#############################################################################

# Create a data frame with the factors
expression.info  <- data.frame(pClass, Sample = samples, row.names = samples)


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
#                        Meta Data Access                               #
#############################################################################

# Obtain a vector of the possible factors
get.factors <- function(pDat, remove.list = c('sample', 'description', 'individual') ) {
    return ( colnames(pDat)[ !colnames(pDat) %in% remove.list ] )
}

# create json files
generate.geo.summary.json <- function(pDat, met, factors)  { 
    if (missing(factors)) {
        factors <- get.factors(pDat)
    }
    
    results.data <- list(Title = met$title,
                         Description = met$description[1],
                         Sample_Count = met$sample_count,
                         Feature_Count = met$feature_count,
                         Sample_Organism = met$sample_organism,
                         Platform = met$platform,
                         Platform_Technology_Type = met$platform_technology_type,
                         Sample_Type = met$sample_type,
                         Factor = lapply(pDat[factors], levels),
                         Pubmed_ID = met$pubmed_id,
                         Reference = met$ref,
                         Reference_Series = met$reference_series,
                         Update_Date = met$update_date)
    data.JSON <- toJSON( results.data )
    return ( data.JSON )
}


#############################################################################
#                        Top Table                                    #
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


filtered.toptable <- function(toptable, gene.names){
    toptable['gene'] <- gene.names[rownames(toptable)]
    
    #Create Sub Data
    X.toptable <- X[rownames(toptable),]
    rownames(X.toptable) <- gene.names[rownames(toptable)]
    
    return(X.toptable)
}


#############################################################################
#                        Graphical Representations                          #
#############################################################################

# Initial Boxplot
samples.boxplot <- function(){
    boxplot <- ggplot(data) + geom_boxplot(aes(x = Var2, y = value, colour = phenotypes)) + theme(axis.text.x = element_text(angle = 90, hjust = 1, colour = as.vector(expression.info$population.colour)), legend.position = 'right')+ labs(x = 'Samples', y = 'Expression Levels')
    # store Boxplot as an .png file in the working directory
    filename <- paste(output.dir,"boxplot.png",sep = "")
    ggsave(filename, plot=boxplot, width = 8, height = 4)
}

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
    dev.off()
}

# Adjusted p-value barplot
adj.p.val.histogram <- function(toptable){
    # store Histogram as an .png file in the working directory
    filename <- paste(output.dir,"histogram.png",sep = "")
    CairoPNG(file = filename, width = 600, height = 600)
    hist(toptable$adj.P.Val, breaks=100, col='skyblue', border='slateblue', xlab = "Adjusted p-values", main=NULL)
    dev.off()
}

# volcano plot (testing)
# volcanoplot1 <- function(toptable){
#     # store Volcano plot as an .png file in the working directory
#     filename <- paste(output.dir,"volcano.png",sep = "")
#     CairoPNG(file = filename, width = 600, height = 600)
#     with(toptable, plot(logFC, -log10(P.Value), pch=20, main="Volcano plot", xlim = c(-max(toptable$logFC)-0.1, max(toptable$logFC)+0.1),ylim = c(0, max(-log10(toptable$P.Value))+0.5)))
#     #volcanoplot(fit, coef=1, highlight=20, names=gene.names, col='steelblue', xlab='Log Fold Change',
#     #            ylab='Log Odds', pch=16, cex=0.5)
#     dev.off()
# }


                                                #Bonferroni cut-off    
volcanoplot2 <- function(toptable,fold.change, t = 0.05/length(gene.names)){
    
    # Highlight genes that have an absolute fold change > 2 and a p-value < Bonferroni cut-off
    toptable$threshold = as.factor(abs(toptable$logFC) > fold.change & toptable$P.Value < t)
    
    # Construct the plot object
    vol = ggplot(data=toptable, aes(x=toptable$logFC, y=-log10(toptable$P.Value), colour=threshold)) +
        geom_point(alpha=0.4, size=1.75)  + xlim(c(-max(toptable$logFC)-0.1, max(toptable$logFC)+0.1)) + ylim(c(0, max(-log10(toptable$P.Value))+0.5)) +
        xlab("log2 fold change") + ylab("-log10 p-value")
    filename <- paste(output.dir,"volcano.png",sep = "")
    ggsave(filename, plot=vol, height = 6, width = 6)
}

get.vol.data <- function(toptable,fold.change, t = 0.05/length(gene.names)){
    
    vol <- data.frame(cbind(round(toptable$logFC,3), round(-log10(toptable$P.Value),3)))
   # colnames(vol) <- c("logFC", "pVal")
    vol.list <- list(logFC = round(toptable$logFC,3),
                     pVal  = round(-log10(toptable$P.Value),3))
    return(vol.list)
}

# Top genes table
top.genes <- function(toptable, n){
    # store Top20 genes as a .csv file in the working directory
    filename <- paste(output.dir,"topgenes.csv",sep = "")
    write.csv(toptable[1:n,], file = filename)
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


# ------- Principal Component Analysis ----------

get.pc.data <- function(X){
    Xpca <- prcomp(t(X), scale= TRUE)
    s <- summary(Xpca)
    
    
    
    # Individual contribution of each princle component
    expVar <- s$importance[2,] * 100   # convert to %
    
    # Cumulative Variance
    cumVar <- s$importance[3,] * 100
    
    pcnames <-names(expVar)
    
    names(expVar) <- NULL
    names(cumVar) <- NULL
    
    results <- list(pcnames = pcnames,
                    expVar = expVar,
                    cumVar = cumVar)
    
    return(results)
}



#############################################################################
#                        Function Calling                                 #
#############################################################################

exportJson <- generate.geo.summary.json(pData(eset), met)
write(exportJson, 'factors.json')

# no.of.top.genes<-250
samples.boxplot()
toptable <- find.toptable(X, newPClass, toptable.sortby, no.of.top.genes, gene.names)
X.toptable <- filtered.toptable(toptable, gene.names)
heatmap(X.toptable,expression.info$population.colour)
adj.p.val.histogram(toptable)
toptable.all <- find.toptable(X, newPClass, toptable.sortby, length(gene.names) , gene.names)
volcanoplot2(toptable.all,fold.change)
top.genes(toptable,no.of.top.genes)
clustering()

filename <- paste(output.dir,"data.json",sep = "")
pcdata <- get.pc.data(X)
volcanoplot.data <- get.vol.data(toptable,fold.change)
json.list <- list(pc = pcdata, vol = volcanoplot.data)
write(toJSON(json.list), filename)
