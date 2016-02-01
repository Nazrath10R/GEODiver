#!/usr/bin/Rsript
# ---------------------------------------------------------
# Filename      : DGEA.R
# Authors       : IsmailM, Nazrath, Suresh, Marian, Anissa
# Description   : Differential Gene Expression Analysis
# Run           : Rscript InputGeo.R --accession GDS5093 --outputdir ~/Desktop/
# ---------------------------------------------------------

#############################################################################
#                        Gene Expression  Analysis                          #
#############################################################################

library('argparser')
library('GEOquery')

#############################################################################
#                        Command Line Arguments                             #
#############################################################################

parser <- arg_parser("Input GEO Dataset")
parser <- add_argument(parser, "--geodbpath", help="GEO Dataset full path")
parser <- add_argument(parser, "--accession", help="Accession Number of the GEO Database")
parser <- add_argument(parser, "--outrdata", help="Full path to the output rData file")
argv   <- parse_args(parser)

accession.id <- argv$accession
geodbpath    <- argv$geodbpath

#############################################################################
#                        GEO Input                                          #
#############################################################################

# import data sets and process into expression data
if (is.null(geodbpath)) {
  gse <- getGEO(accession.id, GSEMatrix = TRUE)         # Automatically Load GEO dataset
} else {
  gse <- getGEO(filename = geodbpath, GSEMatrix = TRUE) # Load data from downloaded file
}
met  <- Meta(gse)                                       # Extract meta data
eset <- GDS2eSet(gse, do.log2=TRUE)                     # Convert into ExpressionSet Object
X    <- exprs(eset)                                     # Get Expression Data

if( ! is.null(argv$outrdata)){
    save.image(file = argv$outrdata )
}