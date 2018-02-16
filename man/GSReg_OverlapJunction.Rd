\name{GSReg.overlapJunction}

\alias{GSReg.overlapJunction}

\title{Generates junction overlap matrices required for SEVA}

\description{GSReg.overlapJunction generates junction overlap matrices required for SEVA. It may also perform the filtering the junctions based on the expression of the gene.}

\usage{
 GSReg.overlapJunction <- function(juncExprs,
                                  GenestoStudy=NULL,
                                  geneexpr=NULL,
                                  minmeanloggeneexp= 3,
                                  alpha =0,
                                  sparse = F,
                                  genesCoordinatesTxDB = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                  geneIDInTxDB = 'ENTREZID', 
                                  geneIDOut = 'SYMBOL',
                                  org=org.Hs.eg.db, ...) 
}

\arguments{
\item{juncExprs}{A matrix containing junction expression whose columns represent samples. Rows correspond to junctions and whose names are formed in the following format, chrN:D-A, N represents chromosome number, D represents the start coordination coordinate and A the acceptor end coordinationcoordinate. If geneexprs is not specified, overlaps and consequently SEVA statistics become independent of the quantification method for expression (assuming they do not affect the order of the junction expression) because SEVA is a rank-based method. Otherwise, junction expression and total gene expression in geneexprs (expected to have the same quantification method as junction expression) would be log transformed for applying optimal filtering.}
\item{geneexpr}{gene expression matrix whose values have been calculated using the same quantification as juncExprs. Columns must contain the same samples and same sample order as juncExprs. It is used for two types of filtering: see minmeanloggeneexp, alpha.If this parameter is missing, both filters will not be applied.}
\item{minmeanloggeneexp}{The parameter for filtering a genes based on the gene expression. If the geneexpr after transformed in log2 does not pass minmeanloggeneexp, all its corresponding junctions would be filtered. }
\item{alpha}{ The parameter for the filter junction expression . By default this filter is off by alpha=0 but we recommend alpha= 0.1 for a default analysis. In fact, if a junction maximum expression in log2 does not pass the average of its corresponding gene expression in log2 times alpha, the junction will be filtered. }
\item{sparse}{Use sparse matrices for junctions overlap. Not recommended unless you run out of memory.}
\item{genesCoordinatesTxDB}{The annotation database used for alignment. The default is hg19. It must be the annotation dataset used for gene and junction expression alignment.}
\item{geneIDInTxDB}{Gene IDs in the database.}
\item{geneIDOut}{Gene IDs used for geneexpr.}
\item{...}{Other parameters to be sent to mapIds function.}



}

\value{
\item{$Rest}{a list for all genestoStudy. Each of them contains a square matrix whose rows and columns corresponds to the genes junction and the value is one if they overlap otherwise zero i.e. $Rest[["genes"]]["junci","juncj" = I(junc i and junc j overlap)}
\item{$genesJunction}{a list for all geensToStudy. Each of them contains a vector of the names of the junctions.}

}

\author{Bahman Afsari}

\keyword{Junction overlap}

\seealso{GSReg.GeneSets.DIRAC,GSReg.GeneSets.EVA}

\examples{


library(GSReg)
require('Homo.sapiens')
require('org.Hs.eg.db')
require('GenomicRanges')

data(juncExprsSimulated)


overlapMat <- GSReg.overlapJunction(juncExprs =  junc.RPM.Simulated,
                                    geneexpr = geneExrsGSReg)
}
