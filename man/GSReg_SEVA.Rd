\name{GSReg.SEVA}

\alias{GSReg.SEVA}

\title{Applies Splice-EVA (SEVA) algorithm}

\description{GSReg.SEVA identifies Differential Spliced genes by assigning p-value by SEVA.}

\usage{
 function(juncExprs,
		  phenoVect, 
          verbose=T,
          sparse =F, ...)
}

\arguments{
\item{juncExprs}{A matrix containing junction expression whose columns represent samples. Rows correspond to junctions and whose names are formed in the following format, chrN:D-A, N represents chromosome number, D represents the start coordination coordinate and A the acceptor end coordinationcoordinate. If geneexprs is not specified, SEVA statistics are independent of the coordinates for junction expression because they are rank-based . Otherwise, junction expression would be log transformed and apply the same quantification method as used for total gene expression in geneexprs for optimal filtering.}
\item{phenoVect}{a factor containing the labels for columns of juncExprs. }
\item{verbose}{If True, reports that if the progress of the analysis after every 100 genes are analyzed.}
\item{sparse}{Use sparse matrices for junctions overlap. Not recommended unless you run out of memory.}
\item{...}{Parameters would be passed to GSReg.OverlapJunction function.}


}

\value{
a list of analysis for all genes.
\item{$E1}{the modified variance on the pathway within the samples from levels(phenotypes)[1].}
\item{$E2}{the modified variance on the pathway within the samples from levels(phenotypes)[2].}
\item{$E12}{the modified variance on the pathway across the samples from levels(phenotypes)[1] tolevels(phenotypes)[2].}
\item{$VarEta1}{the estimation of the modified variance on the pathway within the samples from levels(phenotypes)[1].}
\item{$VarEta2}{the estimation of the modified variance on the pathway within the samples from levels(phenotypes)[2].}
\item{$zscore}{zscore for the modified variance.}
\item{$zscoreD12D1}{zscoreD12D1 for the modified variance.}
\item{$zscoreD12D2}{zscoreD12D2 for the modified variance.}
\item{$pvalue}{theoretical p-value for null E1 = E2. (Standard EVA).}
\item{$pvalueD12D1}{theoretical p-value for null E1 = E12.}
\item{$pvalueD12D2}{theoretical p-value for null E2 = E12.}
\item{$pvalueTotal}{Bonferonni corrected p-value of the three p-values. }


}
\author{Bahman Afsari}

\keyword{Junction overlap}

\seealso{GSReg.OverlapJunction,GSReg.GeneSets.DIRAC,GSReg.GeneSets.EVA}

\examples{
library(GSReg)
require('Homo.sapiens')
require('org.Hs.eg.db')
require('GenomicRanges')



data(juncExprsSimulated)
SEVAjunc <- GSReg.SEVA(juncExprs = junc.RPM.Simulated,
                       phenoVect = phenotypes,
                       geneexpr = geneExrsGSReg)
                       
print(sapply(SEVAjunc,function(x) x$pvalue))
}


