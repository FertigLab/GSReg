\name{GSReg.GeneSets.EVA}

\alias{GSReg.GeneSets.EVA}

\title{Performs Gene Set Analysis using Expression Variation Analysis (EVA).}

\description{GSReg.GeneSets.EVA performs modified version DIRAC papers. Using a theoretical analysis, we can calculate p-value which makes extreme low p-values available.  }

\usage{
 GSReg.GeneSets.EVA(geneexpres, pathways, phenotypes,minGeneNum=5) 
}

\arguments{
\item{geneexpres}{the matrix of gene expressions. The rownames must represent gene names and the columns represent samples. There must not be any missing values. Please use imputation or remove the genes with missing values.}
\item{pathways}{a list containing pathway information. Each element represents a pathway as a character vector. The genes shown in the pathway must be present in geneexpres. geneexpres must have numeric and finite numbers.}
\item{phenotypes}{a binary factor containing the phenotypes for samples in geneexpres; hence, the column number of geneexpres and the length of phenotypes must be equal. }
\item{minGeneNum}{the minimum number of genes required in a pathway.}

}

\value{
a list of analysis for all pathways.
\item{$E1}{the modified variance on the pathway within the samples from levels(phenotypes)[1].}
\item{$E2}{the modified variance on the pathway within the samples from levels(phenotypes)[2].}
\item{$E12}{the modified variance on the pathway across the samples from levels(phenotypes)[1] tolevels(phenotypes)[2].}
\item{$VarEta1}{the estimation of the modified variance on the pathway within the samples from levels(phenotypes)[1].}
\item{$VarEta2}{the estimation of the modified variance on the pathway within the samples from levels(phenotypes)[2].}
\item{$zscore}{zscore for the modified variance.}
\item{$pvalue}{theoretical p-value for null E1 = E2. (Standard EVA).}
\item{$pvalueD12D1}{theoretical p-value for null E1 = E12.}
\item{$pvalueD12D2}{theoretical p-value for null E2 = E12.}
\item{$pvalueTotal}{Bonferonni corrected p-value of the three p-values. }

}

\author{Bahman Afsari}

\keyword{Expression Variation Analysis}

\seealso{GSReg.GeneSets.DIRAC,cor}

\examples{


### loading and pruning the pathways
library(GSBenchMark)
data(diracpathways)
### loading the data
data(leukemia_GSEA)


### removing genes which contain not a number.
if(sum(apply(is.nan(exprsdata),1,sum)>0))
  exprsdata = exprsdata[-which(apply(is.nan(exprsdata),1,sum)>0),];
 
### extracting gene names
genenames = rownames(exprsdata);

### DIRAC analysis
VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=exprsdata,
 pathways=diracpathways, phenotypes=as.factor(phenotypes))
E1 = sapply(VarAnKendallV,function(x) x$E1);
E2 = sapply(VarAnKendallV,function(x) x$E2);
Kpvalues = sapply(VarAnKendallV,function(x) x$pvalue);

dysregulatedpathways = rbind(E1[which(Kpvalues<0.05)],
E2[which(Kpvalues<0.05)],Kpvalues[which(Kpvalues<0.05)]);
rownames(dysregulatedpathways)<-c("E1","E2","pvalues");
print(dysregulatedpathways)
plot(x=dysregulatedpathways["E1",],y=dysregulatedpathways["E2",],
xlim=range(dysregulatedpathways[1:2,]),ylim=range(dysregulatedpathways[1:2,]))
lines(x=c(min(dysregulatedpathways[1:2,]),max(dysregulatedpathways[1:2,])),
y=c(min(dysregulatedpathways[1:2,]),max(dysregulatedpathways[1:2,])),type="l")

}
