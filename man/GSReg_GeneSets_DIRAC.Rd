\name{GSReg.GeneSets.DIRAC}

\alias{GSReg.GeneSets.DIRAC}

\title{Performs DIRAC for gene set analysis from the paper Eddy et al (2010).}

\description{GSReg.GeneSets.DIRAC performs DIRAC for gene set analysis from the paper Eddy et al (2010). In fact, the Null hypothesis is that the conservation index is not significantly different under two phenotypes.
The function calculates the p-value using permutation test; hence, extremely low p-value cannot be reached.}

\usage{
GSReg.GeneSets.DIRAC(geneexpres, pathways, phenotypes, Nperm = 0, alpha = 0.05, minGeneNum=5) 
}

\arguments{
\item{geneexpres}{the matrix of gene expressions. The rownames must represent gene names and the columns represent samples. There must not be any missing values. Please use imputation or remove the genes with missing values.}
\item{pathways}{a list containing pathway information. Each element represents a pathway as a character vector. The genes shown in the pathway must be present in geneexpres. Please prune the genes in pathways using GSReg.Prune before applying this function.}
\item{phenotypes}{a binary factor containing the phenotypes for samples in geneexpres; hence, the column number of geneexpres and the length of phenotypes must be equal. }
\item{Nperm}{The number of permutation tests required for p-value calculation. If Nperm==0 then the function reports a normal approximation of the p-value.}
\item{alpha}{A parameter smoothes the template estimate. The template corresponding to a comparison is 1, if the probability of the comparison of the two genes is bigger than 0.5 +alpha/n (where n is the number of samples in the phenotype), and otherwise zero.}
\item{minGeneNum}{the minimum number of genes required in a pathway.}
}

\value{
IThe output is a list with three elements. Each element of the output list is a vector are named according to the pathway.
\item{$mu1}{a vector containing the variability in DIRAC sense (1- conservation indices in Eddy et al (2010) paper) for all pathways in phenotype == levels(phenotypes)[1].}
\item{$mu2}{a vector containing the variability in DIRAC sense (1- conservation indices in Eddy et al (2010) paper)  for all pathways in phenotype == levels(phenotypes)[2].}
\item{$pvalues}{a vector containing p-values for each pathway. Low p-values means that the gene expressions have different orderings under different phenotypes. }
}

\author{Bahman Afsari}

\keyword{DIRAC Analysis}

\seealso{GSReg.GeneSet.VReg}

\examples{
library(GSBenchMark)
### loading and pruning the pathways
data(diracpathways)
### loading the data
data(leukemia_GSEA)


### extracting gene names
genenames = rownames(exprsdata);


### DIRAC analysis
DIRAna = GSReg.GeneSets.DIRAC(pathways=diracpathways,geneexpres=exprsdata,Nperm=0,phenotypes=phenotypes)
dysregulatedpathways = rbind(DIRAna$mu1[which(DIRAna$pvalues<0.05)],
DIRAna$mu2[which(DIRAna$pvalues<0.05)],DIRAna$pvalues[which(DIRAna$pvalues<0.05)]);
rownames(dysregulatedpathways)<-c("mu1","mu2","pvalues");
print(dysregulatedpathways[,1:5])
plot(x=dysregulatedpathways["mu1",],y=dysregulatedpathways["mu2",],
xlim=range(dysregulatedpathways[1:2,]),ylim=range(dysregulatedpathways[1:2,]))
lines(x=c(min(dysregulatedpathways[1:2,]),max(dysregulatedpathways[1:2,])),
y=c(min(dysregulatedpathways[1:2,]),max(dysregulatedpathways[1:2,])),type="l")

}
\references{
Eddy, James A., et al. "Identifying tightly regulated and variably expressed networks by Differential Rank Conservation (DIRAC)." 
PLoS computational biology 6.5 (2010): e1000792. }