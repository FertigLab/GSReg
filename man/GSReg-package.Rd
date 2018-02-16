\name{GSReg-package}
\alias{GSReg-package}
\alias{GSReg}
\docType{package}


\title{A package for Gene Set Analysis based on the variability of gene expression in different phenotypes.}

\description{
  The GSReg package applies the analysis of variety among phenotypes for each gene set from gene expression as well as alternative splicing from junction expression. The user can also use Differential Rank Conservation (DIRAC) (Eddy et al. 2010) 
  and a modified version which allows for efficient and easy p value calculation. Both DIRAC and its modified version are rank-based methods, i.e. they only consider the ordering of the expressions within the pathway.  
}


\section{GSReg package features}{
  The package contains several utilities enabling to:
  \describe{
   {A)  Prune Gene Sets based on the available genes in the expression data;}

  {B)  Calculate the DIRAC measure and p-value for it based on permutation test;}
  
  {C)  Calculate for a modified DIRAC method and a fast-efficient p-valuebased on U-Statistic theory;}
  
  {D)  Alternative Splicing Events in genes from a phenotype to other phenotype using their junction expression;}

  
  
}}


\author{
  Bahman Afsari \email{bahman.afsari@gmail.com},
  Elana J. Fertig \email{ejfertig@jhmi.edu}
  }


\source{
  \url{http://www.ncbi.nlm.nih.gov/pubmed/20523739}

}

\references{

	Eddy et al., 
	 "Identifying tightly regulated and variably expressed networks by Differential Rank Conservation (DIRAC).",	
 \emph{PLoS Comp. Bio.}, 2010, \bold{6}(5)
}


\keyword{ package }
