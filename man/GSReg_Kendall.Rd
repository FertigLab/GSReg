\name{GSReg.kendall.tau.distance}

\alias{GSReg.kendall.tau.distance}
\alias{GSReg.kendall.tau.template}
\alias{GSReg.kendall.tau.distance.restricted}



\title{Calculate Kendall-tau distance in different forms}

\description{Different types of calculations of kendall-tau distance.}

\usage{
GSReg.kendall.tau.distance(V)
GSReg.kendall.tau.distance.template(V, Temp)
GSReg.kendall.tau.distance.restricted(V, RestMat)
}


\arguments{
\item{V}{A matrix on which the distances will be calculated. The distance will be calculated between any pair of  columns}
\item{Temp}{A square binary template matrix, like DIRAC template, whose columns and rows correspond to the rows of V. In general, Temp provides a template for comparisons in V. Temp[i,j]==1 means that the template expects that V[i,.]<V[j,.] and 0 means the reverse.}
\item{RestMat}{Restricted matrix. It must be a square matrix and symmetric with binary 0 or 1 whose columns and rows correspond to the rows of V. Only comparisons with one will be considered.}

}

\details{
GSReg.kendall.tau.distance returns kendall-tau calculates distance matrix between any pair of columns of V.


GSReg.kendall.tau.distance.template returns kendall-tau calculates distance matrix between any column V and a template. Temp[i,j] represent a comaprison between the i-th and j-th element of a vector. Hence, the out come's k-th element is (V[i,k]<V[j,k]  & RestMat[i,j]==1)/((nrow(RestMat)*(nrow(RestMat)-1))/2).


GSReg.kendall.tau.distance.restricted calculates the Kendall-tau distance and the only considered comparisons are those RestMat[i,j]==1. It is a modified Kendall-tau distance used by SEVA.


}



\value{

Kendall-tau distance.
}
\author{Bahman Afsari}

\keyword{Junction overlap}

\seealso{GSReg.GeneSets.DIRAC,GSReg.GeneSets.EVA}

\examples{


library(GSReg)
V <- cbind(c(1,5,3),c(3,2,1))
rownames(V) <- c("F1","F2","F3")
colnames(V) <- c("S1","S2")


myRest1 <- cbind(c(0,1,1),c(1,0,1),c(1,1,0))
rownames(myRest1) <- rownames(V) 
colnames(myRest1) <- rownames(V) 

GSReg.kendall.tau.distance.restricted(V,myRest1)

GSReg.kendall.tau.distance(V)

myRest2 <- cbind(c(0,0,1),c(0,0,1),c(1,1,0))
rownames(myRest2) <- rownames(V) 
colnames(myRest2) <- rownames(V) 
GSReg.kendall.tau.distance.restricted(V,myRest2)


Temp1 <- cbind(c(0,1,1),c(0,0,0),c(0,1,0))
rownames(Temp1) <- rownames(V)
colnames(Temp1) <- rownames(V)

GSReg.kendall.tau.distance.template(V,Temp = Temp1)
}


