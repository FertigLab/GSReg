##################################################
##################################################
### An internal function to check the inputs.
##################################################

GSReg.Check.input <- function(V,
                              pathways,
                              exprsdata,
                              phenotypes,
                              genenames,
                              prunedpathways,
                              RestMat,
                              Temp)
{
  if( !missing(V))
  {
    if(!is.numeric(V) || !is.matrix(V))
      stop("V must be a numeric matrix. Please transform vectors it to column matrices.")
    if(sum(!is.finite(V))>0)
      stop("V must contain finite numbers. NAs,NaN's, and Inf are unacceptable")
    #if( !missing(RestMat))
    #  if(ncol(RestMat)!= nrow(V) || nrow(RestMat)!=nrow(V))
    #    stop("RestMat must be a square matrix with the row nmumbers of V.")
  }
  if(!missing(Temp)){
    if(sum(!(Temp %in% c(0,1)))>0)
      stop("Elements of the template must be either 0 or 1.")
    
    
      if(length(intersect(rownames(Temp),colnames(Temp))) < nrow(Temp)) #Check if the template 
        stop("Template must be sqaure with the same rownames and colnames.")
#       if(sum(Temp + t(Temp)+ diag(nrow(Temp)) 
#              != matrix(data = 1,nrow = nrow(Temp), ncol = ncol(Temp)))>0)
#         stop(" Temp must be an anti-symetric matrix.")
    
    if(!missing(V))
      if(length(intersect(rownames(Temp),rownames(V)))<nrow(V) || nrow(V)!=nrow(Temp))
      stop("Template and V must have the same rownames.")
  }
  if(!missing(pathways))
  {
    if(!is.list(pathways) || sum(!sapply(pathways,is.character)))
      stop("pathways must be list of vectors of characters.")
    if(length(names(pathways))==0)
      stop("pathways require names. If you do not know the name of the pathway assign numbers as the names.")
  }
  if(!missing(exprsdata))
    if(!is.numeric(exprsdata) || ! is.matrix(exprsdata) || is.null(rownames(exprsdata)) )
      stop("exprsdata must be a numeric matrix with gene names as the rownames");
  if(!missing(phenotypes))
    if(!is.factor(phenotypes) || length(levels(phenotypes))!=2)
      stop("The phenotype must be factor consisting of only two levels.")
  if(!missing(exprsdata) && !missing(phenotypes))
    if(length(phenotypes)!=dim(exprsdata)[2])
      stop("The length of phenotypes and the column number of exprsdata must match.")
  if(!missing(genenames))
    if(!is.character(genenames) || ! is.vector(genenames))
      stop("genenames must be a vector of character names.")
  if(!missing(prunedpathways))
  {
    if(!is.list(prunedpathways) || sum(!sapply(prunedpathways,is.character)))
      stop("pathways must be list of vectors of characters.")
    if(length(names(prunedpathways)) == 0)
      stop("pathways require names. If you do not know the name of the pathway, please assign numbers as the names.")
    if(sum(sapply(prunedpathways,length)<2))
      stop("The pruned pathways must contain at least two genes. Please remove pathways with less than two genes.")
  }
  if(!missing(genenames) && !missing(prunedpathways))
  {
    if(length(setdiff(unlist(prunedpathways),genenames)))
      stop("The genes in the pathways must be a subset of the genes in the data expression. Use GSReg.Prune. ")
    sum(sapply(prunedpathways,length)<2)
  }  
}


##################################################
##################################################
### internal function for calculating the variance
### of only one pathway. Pathwayexp contains 
### expression data. samplesC1 and sampleC2 are the
### indices for class 1 and class 2.
##################################################
GSReg.Variance.Wrap <- function(pathwayexp, incVec, samplesC1, samplesC2, distFunc = GSReg.kendall.tau.distance.Wrap,                           
                                ...  ){
  
  distAll <- distFunc(pathwayexp,incVec,...)
  
  
  dist1 <- distAll[samplesC1,samplesC1]
  n1 <- length(samplesC1)
  E1 <- sum(dist1)/(n1*(n1-1))
  Dis1P2 <- sum(dist1^2)
  E1P2 <- Dis1P2/(n1*(n1-1))
  E1DCross <- (sum(apply(dist1,1,sum)^2)-Dis1P2)/(n1*(n1-1)*(n1-2))
  VarD1 <- E1P2 - E1^2
  CovE1Cross <- E1DCross - E1^2
  VarEta1 <- 4*((n1*(n1-1))/2*VarD1 + n1*(n1-1)*(n1-2)*CovE1Cross)/((n1*(n1-2))^2)
  
  dist2 <- distAll[samplesC2,samplesC2]
  n2 <- length(samplesC2)
  E2 <- sum(dist2)/(n2*(n2-1))
  Dis2P2 <- sum(dist2^2)
  E2P2 <- Dis2P2/(n2*(n2-1))
  E2DCross <- (sum(apply(dist2,1,sum)^2)-Dis2P2)/(n2*(n2-1)*(n2-2))
  VarD2 <- E2P2 - E2^2
  CovE2Cross <- E2DCross - E2^2
  VarEta2 <- 4*((n2*(n2-1))/2*VarD2 + n2*(n2-1)*(n2-2)*CovE2Cross)/((n2*(n2-2))^2)
  
  vartotal <- VarEta2 + VarEta1;
  zscore <- (E1-E2)/sqrt(abs(vartotal+0.0000001));
  
  ################ To be checked
  ##### EVA2 
  dist12 <- distAll[samplesC1,samplesC2]
  E12 = sum(dist12)/n1/n2
  #ED1D12 = mean(apply(dist1,1,mean)*apply(dist12,1,mean))
  ED1D12 = sum(apply(dist12,1,sum)*apply(dist1,1,sum))/n1/(n1-1)/n2
  CovE1E12 = ED1D12 - E1*E12
  
  #ED2D12 = mean(apply(dist2,1,mean)*apply(dist12,2,mean))
  ED2D12 = sum(apply(dist12,2,sum)*apply(dist2,1,sum))/n2/(n2-1)/n1
  CovE2E12 = ED2D12 - E2*E12
  
  ##
  #CovD12D12p = mean(apply(dist12,2,mean)^2)- E12^2
  CovD12D12p <- (sum(apply(dist12,1,sum)^2)-sum(dist12^2))/n1/n2/(n2-1) - E12^2
  #CovD12D1p2 = mean(apply(dist12,1,mean)^2)- E12^2
  CovD12D1p2 <- (sum(apply(dist12,2,sum)^2)-sum(dist12^2))/n2/n1/(n1-1) - E12^2
  VarD1D2 <- var(c(dist12))
  
  ### Test
  nAll = n1 + n2
  EAll <- sum(distAll)/(nAll*(nAll-1))
  DisAllP2 <- sum(distAll^2)
  EAllP2 <- DisAllP2/(nAll*(nAll-1))
  EAllDCross <- (sum(apply(distAll,1,sum)^2)-DisAllP2)/(nAll*(nAll-1)*(nAll-2))
  CovEAllCross <- EAllDCross - EAll^2
  ## endTest
  
  ### D12 - D11
  
  D12D1 = (E12-E1)
  lambda = n1/(n1+n2)
  ### Test
  
  ### updated
  #vartotD12D1 = (4/lambda * CovE1Cross - 4/lambda * CovE1E12 
  #            + 1/lambda * CovD12D12p + 1 / (1-lambda) * CovD12D1p2)/(n1+n2)
  
  ### vartotD12D1 = (1/lambda+1/(1-lambda))*CovEAllCross/nAll
  
  vartotD12D1 = CovD12D12p/n1 + CovD12D1p2/n2  + VarEta1 - 4/n1 *CovE1E12# + VarD1D2/n1/n2
  #vartotal <- 4*(1/lambda+1/(1-lambda))*CovEAllCross/nAll;
  zscore <- (E1-E2)/sqrt(abs(vartotal+0.0000001));
  
  ### endTest
  
  zscoreD12D1 = D12D1/sqrt(abs(vartotD12D1)+0.0000001)
  
  lambda = n2/(n1+n2)
  D12D2 = (E12-E2)
  
  ### Test
  #vartotD12D2 = (4/lambda * CovE2Cross - 4/lambda * CovE2E12 
  #              + 1/lambda * CovD12D1p2 + 1 / (1-lambda) * CovD12D12p)/(n1+n2)
  
  #vartotD12D2 = (1/lambda+1/(1-lambda))*CovEAllCross/nAll
  vartotD12D2 = CovD12D12p/n1 + CovD12D1p2/n2  + VarEta2 - 4/n2 *CovE2E12#+ VarD1D2/n1/n2
  ### endTest
  zscoreD12D2 = D12D2/sqrt(abs(vartotD12D2)+0.0000001)
  
  #Test
  myvartotal = 4*(1/lambda+1/(1-lambda))*CovEAllCross/nAll
  
  pvalue=2*(1-pnorm(abs(zscore)))
  pvalueD12D1 = 2*(1-pnorm(abs(zscoreD12D1)))
  pvalueD12D2 = 2*(1-pnorm(abs(zscoreD12D2)))
  
  
  pvalueTotal = apply(cbind(pvalue,pvalueD12D2,pvalueD12D1),
                      MARGIN = 1,
                      FUN = function(x) min(c(min(x)*3,1)))
  
  return(list(E1=E1,E2=E2,E12 = E12,
              zscore=zscore,  zscoreD12D1=zscoreD12D1, zscoreD12D2 = zscoreD12D2,    
              VarEta1=VarEta1,VarEta2=VarEta2, #EVA 1
              sdtotal=sqrt(abs(vartotal)),#EVA 1 sd 
              CovD12D12p = CovD12D12p, CovD12D1p2 = CovD12D1p2, CovD1D12 = CovE1E12, CovD1D12 = CovE1E12,#EVA 2
              vartotD12D1 = vartotD12D1, vartotD12D2 = vartotD12D2,
              pvalue=pvalue,
              pvalueD12D1 = pvalueD12D1, 
              pvalueD12D2 = pvalueD12D2,
              pvalueTotal = pvalueTotal,
              mysdtotal = sqrt(myvartotal),
              distAll=distAll,
              dist1=dist1,
              dist12=dist12))#, myzscore =  (E1-E2)/sqrt(abs(myvartotal+0.0000001)),
  #mypvalue=2*(1-pnorm(abs(myzscore)))))
  
}


##################################################
##################################################
### internal function for calculating the variance
### of only one pathway. Pathwayexp contains 
### expression data. samplesC1 and sampleC2 are the
### indices for class 1 and class 2.
##################################################
GSReg.Variance <- function(pathwayexp, samplesC1, samplesC2,distFunc = GSReg.kendall.tau.distance,                           
                           ...  ){
  
  distAll <- distFunc(pathwayexp,...)
  
  
  dist1 <- distAll[samplesC1,samplesC1]
  n1 <- length(samplesC1)
  E1 <- sum(dist1)/(n1*(n1-1))
  Dis1P2 <- sum(dist1^2)
  E1P2 <- Dis1P2/(n1*(n1-1))
  E1DCross <- (sum(apply(dist1,1,sum)^2)-Dis1P2)/(n1*(n1-1)*(n1-2))
  VarD1 <- E1P2 - E1^2
  CovE1Cross <- E1DCross - E1^2
  VarEta1 <- 4*((n1*(n1-1))/2*VarD1 + n1*(n1-1)*(n1-2)*CovE1Cross)/((n1*(n1-2))^2)
  
  dist2 <- distAll[samplesC2,samplesC2]
  n2 <- length(samplesC2)
  E2 <- sum(dist2)/(n2*(n2-1))
  Dis2P2 <- sum(dist2^2)
  E2P2 <- Dis2P2/(n2*(n2-1))
  E2DCross <- (sum(apply(dist2,1,sum)^2)-Dis2P2)/(n2*(n2-1)*(n2-2))
  VarD2 <- E2P2 - E2^2
  CovE2Cross <- E2DCross - E2^2
  VarEta2 <- 4*((n2*(n2-1))/2*VarD2 + n2*(n2-1)*(n2-2)*CovE2Cross)/((n2*(n2-2))^2)
  
  vartotal <- VarEta2 + VarEta1;
  zscore <- (E1-E2)/sqrt(abs(vartotal+0.0000001));
  
  ################ To be checked
  ##### EVA2 
  dist12 <- distAll[samplesC1,samplesC2]
  E12 = sum(dist12)/n1/n2
  #ED1D12 = mean(apply(dist1,1,mean)*apply(dist12,1,mean))
  ED1D12 = sum(apply(dist12,1,sum)*apply(dist1,1,sum))/n1/(n1-1)/n2
  CovE1E12 = ED1D12 - E1*E12
  
  #ED2D12 = mean(apply(dist2,1,mean)*apply(dist12,2,mean))
  ED2D12 = sum(apply(dist12,2,sum)*apply(dist2,1,sum))/n2/(n2-1)/n1
  CovE2E12 = ED2D12 - E2*E12
  
  ##
  #CovD12D12p = mean(apply(dist12,2,mean)^2)- E12^2
  CovD12D12p <- (sum(apply(dist12,1,sum)^2)-sum(dist12^2))/n1/n2/(n2-1) - E12^2
  #CovD12D1p2 = mean(apply(dist12,1,mean)^2)- E12^2
  CovD12D1p2 <- (sum(apply(dist12,2,sum)^2)-sum(dist12^2))/n2/n1/(n1-1) - E12^2
  VarD1D2 <- var(c(dist12))
  
  ### Test
  nAll = n1 + n2
  EAll <- sum(distAll)/(nAll*(nAll-1))
  DisAllP2 <- sum(distAll^2)
  EAllP2 <- DisAllP2/(nAll*(nAll-1))
  EAllDCross <- (sum(apply(distAll,1,sum)^2)-DisAllP2)/(nAll*(nAll-1)*(nAll-2))
  CovEAllCross <- EAllDCross - EAll^2
  ## endTest
  
  ### D12 - D11
  
  D12D1 = (E12-E1)
  lambda = n1/(n1+n2)
  ### Test
  
  ### updated
  #vartotD12D1 = (4/lambda * CovE1Cross - 4/lambda * CovE1E12 
  #            + 1/lambda * CovD12D12p + 1 / (1-lambda) * CovD12D1p2)/(n1+n2)
  
  ### vartotD12D1 = (1/lambda+1/(1-lambda))*CovEAllCross/nAll
  
  vartotD12D1 = CovD12D12p/n1 + CovD12D1p2/n2  + VarEta1 - 4/n1 *CovE1E12# + VarD1D2/n1/n2
  #vartotal <- 4*(1/lambda+1/(1-lambda))*CovEAllCross/nAll;
  zscore <- (E1-E2)/sqrt(abs(vartotal+0.0000001));
  
  ### endTest
  
  zscoreD12D1 = D12D1/sqrt(abs(vartotD12D1)+0.0000001)
  
  lambda = n2/(n1+n2)
  D12D2 = (E12-E2)
  
  ### Test
  #vartotD12D2 = (4/lambda * CovE2Cross - 4/lambda * CovE2E12 
  #              + 1/lambda * CovD12D1p2 + 1 / (1-lambda) * CovD12D12p)/(n1+n2)
  
  #vartotD12D2 = (1/lambda+1/(1-lambda))*CovEAllCross/nAll
  vartotD12D2 = CovD12D12p/n1 + CovD12D1p2/n2  + VarEta2 - 4/n2 *CovE2E12#+ VarD1D2/n1/n2
  ### endTest
  zscoreD12D2 = D12D2/sqrt(abs(vartotD12D2)+0.0000001)
  
  #Test
  myvartotal = 4*(1/lambda+1/(1-lambda))*CovEAllCross/nAll
  
  pvalue=2*(1-pnorm(abs(zscore)))
  pvalueD12D1 = 2*(1-pnorm(abs(zscoreD12D1)))
  pvalueD12D2 = 2*(1-pnorm(abs(zscoreD12D2)))
  
  
  pvalueTotal = apply(cbind(pvalue,pvalueD12D2,pvalueD12D1),
                      MARGIN = 1,
                      FUN = function(x) min(c(min(x)*3,1)))
  
  return(list(E1=E1,E2=E2,E12 = E12,
          zscore=zscore,  zscoreD12D1=zscoreD12D1, zscoreD12D2 = zscoreD12D2,    
         VarEta1=VarEta1,VarEta2=VarEta2, #EVA 1
         sdtotal=sqrt(abs(vartotal)),#EVA 1 sd 
         CovD12D12p = CovD12D12p, CovD12D1p2 = CovD12D1p2, CovD1D12 = CovE1E12, CovD1D12 = CovE1E12,#EVA 2
         vartotD12D1 = vartotD12D1, vartotD12D2 = vartotD12D2,
         pvalue=pvalue,
         pvalueD12D1 = pvalueD12D1, 
         pvalueD12D2 = pvalueD12D2,
         pvalueTotal = pvalueTotal,
         mysdtotal = sqrt(myvartotal)))#, myzscore =  (E1-E2)/sqrt(abs(myvartotal+0.0000001)),
         #mypvalue=2*(1-pnorm(abs(myzscore)))))
  
}

# ##################################################
# ##################################################
# ### Calculate the normalized DIRAC
# ### columns of V
# ### Outcome: \frac{\sum_{i,j} max(P(X_i<X_j),P(X_i>X_j))}{\sum_i 1}
# ##################################################
# GSReg.DIRAC.mu<-function(V, alpha = 0.0)
# {
#   #Checking if V is a numeric matrix 
#   GSReg.Check.input(V)
#   n <- dim(V)[1]
#   m <- dim(V)[2]
#   # TODO convert for Rcpp
#   d <-.C("Nij",
#         as.double(V),
#         as.integer(dim(V)[1]),
#         as.integer(dim(V)[2]),
#         as.double(matrix(data= 0 ,ncol=n,nrow=n)))
#   pij = d[[4]]/m   #pij = p(V_i<V_j)
# 
#   dim(pij)<-c(n,n)
#   
#   return(sum(pmin(pij,t(pij)))/(n*(n-1)))
# }
# ####################################
# ####################################
# ##### Implements DIRAC analysis for 
# ##### one pathway.
# ####################################
# ####################################
# GSReg.DIRAC.Pathways<-function(geneexpres,pathways,phenotypes){
# 
#   samplesC1 <- which(phenotypes==levels(phenotypes)[1])
#   samplesC2 <- which(phenotypes==levels(phenotypes)[2])
#   
#   
#   mu1 <- sapply(pathways,function(x) GSReg.DIRAC.mu(geneexpres[x,samplesC1]))
#   mu2 <- sapply(pathways,function(x) GSReg.DIRAC.mu(geneexpres[x,samplesC2]))
#   
#   names(mu1) <- names(pathways)
#   names(mu2) <- names(pathways)
#   
#   return(list(mu1=mu1,mu2=mu2,diffmu=mu1-mu2))
# }

##################################################
##################################################
### Calculate the normalized DIRAC
### columns of V
### Outcome: \frac{\sum_{i,j} max(P(X_i<X_j),P(X_i>X_j))}{\sum_i 1}
##################################################
GSReg.DIRAC.mu<-function(V, alpha = 0.0)
{
  #Checking if V is a numeric matrix 
  GSReg.Check.input(V)
  n <- dim(V)[1]
  m <- dim(V)[2]
  d <- Nij(V)
  pij <- d/m
  
  dim(pij) <- c(n,n)
  rownames(pij) <- rownames(V)
  colnames(pij) <- rownames(V)
  
  mu <- (pij> 0.5 + alpha/m) #Template
  mu <- mu - diag(diag(mu))
  distances <- GSReg.kendall.tau.distance.template(V = V, mu)
  meandists <-  mean(distances)
  vardists <- var(distances)
  return(list( mu = mu, mean = meandists, var = vardists))
}
####################################
####################################
##### Implements DIRAC analysis for 
##### one pathway.
####################################
####################################
GSReg.DIRAC.Pathways<-function(geneexpres,pathways,phenotypes,alpha = 0.0){
  
  samplesC1 <- which(phenotypes==levels(phenotypes)[1])
  samplesC2 <- which(phenotypes==levels(phenotypes)[2])
  
  n1 <- length(samplesC1)
  n2 <- length(samplesC2)
  
  
  inf1 <- lapply(pathways,function(x) GSReg.DIRAC.mu(geneexpres[x,samplesC1],alpha=alpha/n1))
  inf2 <- lapply(pathways,function(x) GSReg.DIRAC.mu(geneexpres[x,samplesC2],alpha=alpha/n2))
  
  names(inf1) <- names(pathways)
  names(inf2) <- names(pathways)
  
  temp1 <- sapply(X = inf1,FUN = function(x) x$mu)
  temp2 <- sapply(X = inf2,FUN = function(x) x$mu)
  
  mu1 <- sapply(X = inf1,FUN = function(x) x$mean)
  mu2 <- sapply(X = inf2,FUN = function(x) x$mean)
  
  var1 <- sapply(X = inf1,FUN = function(x) x$var)
  var2 <- sapply(X = inf2,FUN = function(x) x$var)
  
  zscores <- (mu1-mu2)/sqrt(var1/n1+var2/n2+1e-7/n1/n2)
  pvalues <- 2*(1-pnorm(abs(zscores)))
  
  return(list(mu1=mu1,mu2=mu2,diffmu=mu1-mu2, zscores = zscores, pvalues=pvalues ))
}

##################################################
##################################################
## GSReg.Prune function prunes the pathway list 
## based on the genenames 
####################################################
GSReg.Prune <- function(pathways,genenames, minGeneNum=5)
{
  prunedpathways <- lapply( pathways,function(x) intersect(x,genenames))
  emptypathways <- which(sapply(prunedpathways,length)<minGeneNum)
  if(length(emptypathways)>0)
  {
    message(paste("After pruning the pathways, there exist pathways with zero or one gene!\n Small pathways were deleted. Deleted pathways:","\n",paste(names(emptypathways),collapse='\n')))
    return(prunedpathways[-emptypathways])  
  }else{ return(prunedpathways)}
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### columns of V
### Outcome: Z(i,j) = #(I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))/(n(n-1)/2)
##################################################

GSReg.kendall.tau.distance.internal.Wrap <- function(V,incVec){
  #Checking if V is a numeric matrix 
  #GSReg.Check.input(V)
  
  n <- dim(V)[1]
  m <- dim(V)[2]
  dist <- kendalltaudistWrap(as.vector(V),as.integer(n),as.integer(m),as.vector(incVec))
  #d <- .C("kendalltaudistWrap",
  #        as.double(V),
  #        as.integer(dim(V)[1]),
  #        as.integer(dim(V)[2]),
  #        as.double(incVec),
  #        as.double(matrix(data= 0 ,ncol=m,nrow=m)))
  #dist <- d[[5]]
  dim(dist) <- c(m,m)
  #for (i in 1:m){
  #  for (j in 1:m){
  #    tempI=sum(incVec[,i]&incVec[,j])
  #    tempDenom=tempI*(tempI-1)/2
  #    dist[i,j]=dist[i,j]/tempDenom
  #  }
  #}
  rownames(dist) <- colnames(V)
  colnames(dist) <- colnames(V)
  return(dist)
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### columns of V
### Outcome: Z(i,j) = #(I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))/(n(n-1)/2)
##################################################

GSReg.kendall.tau.distance.internal <- function(V){
  #Checking if V is a numeric matrix 
  #GSReg.Check.input(V)
  
  n <- dim(V)[1]
  m <- dim(V)[2]
  dist <- kendalltaudist(V)
  dim(dist) <- c(m,m)
  rownames(dist) <- colnames(V)
  colnames(dist) <- colnames(V)
  return(2*dist/(n*(n-1)))
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### Restricted with RestMat
### columns of V
### Outcome: Z(i,j) = #((I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))
##  | Rest(i,j)==1)/(#Rest(i,j)=1)
##################################################

GSReg.kendall.tau.distance.Restricted.internal <- function(V, RestMat){
  #Checking if V is a numeric matrix 
  
  #GSReg.Check.input(V)
  
  inter = intersect(rownames(V),rownames(RestMat))
  if( length(inter)< 2)
    stop("No common restrictions.")
  V = V[inter,]
  RestMat = RestMat[inter,inter]
  
  n <- dim(V)[1]
  m <- dim(V)[2]
  
  dist <- kendalltaudistRestricted(V, RestMat)
  dim(dist) <- c(m,m)
  rownames(dist) <- colnames(V)
  colnames(dist) <- colnames(V)
  return(dist/(sum(RestMat)+0.000001))
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### columns of V
### Outcome: Z(i,j) = #(I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))/(n(n-1)/2)
##################################################



GSReg.kendall.tau.distance.Restricted.Sparse.internal <- function(V, RestMat){
  #Checking if V is a numeric matrix 
  
  #GSReg:::GSReg.Check.input(V)
  
  inter = intersect(rownames(V),rownames(RestMat))
  if( length(inter)< 2)
    stop("No common restrictions.")
  V = V[inter,]
  RestMat <- as(RestMat[inter,inter],"matrix")
  
  n <- dim(V)[1]
  m <- dim(V)[2]

  dist <- kendalltaudistRestricted(V, RestMat)
  dim(dist) <- c(m,m)
  rownames(dist) <- colnames(V)
  colnames(dist) <- colnames(V)
  return(dist/(sum(RestMat)+0.000001))
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### columns of V
### Outcome: Z(i,j) = #(I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))/(n(n-1)/2)
##################################################

GSReg.kendall.tau.distance.template.internal <- function(V, Temp){
  #Checking if V is a numeric matrix 
  #GSReg.Check.input(V, Temp = Temp)
  
  n <- dim(V)[1]
  m <- dim(V)[2]

  matchedTemp <- Temp[rownames(V),rownames(V)]
  mydist <- kendalltaudistFromTemp(V, matchedTemp)

  names(mydist) <- colnames(V)
  return(mydist/(sum(Temp)-sum(diag(Temp))+1e-7))
}