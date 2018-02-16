############################################
##########        GSReg Package   ##########
##########        Bahman Afsari   ##########
##########        Elana J. Fertig ##########

###########################################
#### EVA Analysis Function
###########################################

GSReg.GeneSets.EVA <- function(geneexpres,
                               pathways,
                               phenotypes,
                               verbose = T,
                               minGeneNum = 5, 
                               distFunc = GSReg.kendall.tau.distance,
                               distparamPathways, 
                               ... )
{
  #Checking the input data
  GSReg.Check.input(prunedpathways=pathways,exprsdata=geneexpres,phenotypes=phenotypes)
  #pruning the pathways
  pathways <- GSReg.Prune(pathways,rownames(geneexpres), minGeneNum)
  values <- vector(mode="list",length(pathways))
  names(values) <- names(pathways)
  samplesC1 <- which(phenotypes==levels(phenotypes)[1])
  samplesC2 <- which(phenotypes==levels(phenotypes)[2])
  if(missing(distparamPathways))
  {
    for(i in seq_along(pathways))
    {
      values[[i]] <- GSReg.Variance(geneexpres[pathways[[i]],], samplesC1, samplesC2,distFunc , ...)
    }
  }else{
    valuesnames = names(values)
    for(i in seq_along(pathways))
    {
      if(verbose == T & (i %% 100) == 0)
        cat(i,"\n")
      values[[i]] <- GSReg.Variance(geneexpres[pathways[[i]],], 
                                    samplesC1, samplesC2,distFunc,
                                    distparamPathways[[valuesnames[i]]], ...)
    }
  }
  return(values)
}

###########################################
###########################################
#### DIRAC Analysis Function
###########################################
###########################################
GSReg.GeneSets.DIRAC <- function(geneexpres,pathways,phenotypes,Nperm=0, alpha = 0, minGeneNum = 5)
{
  GSReg.Check.input(prunedpathways=pathways,exprsdata=geneexpres,phenotypes=phenotypes)
  #prune genes in the pathway
  pathways <- GSReg.Prune(pathways,rownames(geneexpres), minGeneNum)
  #calculate DIRAC variability measure
  mus <- GSReg.DIRAC.Pathways(geneexpres,pathways,phenotypes)
  #Calculating a p-value
  if(Nperm > 0)
  {
    musperm <- vector(mode="list",Nperm)
    for( i in seq_len(Nperm))
    {
      musperm[[i]] <- GSReg.DIRAC.Pathways(geneexpres=geneexpres, pathways=pathways, 
                                           phenotypes=sample(phenotypes),alpha=alpha)
    }
    pvaluesperm <- vector(mode="numeric",length=length(pathways))
    names(pvaluesperm) <- names(pathways) 
    for( i in seq_along(pathways))
    {
      z <- sapply(musperm,function(x) x$diffmu[i])
      pvaluesperm[i] <- mean(abs(mus$mu1[i]-mus$mu2[i])<=abs(z)) 
    }
    return(list(mu1=mus$mu1,mu2=mus$mu2,pvalues=pvaluesperm, zscores = mus$zscores))
  }else{
    return(list(mu1=mus$mu1,mu2=mus$mu2,pvalues=mus$pvalues, zscores = mus$zscores))
  }
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### columns of V
### Outcome: Z(i,j) = #(I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))/(n(n-1)/2)
##################################################

GSReg.kendall.tau.distance <- function(V){
  #Checking if V is a numeric matrix 
  GSReg.Check.input(V)
  
  return(GSReg.kendall.tau.distance.internal(V))
}

##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### Restricted with RestMat
### columns of V
### Outcome: Z(i,j) = #((I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))
##  | Rest(i,j)==1)/(#Rest(i,j)=1)
## if sparse is T matrix
##################################################

GSReg.kendall.tau.distance.restricted <- function(V, RestMat){
  #Checking if V is a numeric matrix 
  
  GSReg.Check.input(V,RestMat = RestMat)
  
  return(GSReg.kendall.tau.distance.Restricted.Sparse.internal(V, 
                                                               as.matrix(RestMat)))
  
}


##################################################
##################################################
### Calculate the normalized kendall-tau-distance between
### columns of V
### Outcome: Z(i,j) = #(I(I(V_l^i<V_k^i)!=I(V_l^j<V_k^j)))/(n(n-1)/2)
##################################################

GSReg.kendall.tau.distance.template <- function(V, Temp){
  #Checking if V is a numeric matrix 
  GSReg.Check.input(V, Temp = Temp)
  
  return(GSReg.kendall.tau.distance.template.internal(V,as.matrix(Temp)))
}


##################################################
##################################################
### Make junction overlap matrices
### genesCoordination a vector of geneRanges with an additional
### of gene coordinates and genesCoordination$SYMBOL has the gene names
### Outcome a list for all genestoStudy
### Each of them includes
### out[["genes"]]["junci","juncj"] = I(junc i and junc j overlap)
##################################################


GSReg.overlapJunction <- function(juncExprs,
                                  GenestoStudy=NULL,
                                  geneexpr=NULL,
                                  minmeanloggeneexp= 3,
                                  alpha =0,
                                  sparse = F,
                                  genesCoordinatesTxDB = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                                  geneIDInTxDB = 'ENTREZID', 
                                  geneIDOut = 'SYMBOL',
                                  org=org.Hs.eg.db, ...){
  
  ##loading genes
  if (!requireNamespace('GenomicRanges', quietly = TRUE)){
    stop("Please install GenomicRanges package for this function.")
  }  
  
  
  gnAll <-  genes(genesCoordinatesTxDB)
  gSymbol <- mapIds(org, keys=as.character(gnAll$gene_id), 
                    keytype = geneIDInTxDB,
                    column = geneIDOut, ...)
  gnAll$SYMBOL <- gSymbol
  
  
  
  if(is.null(GenestoStudy)){
    gn <- gnAll
  }else{
    GenestoStudyId <- which(gnAll$SYMBOL %in% GenestoStudy)
    gn <- gnAll[GenestoStudyId]
  }
  
  
  #calculating mean of the log2(x+1) of the gene expression if they are provided
  if(!is.null(geneexpr)){
    meanloggeneexpr <- apply(log2(geneexpr+1),MARGIN = 1,FUN = mean)
  }
  

  #Filtering the genes not expressed  
  if(!is.null(geneexpr) & minmeanloggeneexp >0 ){
    genenotexpressed <- names(which(meanloggeneexpr<minmeanloggeneexp))
    gn <- gn[which(!(gn$SYMBOL %in% genenotexpressed))]
    
  }
    
    
  #}
  
  ### finding junction locations and put them in the size
  z <- strsplit(sub('-',':',rownames(juncExprs)),':')
  mychr <- sapply(X=z,FUN = function(x) x[1])
  mystart <- sapply(X=z,FUN = function(x) x[2])
  myend <- sapply(X=z,FUN = function(x) x[3])
  
  ### puting junction in GRanges format
  juncRanges <- as.data.frame(strsplit(sub('-',':',rownames(juncExprs)),':')
                              ,stringsAsFactors = FALSE)
  junctionsGRanges <- GRanges(seqnames = Rle(mychr), 
                              ranges = IRanges(start = as.numeric(mystart),
                                               end = as.numeric(myend)))
  
  ### finding junction for each gene
  overlapJunction <- findOverlaps(junctionsGRanges,gn)
  
  #making overlap matrix
  juncnames <- rownames(juncExprs)
  genesJunction <- tapply(rownames(juncExprs)[queryHits(overlapJunction)],
                          gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
  
  
  if(!is.null(geneexpr) & alpha>0){
    
    juncmax <-  apply(log2(juncExprs+1),MARGIN = 1, max)
    junc2genes <- gn$SYMBOL[subjectHits(overlapJunction)]
    names(junc2genes)<-rownames(juncExprs)[queryHits(overlapJunction)]
    

    juncexpressed <- names(which(juncmax[names(junc2genes)]>alpha*meanloggeneexpr[junc2genes]))
  }else{
    juncexpressed <- rownames(juncExprs)
  }
  
  
  genesJunctionInd <- tapply(queryHits(overlapJunction),
                             gn$SYMBOL[subjectHits(overlapJunction)],function(x) x)
  
  if(sparse){
    MyRest <- sapply(X = genesJunctionInd, function(x) { 
      y <- junctionsGRanges[x];
      w <- findOverlaps(y,y);
      mymat <- Matrix(0,nrow = length(y), ncol = length(y),
                      dimnames = list(juncnames[x],juncnames[x]),
                      sparse = T);
      mymat[queryHits(w)+ length(x)*(subjectHits(w)-1)] <- 1
      return(mymat)
    })
    
  }else{
  
    MyRest <- sapply(X = genesJunctionInd, function(x) { 
      y <- junctionsGRanges[x];
      w <- findOverlaps(y,y);
      mymat <- matrix(0,nrow = length(y), ncol = length(y),
                      dimnames = list(juncnames[x],juncnames[x]));
      mymat[queryHits(w)+ length(x)*(subjectHits(w)-1)] <- 1
      return(mymat)
    })
  }
  return(list(Rest = MyRest,
              genesJunction = genesJunction,
              juncexpressed=juncexpressed))
}


##################################################
##################################################
### Splice-EVA (SEVA) algorithm
### There is a filter for junctions if 
### (max junc)*alpha < mean(gene)
##################################################


### You can call this function to filter the junctions that whose log2(expression+1) is 
# less than alpha* log2(its gene expression + 1). The function returns the junction names
GSReg.SEVA <- function(juncExprs, 
                       phenoVect,
                       verbose=T,
                       sparse =F, ...){
  
  
  
  
  
  juncMatrices <- GSReg.overlapJunction(juncExprs,
                                        sparse = F, ...)
  
  
  MyRest <- juncMatrices$Res 
  genesJunction <- juncMatrices$genesJunction
  juncexpressed <- juncMatrices$juncexpressed
  
  if(sparse){
  junctionPValue <- GSReg.GeneSets.EVA(geneexpres = juncExprs[juncexpressed,names(phenoVect)],
                                       phenotypes = phenoVect,
                                       verbose = verbose,
                                       minGeneNum = 2,
                                       pathways = genesJunction[names(which(sapply(genesJunction,length) >2))],
                                       distFunc = GSReg.kendall.tau.distance.Restricted.internal,
                                       distparamPathways = MyRest  )
  }else{
    junctionPValue <- GSReg.GeneSets.EVA(geneexpres = juncExprs[juncexpressed,names(phenoVect)],
                                         phenotypes = phenoVect,
                                         verbose = verbose,
                                         minGeneNum = 2,
                                         pathways = genesJunction[names(which(sapply(genesJunction,length) >2))],
                                         distFunc = GSReg.kendall.tau.distance.Restricted.Sparse.internal,
                                         distparamPathways = MyRest  )
  }
  return(junctionPValue)
  
  
}