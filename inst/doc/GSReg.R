### R code from vignette source 'GSReg.Rnw'

###################################################
### code chunk number 1: start
###################################################
options(width=85)
options(continue=" ")
#rm(list=ls())


###################################################
### code chunk number 2: GSReg.Rnw:165-170
###################################################
library(GSBenchMark)
data(diracpathways)
class(diracpathways)
names(diracpathways)[1:5]
class(diracpathways[[1]])


###################################################
### code chunk number 3: GSReg.Rnw:179-181
###################################################
data(GSBenchMarkDatasets)
print(GSBenchMark.Dataset.names)


###################################################
### code chunk number 4: GSReg.Rnw:187-190
###################################################
DataSetStudy = GSBenchMark.Dataset.names[[9]]
print(DataSetStudy)
data(list=DataSetStudy)


###################################################
### code chunk number 5: GSReg.Rnw:199-201
###################################################
if(sum(apply(is.nan(exprsdata),1,sum)>0))
  exprsdata = exprsdata[-which(apply(is.nan(exprsdata),1,sum)>0),];


###################################################
### code chunk number 6: GSReg.Rnw:205-207
###################################################
genenames = rownames(exprsdata);
genenames[1:10]


###################################################
### code chunk number 7: GSReg.Rnw:237-238
###################################################
library(GSReg)


###################################################
### code chunk number 8: GSReg.Rnw:244-248
###################################################
Nperm = 10
system.time({DIRACperm =GSReg.GeneSets.DIRAC(exprsdata,diracpathways,phenotypes,Nperm=Nperm)})
system.time({DIRACAn =GSReg.GeneSets.DIRAC(exprsdata,diracpathways,phenotypes)})



###################################################
### code chunk number 9: GSReg.Rnw:251-252
###################################################
hist(DIRACAn$pvalues,xlab="pvalue",main="Hist of pvalues applying DIRAC Analysis.")


###################################################
### code chunk number 10: GSReg.Rnw:259-268
###################################################

plot(x=abs(DIRACAn$zscores),y=DIRACperm$pvalues,xlab="|Z-score|",
     ylab="p-value",col="red1",main="DIRAC p-value comparisons")
zscorelin <- seq(min(abs(DIRACAn$zscores)),max(abs(DIRACAn$zscores)),by = 0.1)
pvaltheoretic = (1-pnorm(zscorelin))*2
lines(x=zscorelin,y=pvaltheoretic,type="l",pch=50,lty=5,col="darkblue")
legend("topright",legend=c("permutation test","Normal Approx."),
       col=c("red1","blue"),text.col=c("red1","blue"),
       lty=c(NA,1),lwd=c(NA,2.5),pch=c(21,NA))


###################################################
### code chunk number 11: GSReg.Rnw:300-309
###################################################
#Calculating the variance for the pathways
#Calculate how much it takes to calculate the statistics and their p-value for all pathways

system.time({VarAnKendallV = GSReg.GeneSets.EVA(geneexpres=exprsdata,
  pathways=diracpathways, phenotypes=as.factor(phenotypes)) })


names(VarAnKendallV)[[1]]
VarAnKendallV[[1]]


###################################################
### code chunk number 12: GSReg.Rnw:320-337
###################################################
Nperm = 10;
VarAnPerm = vector(mode="list",length=Nperm)
for( i in seq_len(Nperm))
{
  VarAnPerm[[i]] = GSReg.GeneSets.EVA(geneexpres=exprsdata, pathways=diracpathways, 
                                           phenotypes=sample(phenotypes))
}

pvaluesperm = vector(mode="numeric",length=length(VarAnPerm[[1]]))

for( i in seq_along(VarAnPerm[[1]]))
{
  z = sapply(VarAnPerm,function(x) x[[i]]$E1 - x[[i]]$E2)
  pvaluesperm[i] = mean(abs(VarAnKendallV[[i]]$E1-VarAnKendallV[[i]]$E2)<abs(z)) 
}
 zscore = sapply(VarAnKendallV,function(x) x$zscore);
 pvalustat = sapply(VarAnKendallV,function(x) x$pvalue);


###################################################
### code chunk number 13: GSReg.Rnw:342-351
###################################################
 plot(x=abs(zscore),y=pvaluesperm,xlab="|Z-score|",
  ylab="p-value",col="red1",main="p-value comparisons")
 zscorelin = seq(0,6,0.1);
 pvaltheoretic = (1-pnorm(zscorelin))*2
 lines(x=zscorelin,y=pvaltheoretic,type="l",pch=50,lty=5,col="darkblue")

  legend("topright",legend=c("permutation test","U-Stat Estimation"),
         col=c("red","blue"),text.col=c("red","blue"),
  lty=c(NA,1),lwd=c(NA,2.5),pch=c(21,NA))


###################################################
### code chunk number 14: GSReg.Rnw:358-359
###################################################
hist(x=pvalustat,breaks=20,main="P-value Hist of U-Stat",xlim=c(0,1))


###################################################
### code chunk number 15: GSReg.Rnw:372-378
###################################################
plot(x=DIRACAn$pvalues,y=pvalustat,xlab ="DIRAC",
        ylab="EVA",main=sprintf("P-value Comparison corr=%2.2g",cor(x=DIRACAn$pvalues,y=pvalustat)))

lmfit = lm(pvalustat~DIRACAn$pvalues-1)
abline(lmfit)
cor.test(x=DIRACAn$pvalues,y=pvalustat)


###################################################
### code chunk number 16: GSReg.Rnw:382-383
###################################################
cor(x=DIRACAn$pvalues,y=pvalustat)


###################################################
### code chunk number 17: GSReg.Rnw:389-397
###################################################
load("DIRACAn.rda")
significantPathwaysDIRAC = names(DIRACAn$mu1)[which(DIRACAn$pvalues<0.05)];
print(significantPathwaysDIRAC)
mu1 = DIRACAn$mu1[significantPathwaysDIRAC];
mu2 = DIRACAn$mu2[significantPathwaysDIRAC];
significantPathwaysGSV = names(which(pvalustat<0.05));
eta1 = sapply(VarAnKendallV,function(x) x$E1)[significantPathwaysGSV];
eta2 = sapply(VarAnKendallV,function(x) x$E2)[significantPathwaysGSV];


###################################################
### code chunk number 18: GSReg.Rnw:405-409
###################################################
  plot(x=DIRACAn$pvalues,y=pvalustat,xlab ="DIRAC",
        ylab="EVA",main=sprintf("P-value Comparison corr=%2.2g",cor(x=DIRACAn$pvalues,y=pvalustat)))
  lmfit = lm(pvalustat~DIRACAn$pvalues-1)
abline(lmfit)


###################################################
### code chunk number 19: GSReg.Rnw:425-426 (eval = FALSE)
###################################################
## DIRACAn =GSReg.GeneSets.DIRAC(exprsdata,diracpathways,phenotypes,Nperm=1000)


###################################################
### code chunk number 20: GSReg.Rnw:429-439
###################################################
significantPathwaysDIRAC = names(DIRACAn$mu1)[which(DIRACAn$pvalues<0.05)];
mu1 = DIRACAn$mu1[significantPathwaysDIRAC];
mu2 = DIRACAn$mu2[significantPathwaysDIRAC];

#The dysregulated pathways
names(mu1)
plot(x=mu1,y=mu2,
     xlim=c(0,max(mu1,mu2)),ylim=c(0,max(mu1,mu2)),xlab="normal",ylab="disease",
     main="(a)  DIRAC significantly dysregulated pathways")
lines(x=c(0,max(mu1,mu2)),y=c(0,max(mu1,mu2)))


###################################################
### code chunk number 21: GSReg.Rnw:454-463
###################################################
significantPathwaysGSV = names(which(pvalustat<0.05));
eta1 = sapply(VarAnKendallV,function(x) x$E1)[significantPathwaysGSV];
eta2 = sapply(VarAnKendallV,function(x) x$E2)[significantPathwaysGSV];

#The dysregulated pathways
names(eta1)
plot(x=eta1,y=eta2,xlim=c(0,max(eta1,eta2)),ylim=c(0,max(eta1,eta2)),xlab="normal",ylab="disease",
     main="(b) EVA: Dysregulated pathways")
lines(x=c(0,max(eta1,eta2)),y=c(0,max(eta1,eta2)))


###################################################
### code chunk number 22: GSReg.Rnw:469-471
###################################################
print(significantPathwaysGSV)
print(significantPathwaysDIRAC)


###################################################
### code chunk number 23: GSReg.Rnw:547-558
###################################################
require('Homo.sapiens')
require('org.Hs.eg.db')
require('GenomicRanges')

data(juncExprsSimulated)


overlapMat <- GSReg.overlapJunction(juncExprs =  junc.RPM.Simulated,
                                    geneexpr = geneExrsGSReg)
print(overlapMat[["Rest"]][["CMPK1"]])



###################################################
### code chunk number 24: GSReg.Rnw:562-586
###################################################
V <- cbind(c(1,5,3),c(3,2,1))
rownames(V) <- c("F1","F2","F3")
colnames(V) <- c("S1","S2")
GSReg.kendall.tau.distance(V)

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


###################################################
### code chunk number 25: GSReg.Rnw:589-599
###################################################

data(juncExprsSimulated)
SEVAjunc <- GSReg.SEVA(juncExprs = junc.RPM.Simulated,
                       phenoVect = phenotypes,
                       geneexpr = geneExrsGSReg)
                       
print(sapply(SEVAjunc,function(x) x$pvalue))
#if you want to check Translational as well you can use 2 other p-values
print(sapply(SEVAjunc,function(x) x$pvalueTotal))



###################################################
### code chunk number 26: sessioInfo
###################################################
toLatex(sessionInfo())


