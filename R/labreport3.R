www = "http://www.ida.liu.se/~732A37/T1-9.dat"
data <- read.delim(www, header = FALSE, sep="\t")
colnames(data) <- c("NAT","100m(s)","200m(s)","400m(s)",
                    "800m(min)","1500m(min)","3000m(min)","Mara(min)")
#head(data)
nations <- as.character(data[,1])
data <- data[,-1]
col_mu <- colMeans(data)
col_sigmasquared<-apply(data,2,var)

corrmatr <- cor(data)
corrmatr
eigenstuff <-eigen(corrmatr)
eigenstuff
data2 <- scale(data)
std_lambdas <- eigen(cor(data2))$values #Same as eigenstuff$values
pcastuff <- prcomp(data2,center = FALSE,scale.=FALSE)

pcastuff$rotation[,1:2]

pccorrstdvars <- matrix(0,nrow=2,ncol=7)

for(i in 1:2){
  for(j in 1:7){
    pccorrstdvars[i,j] <- pcastuff$rotation[j,i] * sqrt(std_lambdas[i]) 
  }
}
paste("Correlations between the two first PCs and the standardized data variables:")
rownames(pccorrstdvars) <- c("PC1","PC2")
colnames(pccorrstdvars) <- c("100m(s)","200m(s)","400m(s)",
                       "800m(min)","1500m(min)","3000m(min)","Mara(min)")
pccorrstdvars
paste("Percent of total sample variance explained by first PC:",
      sprintf("%2.3f",eigenstuff$values[1]/sum(eigenstuff$values)*100))
paste("Percent of total sample variance explained by second PC:",
      sprintf("%2.3f",eigenstuff$values[2]/sum(eigenstuff$values)*100))
paste("Cumulative Percent of total sample variance explained by first two PCs:",
      sprintf("%2.3f",(eigenstuff$values[2]+eigenstuff$values[1])
              /sum(eigenstuff$values)*100))
pcscores <- data.frame(PC1=pcastuff$x[,1], PC2=pcastuff$x[,2])
plot(pcscores,pch="",main=c("Score plot for PC1 and PC2",
     "of national track records for women"))
text(pcscores[,1],pcscores[,2],labels=nations,cex=0.7)
MLfactanal <- factanal(data,factors = 2,scores = "Bartlett")
MLfactanal
plot(MLfactanal$scores, pch="",main="Factor scores for ML factor analysis")
text(MLfactanal$scores, labels=nations,cex=0.7)
par(mfrow=c(2,2))
for(i in 1:7){
  qqnorm(data[,i])
  qqline(data[,i])
}
par(mfrow=c(1,1))
eigenthings <- eigen(cov(data))
estsqrteigenvalS <- sqrt(eigenthings$values)
paste("percentage of sample variance explained by first common factor:",
      signif(eigenthings$values[1] / sum(diag(cov(data))),3))
L <- estsqrteigenvalS[1] * eigenthings$vectors[,1]
L <- as.matrix(L,drop=FALSE)
paste("loadings:")
L
loadings <- L %*% t(L)
communalities <- diag(loadings)
paste("communalities:")
communalities
specificfactors <- diag(x=diag(cov(data) - loadings))
paste("Uniquenesses:")
diag(specificfactors)
residualmatr <- cov(data) - loadings - specificfactors
paste("Residual Matrix:")
residualmatr

centereddata <- t(as.matrix(data - col_mu))
scores <- as.vector(solve(t(L) %*% L) %*% t(L) %*% centereddata)
scores <- as.data.frame(cbind(scores,rep(0,length(scores))))
plot(scores,pch="",ylim=c(-0.1,0.1),main="Factor 1 score plot for sample covariance matrix S",
     ylab="",xlab="Factor 1 score")
text(scores,labels=nations,cex=0.7)
estsqrteigenvalS2 <- sqrt(eigenstuff$values)
paste("percentage of sample correlation explained by first common factor:",
      signif(eigenstuff$values[1] / 7 * 100,3))
paste("percentage of sample correlation explained by second common factor:",
      signif(eigenstuff$values[2] / 7 * 100,3))
L2 <- estsqrteigenvalS2[1:2] * eigenstuff$vectors[,1:2]
L2 <- as.matrix(L2,drop=FALSE)
paste("loadings:")
L2
loadings2 <- L2 %*% t(L2)
specificfactors2 <- diag(x=diag(corrmatr - loadings2))
paste("Uniquenesses:")
diag(specificfactors2)
residualmatr2 <- corrmatr - loadings2 - specificfactors2
paste("Residual Matrix:")
residualmatr2
scores2 <- (solve(t(L2) %*% L2) %*% t(L2) %*% t(data2))
plot(scores2[1,],scores2[2,], pch = "",main="Factor score plot for sample correlation matrix R",
     xlab="Factor 1",ylab="Factor 2")
text(scores2[1,],scores2[2,],labels=nations,cex=0.7)
## NA
