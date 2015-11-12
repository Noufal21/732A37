library(rgl)
library(scatterplot3d)
library(ggplot2)
www = "http://www.ida.liu.se/~732A37/T1-9.dat"
data <- read.delim(www, header = FALSE, sep="\t")
nations <- as.character(data[,1])
data <- data[,-1]
colnamez <- c("100m(s)","200m(s)","400m(s)",
                    "800m(min)","1500m(min)","3000m(min)","Mara(min)")
colnames(data) <- colnamez
col_mu <- round(colMeans(data),3)

col_sigmasquared<- apply(data,2,var)

mins <- apply(data,2,min)
maxs <- apply(data,2,max)
print('Means')
col_mu
paste("standard deviations:" )
round(sqrt(col_sigmasquared),3)
paste("Min values:" )
mins
paste("Max values:" )
maxs
par(mfrow = c(2,4))
for(i in 1:3){
boxplot(data[,i],ylab="seconds",main=c("Records for distance:",
                                           colnames(data)[i]))
}
for(j in 4:7){
boxplot(data[,j],ylab="minutes",main=c("Records for distance:",
                                           colnames(data)[j]))
}
plot.new()

par(mfrow = c(2,4))
for(i in 1:length(data)){
  plot(seq_along(data[,i]),data[,i],main=c("Records for distance:",
                                           colnames(data)[i]),xlab="",ylab="time unit")
}

plot.new()
for(i in 1:length(data)){
  qqnorm(data[,i],main=c("Q-Q plot of records",
                                           colnames(data)[i]))
  qqline(data[,i])
}
par(mfrow = c(1,1))
covmatr <- var(data)
corrmatr <- cor(data)

print("Covariance matrix:")
signif(covmatr,3)
print("Correlation matrix:")
signif(corrmatr,3)
pairs(data)
## for(i in 1:(length(data)-2)){
##   plot3d(data[,i],data[,i+1],data[,i+2],
##                 xlab="even smaller distance record",ylab="smaller distance record",
##          zlab="record",col="red",size=3,type="s")
## }
## 
par(mfrow = c(1,1))
plot(data[,3],data[,4],pch="", main="Country comparison",xlab=colnames(data)[3],ylab=colnames(data)[4])
text(data[,3],data[,4],labels=nations)
plot(data[,1],data[,7],pch="", main="Country comparison",xlab=colnames(data)[1],ylab=colnames(data)[7])
text(data[,1],data[,7],labels=nations)
plot(data[,5],data[,6],pch="", main="Country comparison",xlab=colnames(data)[5],ylab=colnames(data)[6])
text(data[,5],data[,6],labels=nations)
centerdata <- data - matrix(rep(col_mu, dim(data)[1]),ncol=7,byrow=TRUE)
centerdata <- as.matrix(centerdata)
euclideandist <- c()
varianceadjusteddists <- c()
V <- diag(col_sigmasquared)
covarianceadjusteddists <-c()
for(i in 1:dim(data)[1]){
  euclideandist <- c( euclideandist, sqrt(sum(centerdata[i,] * centerdata[i,])))
  varianceadjusteddists <- c( varianceadjusteddists,
            sqrt(centerdata[i,] %*% solve(V) %*% centerdata[i,]))
  covarianceadjusteddists <- c( covarianceadjusteddists,
            sqrt(centerdata[i,] %*% solve(covmatr) %*% centerdata[i,]))
}

plot(euclideandist,pch="", main="Country comparison",ylab="Euclidean distance from center")
text(euclideandist,labels=nations)
plot(varianceadjusteddists,pch="",
     main="Country comparison",ylab="Variance adjusted distance from center")
text(varianceadjusteddists,labels=nations)
plot(covarianceadjusteddists,pch="",
     main="Country comparison",ylab="Covariance adjusted distance from center")
text(covarianceadjusteddists,labels=nations)
## 
## 
