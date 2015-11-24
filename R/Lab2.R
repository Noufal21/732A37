library(outliers)
library(car)
www = "http://www.ida.liu.se/~732A37/T1-9.dat"
www2 = "http://www.ida.liu.se/~732A37/T5-12.dat"
www3 = "http://www.ida.liu.se/~732A37/T6-13.dat" 
data <- read.delim(www, header = FALSE, sep="\t")
data2 <- read.delim(www2, header = FALSE, sep="")
data3 <- read.delim(www3, header = FALSE, sep="")
colnames(data) <- c("NAT","100m(s)","200m(s)","400m(s)",
                    "800m(min)","1500m(min)","3000m(min)","Mara(min)")
colnames(data2) <- c("Tail length","Wing length")
colnames(data3) <- c("MaxBreath","BasHeight","BasLength","NasHeight","Time period")
# head(data)
# head(data2)
# head(data3)

nations <- as.character(data[,1])
data <- data[,-1]
col_mu <- colMeans(data)
covmatr <- var(data)

centerdata <- data - matrix(rep(col_mu, dim(data)[1]),ncol=7,byrow=TRUE)
centerdata <- as.matrix(centerdata)
covarianceadjusteddists <-c()

for(i in 1:dim(data)[1]){
  covarianceadjusteddists <- c( covarianceadjusteddists, sqrt(centerdata[i,] %*% solve(covmatr) %*% centerdata[i,]))
}

#We calculate the chi-square test statistic for seven degrees of freedom
#at 0.1% significance level and compare them to the squared Mahalanobis distances from the mean vector
#in order to determine which countries can be considered outliers.
#
for(i in 1:dim(data)[1]){
  if(covarianceadjusteddists[i]^2 > qchisq(0.999,7)){
    print(paste(nations[i]," is an outlier according to chisq test at 0.1% significance level"))
  }
}
#N Korea probably do not perform very extremely, but their records must not exhibit the same
#covariance with each other as that of other nations. One interpretation is that their training
#regimen produced vastly different results from that of most countries. Another interpretation could be
#that the records are fabricated by somebody without knowledge of running record covariances.

plot(data2[,1],data2[,2],main="Wing length vs tail length of female hook-billed kites",
     xlab=names(data2)[1],ylab=names(data2)[2],pch="+",ylim=c(230,320),xlim=c(170,220))
qqnorm(data2[,1],main="Q-Q Plot of tail length")
qqline(data2[,1])
qqnorm(data2[,2],main="Q-Q Plot of wing length")
qqline(data2[,2])
#The scatter-plot as well as the qq-plots suggest that this data is well described by a
#bivariate normal distribution is aviable population model. One can easily imagine an oval shapee
#in the scatter-plot encapsulating all the data points, and the Q-Q plot indicate that the data
#variables are each approximately normally distributed.

alpha <- 0.95
col_mu2 <- colMeans(data2)
vardata2 <-var(data2)
lambda<-eigen(vardata2)
n <- dim(data2)[1]
p <- 2
aa <- c(1,0)
ab <- c(0,1)
quadraa <- aa %*% vardata2 %*% aa
quadrab <- ab %*% vardata2 %*% ab
csquared <- p*(n - 1) /(n * (n- p)) * qf(alpha,p,n-p)
simuldista <- sqrt(csquared * quadraa)
simuldistb <- sqrt(csquared * quadrab)
t_val <- -qt(0.05/(2*p),n-1)
bonfdista <- t_val * sqrt(quadraa/n)
bonfdistb <- t_val * sqrt(quadrab / n)

paste("Simultaneous confidence interval for Tail length: (",round(col_mu2[1] - simuldista,2),
      ",",round(col_mu2[1] + simuldista,2),")")
paste("Simultaneous confidence interval for Wing length: (",round(col_mu2[2] - simuldistb,2),
      ",",round(col_mu2[2] + simuldistb,2),")")
paste("Bonferroni confidence interval for Tail length: (",round(col_mu2[1] - bonfdista,2),
      ",",round(col_mu2[1] + bonfdista,2),")")
paste("Bonferroni confidence interval for Wing length: (",round(col_mu2[2] - bonfdistb,2),
      ",",round(col_mu2[2] + bonfdistb,2),")")


# dataEllipse(data2$`Tail length`,data2$`Wing length`,levels=0.95,
#             ylim=c(230,320),xlim=c(165,220),center.pch = NULL,center.cex = 0.8,
#             main=c("Wing length vs tail length of female hook-billed kites",
#                    "95% mean vector confidence ellipse"),
#             xlab=names(data2)[1],ylab=names(data2)[2],lwd=1)


evs <- sqrt(lambda$values * csquared)
evecs <- lambda$vectors

a <- evs[1]
b <- evs[2]
x0 <- col_mu2[1]
y0 <- col_mu2[2]
alpha <- atan(evecs[ , 1][2] / evecs[ , 1][1])
theta <- seq(0, 2 * pi, length=(1000))

x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)


plot(x, y, type = "l",main=c("Wing length (mm) vs tail length (mm) of female hook-billed kites",
                             "Mean vector 95% confidence ellipse"),
     xlab=names(data2)[1],ylab=names(data2)[2]
     )
points(col_mu2[1],col_mu2[2],cex=2,pch="+")
#according to wikipedia, Most accipitrids exhibit sexual dimorphism in size, 
#, unusually for birds, it is the females that are larger than the males


