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
  covarianceadjusteddists <- c( covarianceadjusteddists,
                                sqrt(centerdata[i,] %*% solve(covmatr)
                                     %*% centerdata[i,]))
}

for(i in 1:dim(data)[1]){
  if(covarianceadjusteddists[i]^2 > qchisq(0.999,7)){
    print(paste(nations[i]," is an outlier according to chisq test at 0.1% significance level"))
  }
}

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
plot(data2[,1],data2[,2],main="Wing length vs tail length of female hook-billed kites",
     xlab=names(data2)[1],ylab=names(data2)[2],pch="+",ylim=c(230,320),xlim=c(170,220))
qqnorm(data2[,1],main="Q-Q Plot of tail length")
qqline(data2[,1])
qqnorm(data2[,2],main="Q-Q Plot of wing length")
qqline(data2[,2])

old <- data3[1:30,1:4]
oldmu <- colMeans(old)
mid <- data3[31:60,1:4]
midmu <- colMeans(mid)
new <- data3[61:90,1:4]
newmu <- colMeans(new)

totmu <- colMeans(data3[,1:4])

betwold<-tcrossprod((oldmu - totmu),(oldmu - totmu))
betwmid<-tcrossprod((midmu - totmu),(midmu - totmu))
betwnew<-tcrossprod((newmu - totmu),(newmu - totmu))

Sold <- var(old)
Smid <- var(mid)
Snew <- var(new)


W <- 29 * (Sold + Smid + Snew)
B <- 30 * (betwold + betwmid + betwnew)
WilkLambda <- det(W) / det(B + W)

weirdexpr <- ((90-4-2)/(4)) * (1-sqrt(WilkLambda)) / sqrt(WilkLambda)
limit <- qf(0.95,8,168)
paste("Wilks Lambda:",signif(WilkLambda,3))
cat("Test statistic value for three groups and one or more variables for\n",
      "this Wilk's Lambda:",signif(weirdexpr,3))
cat("F-statistic for above situation at 5% significance level:",
    signif(limit,3))
p <- 4
g <- 3
m <- p * g * (g - 1) / 2
xbarmatr <- rbind(oldmu,midmu,newmu)

tterm <- -qt(0.05/(2 * m),87)
wterms <- sqrt(diag(W)*(2/30)*(1/87))

diff31 <- xbarmatr[3,] - xbarmatr[1,]
diff32 <- xbarmatr[3,] - xbarmatr[2,]
diff21 <- xbarmatr[2,] - xbarmatr[1,]

simulconfints31 <-c()
simulconfints32 <-c()
simulconfints21 <-c()
for(i in 1:4){
  simulconfints31[c(i,i+4)] <- c(diff31[i]-tterm*wterms[i],
                                 diff31[i]+tterm*wterms[i])
  simulconfints32[c(i,i+4)] <- c(diff32[i]-tterm*wterms[i],
                                 diff32[i]+tterm*wterms[i])
  simulconfints21[c(i,i+4)] <- c(diff21[i]-tterm*wterms[i],
                                 diff21[i]+tterm*wterms[i])
}

  paste("Simultaneous 95% confidence intervals for treatment differences")
  
for(i in 1:4){

  print(paste("between periods 3 and 1 for variable",i,": (",
        signif(simulconfints31[i],3),",",signif(simulconfints31[i+4],3),")"))
  print(paste("between periods 3 and 2 for variable",i,": (",
              signif(simulconfints32[i],3),",",signif(simulconfints32[i+4],3),")"))
  print(paste("between periods 2 and 1 for variable",i,": (",
              signif(simulconfints21[i],3),",",signif(simulconfints21[i+4],3),")"))
}
## NA
