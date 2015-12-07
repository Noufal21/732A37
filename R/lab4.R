diag(smpcovmatr) <- c( 1106,2382,2136,.016,70.56)
diag2 <- c(396.7,1143,2.189,.216)
for(i in 1:4){
smpcovmatr[i,i+1] <- diag2[i]
}
diag3 <- c(108.4,-.214,-20.84)
for(i in 1:3){
  smpcovmatr[i,i+2] <- diag3[i]
}
diag4 <-c(.787,-23.96)
for(i in 1:2){
  smpcovmatr[i,i+3] <- diag4[i]
}
smpcovmatr[1,5] <- c(26.23)
smpcovmatr[lower.tri(smpcovmatr)] = t(smpcovmatr)[lower.tri(smpcovmatr)]


