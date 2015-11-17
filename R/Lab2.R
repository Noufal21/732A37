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
head(data)
head(data2)
head(data3)
