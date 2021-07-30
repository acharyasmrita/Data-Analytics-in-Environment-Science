setwd("C:/Users/sacharya/Documents/DSCI607")

rawdata<-matrix(scan("unif100marked-geoEAS (1).txt", skip=5),ncol=3,byrow=T)
xyv=as.data.frame(rawdata)
names(xyv) <- c("x","y","z" )

xyv 

install.packages("spatstat")
install.packages("car")
install.packages("sgeostat")

install.packages("C:/Users/sacharya/Downloads/seeg_1.0.tar.gz", repos=NULL, type="source")


library(seeg)

xyv<-scan.geoeas.ppp("unif100marked-geoEAS (1).txt")
xyv.v <- vario(xyv,num.lags=10,type='isotropic', maxdist=0.45)
xyv.v
var(xyv$v)
model.semivar.cov(var=xyv.v,nlags = 10,n0=0.04,c0=0.088,a=0.1)
