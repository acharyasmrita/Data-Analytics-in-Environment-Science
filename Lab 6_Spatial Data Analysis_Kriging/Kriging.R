setwd("C:/Users/sacharya/Documents/DSCI607")
library(sgeostat)
library(ggplot2)
data("maas")
mass.point<- point(maas)
plot.point(mass.point,v='zinc',xlab='Easting',ylab='Northing',legend.pos=2,pch=c(21:24),main="",cex=0.7)

library(seeg)
library(car)
mass.v <-vario(maas,num.lags=10,type='isotropic',maxdist=2000)
var(maas$zinc)
m.mass.v<-model.semivar.cov(var=mass.v,nlags=10,n0=50000,c0=150000,a=1000)
maas.vsph<-fit.variogram(model="spherical",mass.v,nugget=50000,sill=150000,range=1000,iterations=30)
maas.ok<-Okriging(maas,maas.vsph,step=100,maxdist=1000)
maas.ok
library(seeg)
plotkriged(maas, maas.ok,outpdf="dataset-kriged.pdf")



#solving matrices
install.packages("matlib")
library(matlib)
b=matrix(c(0.8,0.01,0.01,1,0.01,0.8,0.04,1,0.01,0.04,0.8,1,1,1,1,0), 4, 4)
a=c(0.56,0.09,0.09,1)
Solve(b, a, fractions = FALSE)

