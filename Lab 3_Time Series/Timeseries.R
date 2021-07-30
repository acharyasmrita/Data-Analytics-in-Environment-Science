setwd("C:/Users/sacharya/Documents/DSCI607")

sstoi.pa <-read.table("sstoi_pa.txt", header=T)
mo.sstoi <-length(sstoi.pa$YR)-11
yr.sstoi <-mo.sstoi/12

#converting to time series
sstoi.pa12 <- ts(sstoi.pa[,4], start=c(1951,1),end=c(1996,12),frequency=12)
sstoi.pa12 <- window(sstoi.pa12, start=c(1951,1), end=c(1996,12))

sstoi.pa34 <- ts(sstoi.pa[,10], start=c(1950,1), end=c(1996,12), frequency = 12)
sstoi.pa34 <- window(sstoi.pa34, start=c(1951,1), end=c(1996,12))

par(mfrow=c(2,1))
ts.plot(sstoi.pa12)
ts.plot(sstoi.pa34)
acf(sstoi.pa12)
acf(sstoi.pa34)
spec.pgram(sstoi.pa12, spans=c(3,3,3), demean = T, plot=T,xlim=c(0,2))
spec.pgram(sstoi.pa34, spans=c(3,3,3), demean = T, plot = T,xlim=c(0,2))

# 
TB.df <- read.table("TB-Flow.csv", sep=",",header=T)
TB.df <- TB.df[complete.cases(TB.df),]
names(TB.df) <-c("TBday","TBQ")
attach(TB.df)

qtb <- ts(TBQ, freq=365, start=1952, end=2010)
ts.plot(qtb)
acf(qtb,lag.max = 365)
spec.pgram(qtb,spans=c(3,3,3),demean=T,plot=T,xlim=c(0,2),ylim=c(10^6,10^8))
