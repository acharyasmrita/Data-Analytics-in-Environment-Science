
setwd("C:/Users/sacharya/Documents/DSCI607")

#read.table for txt files
Data<-read.table("exercise 4 data.txt",header=T) 
View(Data)

scatter.smooth(x=Data$x, y=Data$y, main="y ~ x")
plot(Data$x,Data$y)
lm(Data$y~Data$x)
lm(formula = Data$y~Data$x)
abline(lm(Data$y ~ Data$x))
confint.lm(lm(Data$y ~Data$x), 0.05)
par(mfrow=c(2,2));plot(lm(Data$y~Data$x))

#total standard error, SST
sst<-sum((Data$y-mean(Data$y))^2)
sst
re<-lm(Data$y~Data$x)$coeff[1]+lm(Data$y~Data$x)$coeff[2]*Data$x

#total residual errot, SSE
sse<-sum((Data$y-re)^2)

#total explained errot, SSM
ssm<-sst-sse
ssm

#variance
R2<-ssm/sst
R2
mod<-lm(Data$y~Data$x)
mod
summary(mod)
anova(mod)


