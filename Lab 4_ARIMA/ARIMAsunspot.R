setwd("C:/Users/sacharya/Documents/DSCI607")

library("forecast")
sunspotyr <- read.table("year-spot1700-2011.txt", header=T)


yrspots <- matrix(scan("year-spot1700-2011.txt",skip=6),ncol=2,byrow=T)
yrspots

#raster time series .rts
yrspots.rts <- ts(yrspots[,2], start=1700, deltat=1)
ts.plot(yrspots.rts,ylab="Number Sunspots", xlab="Years")
winspots.rts <- window(yrspots.rts, start=1700, end=1750)
ts.plot(winspots.rts)
lag.plot(yrspots.rts,lags=11,layout=c(2,3))
acf(yrspots.rts)
acf(yrspots.rts, type="covariance")
var(yrspots[,2])
max(acf(yrspots.rts, type="covariance")$acf)

#AR(9), to predict the sunspots 22 year ahead from 2011 i.e. from 2012-2033

#predict
ar.yrspots <- ar.yw(yrspots.rts)
ar.yrspots$order.max
ar.yrspots$order
yrspots.pred <- predict(ar.yrspots, yrspots.rts, n.ahead = 22)

#predict
up <- yrspots.pred$pred + 2*yrspots.pred$se
low <- yrspots.pred$pred - 2*yrspots.pred$se
min <- min(yrspots.rts, low)
max <- max(yrspots.rts, up)

#plot
ts.plot(yrspots.rts, yrspots.pred$pred, col=1, lty=c(1,2), xlim=c(1700,2033),
        ylim=c(min, max),ylab="X")
lines(up, col=1, lty=3)
lines(low, col=1, lty=3)
legend("top", leg=c("Data", "Pred","Upper & Lower"), lty=c(1,2,3))


ts.plot(yrspots.rts, yrspots.pred$pred, col=1, lty=c(1,2), xlim=c(1960,2033),
        ylim=c(min, max),ylab="X")

lines(up, col=1, lty=3)
lines(low, col=1, lty=3)
legend("top", leg=c("Data", "Pred","Upper & Lower"), lty=c(1,2,3))
ar.yrspots$aic

#Predict sunspots 22 year ahead from year 1980 (1981-2003)

graphics.off()
x.pred <- predict(ar.yrspots, xlim = c(1980), n.ahead=22)
up <- yrspots.pred$pred + 2*yrspots.pred$se
low <- yrspots.pred$pred - 2*yrspots.pred$se

minx <- min(yrspots.rts, low)
maxx <- max(yrspots.rts, up)
ts.plot(yrspots.rts,yrspots.pred$pred, col=1, lty=c(1,2), xlim=c(1960,2003), ylim=c(minx, maxx),ylab="X")
lines(up, col=1, lty=3)
lines(low, col=1, lty=3)
legend("top", leg=c("Data", "Pred","Upper & Lower"), lty=c(1,2,3))


#
library(forecast)
ar.yrspots<-ar(yrspots.rts, order=c(9,1,0))
print(ar.yrspots)
predict_ar$pred[1]
graphics.off()
predict(ar.yrspots, n.ahead=22)
par(mfrow=c(1,1))
ts.plot(yrspots.rts, xlim=c(1970,1980))
AR_forecast <- predict(ar.yrspots, n.ahead=22)$pred
AR_forecast_se <- predict(ar.yrspots, n.ahead=22)$se
points(AR_forecast, type ="l",col = 2, lwd = 3)
points(AR_forecast - 2*AR_forecast_se, type = "l", col= 2, lty = 2, lwd = 2)
points(AR_forecast + 2*AR_forecast_se, type = "l", col= 2, lty = 2, lwd = 2)


#ARIMA(9,1,0)

arima.yrspots<-arima(yrspots.rts, order=c(9,1,0))
yrspots.predaar <-predict(arima.yrspots, n.ahead=22)

#predict
up <-yrspots.predaar$pred + 2*yrspots.pred$se
low<-yrspots.predaar$pred - 2*yrspots.pred$se
min <- min(yrspots.rts, low)
max <- max(yrspots.rts, up)

#plot
ts.plot(yrspots.rts,yrspots.predaar$pred,col=1,lty=c(1,2),xlim=c(1700,2033),ylim=c(min,max),ylab="X")
lines(up,col=1,lty=3)
lines(low,col=1, lty=3)
legend("top",leg=c("Data","Pred","Upper & Lower"), lty=c(1,2,3))

ts.plot(yrspots.rts,yrspots.predaar$pred, col=1,lty=c(1,2),xlim=c(1960,2033),ylim=c(min,max),ylab="X")
lines(up,col=1,lty=3)
lines(low,col=1, lty=3)
legend("top",leg=c("Data","Pred","Upper & Lower"), lty=c(1,2,3))
arima.yrspots$aic

#Predict sunspots 22 year ahead from year 1980 (1981-2003)

library(forecast)
arima.yrspots<-arima(yrspots.rts, order=c(9,1,0))
print(arima.yrspots)
predict_arima$pred[1]
graphics.off()
predict(arima.yrspots, n.ahead=22)

par(mfrow=c(1,1))
ts.plot(yrspots.rts)
arima_fit <- yrspots.rts - residuals(arima.yrspots)
points(arima_fit, type ="l", col=2, lty=2)

ts.plot(yrspots.rts, xlim = c(1960,2003))
AR_forecast <- predict(arima.yrspots, n.ahead=22)$pred
AR_forecast_se <- predict(arima.yrspots, n.ahead=22)$se
points(AR_forecast, type ="l",col = 2, lwd = 3)
points(AR_forecast - 2*AR_forecast_se, type = "l", col= 2, lty = 2, lwd = 2)
points(AR_forecast + 2*AR_forecast_se, type = "l", col= 2, lty = 2, lwd = 2)
