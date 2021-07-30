

setwd("C:/Users/sacharya/Documents/DSCI607")
install.packages("spdep")
install.packages("maptools")
install.packages("spatialreg")
library(spdep)
library(maptools)
library(spatialreg)


data(nc.sids)
help(nc.sids) 
nc.sids
names(nc.sids)
sidspolys <- readShapePoly(system.file("shapes/sids.shp", package="maptools"))
plot(sidspolys,axes=T)

#coordinates(), set or retrieve the spatial coordinates. For spatial polygons it returns the centroids.
text(coordinates(sidspolys), labels=row.names(sidspolys), cex=0.6)
# We can plot a map of counties coded according to the value of a variable, say SID79. First, we
# examine the values of SID79 cut the values of the variable in intervals
nc.sids$SID79  #SID deaths, 1979-84

Sid79 <- as.ordered(cut(nc.sids$SID79, breaks=c(0, 5, 10, 15, 20, 30, 40, 50, 60),
                        include.lowest=TRUE))
unclass(Sid79)
#we see how each county is assigned a code from the numbers 1–7. Then we assign a gray tone to
# each interval and plot the map of polygons according to the intervals and color codes selected

cols <- grey(seq(0.3,1,0.1))
plot(sidspolys, col=cols[unclass(Sid79)], border = par("fg"),axes=T)
legend(c(-84,-82 ), c(32, 34), legend=paste("Sids 79", levels(Sid79)), fill=cols, bty="n")

# As explained in the help of spdep, we can also view the map as probabilities of observing these
# values if rates were to follow a Poisson distribution. First, compute the parameter for the Poisson
sids.phat <- sum(nc.sids$SID79) / sum(nc.sids$BIR79)


# apply the cumulative of the Poisson (ppois)
pm <- ppois(nc.sids$SID79, sids.phat*nc.sids$BIR79)
pm.f <- as.ordered(cut(pm, breaks=c(0.0, 0.01, 0.05, 0.1, 0.9, 0.95, 0.99, 1),
                       include.lowest=TRUE))

# as before cut, assign colors and map
cols <- grey(seq(0.3,1,0.1))
plot(sidspolys, col=cols[unclass(pm.f)],border = par("fg"),axes=T)
legend(c(-84,-82 ), c(32, 34), legend=paste("prob.", levels(pm.f)), fill=cols, bty="n")     


# For auto-correlation and auto-regression, we need to build a neighbor structure to be stored as
# object of class nb. Let us find the weighted neighbor matrix. Neighbors can be defined according
# to distance between county seats or amount of shared borders. The weights would represent
# intensity of neighbor relationship (e.g. extent of common boundary, closeness of centroids).
# Assigning weights is critical since spatial correlation and spatial regression models will
# eventually depend on the weights.
# First, let us use the distance. Bind the coordinates of county seats as a 100×2 matrix. This can be
# done using various coordinate systems: eastings-northings, UTM and long-lat

sids.coords <- cbind(nc.sids$east, nc.sids$north)
sids.utm <- cbind(nc.sids$x, nc.sids$y)
sids.lonlat <- cbind(nc.sids$lon, nc.sids$lat)



# Let us use the eastings and northings and then use function dnearneigh to find neighbors (start
# with county seats < 30 miles of each other). This function will exclude the redundant cases when
# a region is neighbor with itself.
sids.nb <- dnearneigh(sids.coords, 0, 30, row.names = rownames(nc.sids))
sids.nb
# you can find the link for each county by addressing the id number, for example for counties 1
# and 2
sids.nb[1:2]


# This means for example that county 1 (Alamance) has neighbors 17 19 32 41 68. Note that selfneighbors
# are excluded.
# However, as you go through sids.nb you realize that counties 28 and 48 have no neighbors at
# this distance. This is more evident if you plot the links
plot(sids.nb, sids.coords)

# where you see some isolated county seats. Increase the neighbor cutoff to 35 miles to see if we
# include these counties

sids.nb <- dnearneigh(sids.coords, 0, 35, row.names = rownames(nc.sids))
plot.nb(sids.nb, sids.coords)


plot(sidspolys, border = "grey")
plot(sids.nb, sids.lonlat, add=T)

# Next, we need weights for neighborhood structure. In this case, a simple calculation of weights is
# the inverse of the distance; i.e. the closest the regions the larger the weight (the closest two                                                                                 regions are, the more intense the neighbor effect).
# We can get distances between each pair

sids.dists <- nbdists(sids.nb, sids.coords)
sids.dists
# and also we can get an intensity of neighborhood using the inverse of the distance applied to all
# pairs

inten <- lapply(sids.dists, function(x) 1/x)
inten

# Next we produce a neighbor matrix with these weights. This is a 100x100 matrix. Style “W” is
# row standardized, i.e., divide by number of neighbors in row. For brevity in this guide, we only
# show the first 4 rows and 10 columns
sids.nbmat <- nb2mat(sids.nb, glist=inten, style="W")
sids.nbmat[1:4,1:10]

# Same result can be put in a list because the function to compute Moran statistics will require a
# list. Shown here are only the first four rows
sids.nblsw <- nb2listw(sids.nb, glist=inten, style="W")
sids.nblsw

# 
# We can build a neighbor object using the polygons instead of centroids. Neighboring relations
# are established based on shared borders instead ot distance between county seats.
sids.nb.pol <- poly2nb(sidspolys, row.names = rownames(nc.sids))
# having similar structure as the sids.nb object we built above using distance between county seats.
# We can examine the differences
diffnb(sids.nb, sids.nb.pol, verbose=TRUE)

# Next, we need weights for neighborhood structure. The weights would represent intensity of
# neighbor relationship (in this case extent of common boundary). Results are put in a list because
# later the function to compute Moran’s I statistic will require a list. Shown here only the first four
# rows
sids.nblsw.pol <- nb2listw(sids.nb.pol, glist=NULL, style="W")
sids.nblsw.pol$weights[1:4]

# Now, with these results we are ready to apply spatial auto-correlation (e.g., using the Moran’s I
#                                                                         statistic) and auto-regression (SAR).
# Apply the Moran’s I spatial auto-correlation test to the SID79 variable of nc.sids data set with
# the neighbor structures built in the previous section. The arguments are the variable, the
# neighborhood structure, and a decision on whether we assume normality or randomization. Let
# us assume normality.

moran.test(nc.sids$SID79, sids.nblsw, randomisation=F)


# Here although the value for the Moran statistic is low and would seem to indicate no spatial
# pattern, the low variance makes the Z value high and consequently the p-value is low enough to
# conclude that there is spatial pattern.
# We have assumed that the mean and variance are the same for all regions. As discussed in
# Kaluzny et al. (1996) the variance of sids rate increases for counties with low birth rates and
# therefore the data needs to be transformed using the Freeman-Tukey (FT) transform of the sids
# rate and multiplied by the square root of births, to achieve constant variance. The FT transform is

ft.SID79 <- sqrt(1000)*(sqrt(nc.sids$SID79/nc.sids$BIR79) +
                          sqrt((nc.sids$SID79+1)/nc.sids$BIR79))
# and then multiply by square root of births
tr.SID79 <- ft.SID79*sqrt(nc.sids$BIR79)
names(tr.SID79) <- rownames(nc.sids)
# So, modify the Moran test as follows
moran.test(tr.SID79, sids.nblsw, randomisation=F)
# 
# This produces an increase in Moran’s I, and consequently an improvement in the p-value.
# Let us try randomization
moran.test(nc.sids$SID79, sids.nblsw, randomisation=T)

# This also represents an improvement in p-value compared to normality. So far, results indicate
# that there is spatial auto-correlation of the Sids rate 1979-1984 in North Carolina counties. We
# can confirm with Monte Carlo simulations (calculating Moran’s I many times.
moran.mc(nc.sids$SID79, sids.nblsw, nsim=1000)
# 
# We also have a low p-value.
# You can repeat with the Geary statistic, by just changing “moran” for "geary".

geary.test(nc.sids$SID79, sids.nblsw, randomisation=T)
# We also have a low p-value.
# You can repeat with the Geary statistic, by just changing “moran” for "geary".
geary.mc(nc.sids$SID79, sids.nblsw, nsim=1000)
# We can plot the spatial auto-correlation as a function of lag distance or spatial correlogram
plot(sp.correlogram(sids.nb, nc.sids$SID79, order = 10, method = "corr", style = "W"))

# In the top panel we can appreciate that values are positively correlated at short lags indicating
# spatial pattern. Another option for method is “I” and this will do Moran’s I at different lags
# (bottom panel).
plot(sp.correlogram(sids.nb, nc.sids$SID79, order = 10, method = "I", style = "W"))

# Spatial auto-correlation depends on the weights, therefore it is good practice to perform several
# runs of Moran and Geary with different weights.
# The spatial auto-regression (SAR) linear model (SLM) object is sarslm. Package spdep provides
# functions lagsarslm for the spatial lag model and errorsarslm for the spatial error model, to
# 180
# perform SAR auto-regression. We will use lagsarslm to estimate the coefficients of a predictor
# and the autocorrelation parameter ρ.
# First, assume that race is related to SIDS rates, and therefore SID is modeled as a function of
# nonwhite birth rates NWBIR.
install.packages("spatialreg")
sids.lag <- lagsarlm(SID79 ~ NWBIR79, data=nc.sids, sids.nblsw,tol.solve=1e-9)
sids.lag
# note that the results include estimates of the coefficients (intercept and slope) and the parameter
# rho for ρ. This model is then used to predict the sids rate by region (county). It was necessary to
# decrease the tolerance to 10-9. By default, it is 10-7. Using summary, we get more information
summary(sids.lag)

# The coefficient estimates have low p-values. A special function is the Likelihood Ratio Test
# (LR) to check whether the parameter rho is nonzero. In this case, the estimate is still relatively
# good with a p-value of 0.08. This tells us that ρ is significantly different from zero. However, the
# Lagrange Multiplier (LM) test for lack of serial correlation of the residuals yields high p-value
# and therefore the null hypothesis of serially correlated residuals cannot be rejected.
# Predicted values and residuals are available as part of the slm object type as fitted and resid. This
# model can be diagnosed much in the same manner that we evaluated linear regression models:
#   scatter plots, qq plots and residual-predicted plots. For example,
lim<- max(nc.sids$SID79, fitted(sids.lag))
split.screen(c(2,1))
screen(1)
par(mar=c(4,4,1,.5),xaxs=r, yaxs=r)
plot(fitted(sids.lag), resid(sids.lag),xlab="Estimated",ylab="Residuals")
abline(h=0)
split.screen(c(1,2), screen = 2)
screen(3)
plot(nc.sids$SID79,
     fitted(sids.lag),xlab="Observed",ylab="Estimated",xlim=c(0,lim),ylim=c(0,lim))
abline(a=0,b=1)
screen(4)
qqnorm(resid(sids.lag))
qqline(resid(sids.lag))


# 
# where we see that there are some outliers but it is generally fine. The residuals do not behave like
# a normal distribution.
# We know that the Freeman-Tukey (FT) transform helps to stabilize the variance of sids. Let us
# perform SAR on the FT transform of SID79 and NWBIR79. We already have tr.SID79. Let us
# calculate the one for NWBIR79 
ft.NWBIR79 <- sqrt(1000)*(sqrt(nc.sids$NWBIR79/nc.sids$BIR79) +
                            sqrt((nc.sids$NWBIR79+1)/nc.sids$BIR79))
# and then multiply by square root of births   
tr.NWBIR79 <- ft.NWBIR79*sqrt(nc.sids$BIR79)
names(tr.NWBIR79) <- rownames(nc.sids)


tr.sids <- data.frame(tr.SID79, tr.NWBIR79)

# Run SAR on the transformed data
sids.tr.lag <- lagsarlm(tr.SID79 ~ tr.NWBIR79, data=tr.sids, sids.nblsw, tol.solve=1e-12)
summary(sids.tr.lag)

# Note that we have improved the estimation of rho using this transform and the LM test of
# residuals yields low p-value, allowing us to reject the null of serially correlated residuals.
par(mfrow=c(2,2))
plot(tr.sids$tr.SID79, fitted(sids.tr.lag))
plot(fitted(sids.tr.lag), resid(sids.tr.lag))
qqnorm(resid(sids.tr.lag))
# As we can see we now have achieved a better behavior of the residuals.
