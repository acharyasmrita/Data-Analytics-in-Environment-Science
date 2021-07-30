setwd("C:/Users/sacharya/Documents/DSCI607")
library("FactoMineR")
library("factoextra")
library("ggplot2")

watershed <- read.table("watsheds.txt", header = T)
watershed.active <- watershed[1:28, 2:9]
watershed.active 
#select first two columns and all the rows to look subset data

head(watershed.active[,1:2],2)

#PCA (x, scale.unit = TRUE, ncp = 5, graph = TRUE)

was.pca <- PCA(watershed.active,scale.unit = TRUE, graph = FALSE)
print(was.pca)

eig.val <- get_eigenvalue(was.pca)
eig.val

wa.pca <-princomp(watershed.active, cor = T)
summary(wa.pca)

#It takes 4 components to explain more than 90% of variance i.e. comp.1 to 4
#plotting loadings, variances and biplots
library(ggplot2)

plot(wa.pca)
barplot(loadings(wa.pca),beside = T)
par(mfrow=c(1,1))
win.graph();biplot(wa.pca, choices=c(1,2))
win.graph();biplot(wa.pca, choices=c(1,3))
win.graph();biplot(wa.pca, choices=c(1,4))
win.graph();biplot(wa.pca, choices=c(2,3))
win.graph();biplot(wa.pca, choices=c(3,4))

fviz_eig(was.pca, addlables = TRUE, ylim = c(0,70))
#chose first two with higher cumulative variance percentage
var <- get_pca_var(was.pca)
var
ind<-get_pca_ind(was.pca) #for individual datasets
ind

#coordinates
head(var$coord) #cordinate information

#cos2 : quality on the factor map
head(var$cos2)

#contributions to the principal components
head(var$contrib)
head(var$coord.3)

#to show variables

fviz_pca_var(was.pca, col.var = "black") #visualize the result of the variable information
library("corrplot")
corrplot(var$cos2, is.corr=FALSE) #darker one ligand represents high quality to represent correlation to contribute to principal component
fviz_cos2(was.pca, choice = "var", axes = 1:2)



#the function dimdesc() [in FactoMineR], for dimension description,
#can be used to identify the most significant associated variables
#with a given principal component it can be used as follow:

was.desc<- dimdesc(was.pca, axes = c(1,2), proba = 0.05)
#description of dimension 1
was.desc$Dim.1

#description of dimension 2
was.desc$Dim.2

#to get individual information 
ind<- get_pca_ind(was.pca)
ind

#coordinates of individuals
head(ind$coord)

#quality of individuals
head(ind$cos2)

#to show contributions of individuals
head(ind$contrib)

fviz_pca_ind(was.pca)

fviz_pca_ind(was.pca, col.ind = "cos2", gradient.cols=c("#00AFBB", "#E7B800", "#FC4E07"),
             repel= TRUE)  #avoid text overlapping (slow if many points)

fviz_cos2(was.pca, choice = "ind")
fviz_contrib(was.pca, choice = "ind", axes = 1:2)

fviz_pca_biplot(was.pca, repel = TRUE, 
                col.var = "#2E9FDF",
                col.ind = "#696969")
fviz_pca_var(was.pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

