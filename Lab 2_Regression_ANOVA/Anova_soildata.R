getwd()

#install.packages(c("ggplot2","ggpubr","tidyverse","broom","AICmodvag"))
#library(ggplot2)  #ggplot2 is a data visualization package
#library(ggpubr)   #ggpubr for creating easily publication ready plots
#library(tidyverse)  #tidyverse for data manipulation and visualization
#library(broom)   #broom for printing a nice summary of statistical tests as data frames
#library(AICcmodavg)   #function to implement model selection and multimodel inference based
#on Akaike's information criterion (AIC) and the second-order (AICc), as well as their
#quasi-likelihood counterparts (QAIC, QAICC) from various model object classes

#Step 1: Load the data into R
setwd("C:/Users/sacharya/Documents/DSCI607")

read.csv("FESTUCA.csv")
summary("FESTUCA.csv")

soil.data<-read.csv("FESTUCA.csv")
summary(soil.data)

#Step 2: Perform ANOVA test
#One-way ANOVA using Weight and pH  
#H0=there is no difference in average weight and pH; H1=there is a difference in average weight and pH
one.way <-aov(Weight~pH, data = soil.data)
summary(one.way)

#compute a set of confidence intervals on the differences between the means of the levels of a factor with the specified family-wise probability of coverage.
TukeyHSD(one.way)


#Two way ANOVA using Weight with pH and Calluna
#Two way ANOVA without interaction of two independent variables

two.way<- aov(Weight ~ pH + Calluna, data = soil.data)
summary(two.way)

#Two way ANOVA with interaction of two independent variables (Adding interactions between variable)

interaction <- aov(Weight ~ pH*Calluna, data = soil.data)
summary (interaction)
summary(soil.data)

#Step 3: setting models for three types of test
model.set <- list(one.way, two.way, interaction)
model.names <- c("one.way","two.way","interaction")

#make a table for models
aictab(model.set, modnames = model.names)

#Step 5: Check for homosedasticity
#To check wether the model fits the assumption of homosedasticity, look at the model plot the best fit model

par(mfrow=c(2,2))
plot(interaction)

#Step 6: Do a post-hoc test

tukey.interaction <- TukeyHSD(interaction)
tukey.interaction

#Step 7: Plot the results in a graph
#Find the groupwise differences
#To do this we can run another ANOVA + TukeyHSD test

#Plot the raw data
interaction.plot <- ggplot(soil.data, aes(x = pH, y = Weight, group = Calluna)) + geom_point(cex = 1.5, pch = 1.0, position = position_jitter(w = 0.1, h = 0))
interaction.plot
