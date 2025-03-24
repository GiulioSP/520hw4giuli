# Installs and Imports of packages
#install.packages("mltools")
#install.packages("tidyverse", dependencies =TRUE)
#install.packages("curl", dependencies =TRUE)
#install.packages("class")
#install.packages("ggfortify")

library(mltools)
library(tidyverse)
library(data.table)
library(class)
library(ggfortify)

# Importing data
knowles <- readRDS("/home/giuliosp/Desktop/Joy Lab/520 hw4/knowles_matched_TaLG_final.rds")
uromol <- readRDS("/home/giuliosp/Desktop/Joy Lab/520 hw4/UROMOL_TaLG.teachingcohort.rds")

# Note any relationships of correspondence or parent-child between columns using gini impurities
uromolDT <- data.table(uromol)
knowlesDT <- data.table(knowles)
uromol_gini <- gini_impurities(uromolDT, wide=TRUE)
knowles_gini <- gini_impurities(knowlesDT, wide=TRUE)

# Dropping columns not present in both datasets, as it makes testing invalid
colsIntersect <- intersect(colnames(uromol), colnames(knowles))
exprColsIntersect <- intersect(colnames(uromol$exprs), colnames(knowles$exprs))

knowlesV1 <- knowles[, colsIntersect]
knowlesV1$exprs <- knowlesV1$exprs[, exprColsIntersect]

uromolV1 <- uromol[, colsIntersect]
uromolV1$exprs <- uromolV1$exprs[, exprColsIntersect]

# Drop unnecessary data



# Visualize expression data as PCA
pca_res <- prcomp(uromolV1$exprs, scale. = TRUE)

autoplot(pca_res)
autoplot(pca_res, data = uromolV1$exprs, colour = 'Progression')


# Prepare categorical data as factors







# Train classifier
knnres <- knn()
head(knnres)

