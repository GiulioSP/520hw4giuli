# Installs and Imports of packages
install.packages("mltools")
install.packages("tidyverse", dependencies =TRUE)
install.packages("curl", dependencies =TRUE)
install.packages("class")
install.packages("ggfortify")

library(mltools)
library(tidyverse)
library(data.table)
library(class)
library(ggfortify)

# Importing data
wd <- getwd()
knowles <- readRDS(paste(wd, "/knowles_matched_TaLG_final.rds", sep=""))
uromol <- readRDS(paste(wd, "/UROMOL_TaLG.teachingcohort.rds", sep=""))

# Note any relationships of correspondence or parent-child between columns using gini impurities
uromolDT <- data.table(uromol)
knowlesDT <- data.table(knowles)
uromol_gini <- gini_impurities(uromolDT, wide=TRUE)
knowles_gini <- gini_impurities(knowlesDT, wide=TRUE)

write_csv(uromol_gini, path = "uromol_gini.csv")
write_csv(knowles_gini, path = "knowles_gini.csv")

# Dropping columns not present in both datasets, as it makes testing invalid
colsIntersect <- intersect(colnames(uromol), colnames(knowles))
exprColsIntersect <- intersect(colnames(uromol$exprs), colnames(knowles$exprs))

knowlesV1 <- knowles[, colsIntersect]
knowlesV1$exprs <- knowlesV1$exprs[, exprColsIntersect]

uromolV1 <- uromol[, colsIntersect]
uromolV1$exprs <- uromolV1$exprs[, exprColsIntersect]

# Drop rows with N/A for key columns
knowlesV2 <- drop_na(knowlesV1, "Progression", "Recurrence")
uromolV2 <- drop_na(uromolV1, "Progression", "Recurrence")

# Drop unnecessary data
knowlesV2 <- subset(knowlesV2, select = -PFS_time.)
uromolV2 <- subset(uromolV2, select = -PFS_time.)

# Visualize expression data as PCA
knowlesV2$ProgressionBoolean <- as.logical(knowlesV2$Progression)
knowlesV2$RecurrenceBoolean <- as.logical(knowlesV2$Recurrence)

uromolV2$ProgressionBoolean <- as.logical(uromolV2$Progression)
uromolV2$RecurrenceBoolean <- as.logical(uromolV2$Recurrence)

pcaUromol <- prcomp(uromolV2$exprs, scale. = TRUE)
pcaKnowles <- prcomp(knowlesV2$exprs, scale. = TRUE)

autoplot(pcaUromol, data = uromolV2, colour = 'ProgressionBoolean') + ggtitle("Uromol PCA labelled by Progression")
autoplot(pcaUromol, data = uromolV2, colour = 'RecurrenceBoolean') + ggtitle("Uromol PCA labelled by Recurrence")
autoplot(pcaKnowles, data = knowlesV2, colour = 'ProgressionBoolean') + ggtitle("Knowles PCA labelled by Progression")
autoplot(pcaKnowles, data = knowlesV2, colour = 'RecurrenceBoolean') + ggtitle("Knowles PCA labelled by Recurrence")
autoplot(pcaUromol, data = uromolV2, colour = 'UROMOL2021.classification') + ggtitle("Uromol PCA labelled by UROMOL2021 classification")
autoplot(pcaKnowles, data = knowlesV2, colour = 'UROMOL2021.classification') + ggtitle("Knowles PCA labelled by UROMOL2021 classification")


# Visualize RFS time
# remove rows with no recurrence, since the RFS in those cases is due to followup rate, not prognosis
uromolV3 <- subset(uromolV2, uromolV2$Recurrence==1)
knowlesV3 <- subset(knowlesV2, knowlesV2$Recurrence==1)

ggplot(uromolV3, aes(RFS_time)) + geom_density() + ggtitle("Density of uromol RFS time")
ggplot(knowlesV3, aes(RFS_time)) + geom_density() + ggtitle("Density of knowles RFS time")

# Prepare categories of RFS to be used as outcome variable for classification
uromolV2 <- uromolV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))
knowlesV2 <- knowlesV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))

#replacing N/A entries with 0 in RFS_time for cases with no recurrence in knowles dataset. 
knowlesV2 <- knowlesV2 %>% mutate(RFS_time = ifelse(Recurrence == 0, 0, RFS_time))

knowlesV4 <- drop_na(knowlesV2, "RFS_category")
uromolV4 <- drop_na(uromolV2, "RFS_category")

uromolV4 %>% group_by(RFS_category) %>% count()
knowlesV4 %>% group_by(RFS_category) %>% count()

# Train classifier
#knowlesV4 <- subset(knowlesV4, select = -RFS_time)
#uromolV4 <- subset(uromolV4, select = -RFS_time)
knnres <- knn(train = uromolV4, test= knowlesV4, cl = train$RFS_category, k=10 )
head(knnres)

