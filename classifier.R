# Installs and Imports of packages
#install.packages("mltools")
#install.packages("tidyverse", dependencies =TRUE)
#install.packages("curl", dependencies =TRUE)
#install.packages("class")
#install.packages("ggfortify")
#install.packages("caTools")
#install.packages("caret")


library(caret)
library(mltools)
library(tidyverse)
library(data.table)
library(class)
library(ggfortify)
library(caTools)

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
knowlesV2 <- subset(knowlesV2, select = -c(PFS_time., Tumor.stage, Tumor.grade))
uromolV2 <- subset(uromolV2, select = -c(PFS_time., Tumor.stage, Tumor.grade))

# Visualize expression data as PCA
#uromolV3 <- uromolV2
#knowlesV3 <- knowlesV2

#knowlesV3$ProgressionBoolean <- as.logical(knowlesV3$Progression)
#knowlesV3$RecurrenceBoolean <- as.logical(knowlesV3$Recurrence)

#uromolV3$ProgressionBoolean <- as.logical(uromolV3$Progression)
#uromolV3$RecurrenceBoolean <- as.logical(uromolV3$Recurrence)

#pcaUromol <- prcomp(uromolV3$exprs, scale. = TRUE)
#pcaKnowles <- prcomp(knowlesV3$exprs, scale. = TRUE)

#autoplot(pcaUromol, data = uromolV3, colour = 'ProgressionBoolean') + ggtitle("Uromol PCA labelled by Progression")
#autoplot(pcaUromol, data = uromolV3, colour = 'RecurrenceBoolean') + ggtitle("Uromol PCA labelled by Recurrence")
#autoplot(pcaKnowles, data = knowlesV3, colour = 'ProgressionBoolean') + ggtitle("Knowles PCA labelled by Progression")
#autoplot(pcaKnowles, data = knowlesV3, colour = 'RecurrenceBoolean') + ggtitle("Knowles PCA labelled by Recurrence")
#autoplot(pcaUromol, data = uromolV3, colour = 'UROMOL2021.classification') + ggtitle("Uromol PCA labelled by UROMOL2021 classification")
#autoplot(pcaKnowles, data = knowlesV3, colour = 'UROMOL2021.classification') + ggtitle("Knowles PCA labelled by UROMOL2021 classification")

# Visualize RFS time
# remove rows with no recurrence, since the RFS in those cases is due to followup rate, not prognosis
#uromolV3 <- subset(uromolV3, uromolV3$Recurrence==1)
#knowlesV3 <- subset(knowlesV3, knowlesV3$Recurrence==1)

#ggplot(uromolV3, aes(RFS_time)) + geom_density() + ggtitle("Density of uromol RFS time")
#ggplot(knowlesV3, aes(RFS_time)) + geom_density() + ggtitle("Density of knowles RFS time")

# Prepare categories of RFS to be used as outcome variable for classification
uromolV2 <- uromolV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))
knowlesV2 <- knowlesV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))

# Make factor variables be read as factors
uromolV2$RFS_category <- factor(uromolV2$RFS_category)
uromolV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "UROMOL2021.classification",  "RFS_category")] <- lapply(uromolV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "UROMOL2021.classification",  "RFS_category")], factor)
knowlesV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "UROMOL2021.classification",  "RFS_category")] <- lapply(knowlesV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "UROMOL2021.classification",  "RFS_category")], factor)


#replacing N/A entries with 0 in RFS_time for cases with no recurrence in knowles dataset. 
knowlesV2 <- knowlesV2 %>% mutate(RFS_time = ifelse(Recurrence == 0, 0, RFS_time))

knowlesV4 <- drop_na(knowlesV2)
uromolV4 <- drop_na(uromolV2)


uromolV4 %>% group_by(RFS_category) %>% count()
knowlesV4 %>% group_by(RFS_category) %>% count()

# Train classifier
#knnPred <- knn(train = uromolV4$exprs, test= knowlesV4$exprs, cl = uromolV4$RFS_category, k=4)
#cm <- table(knowlesV4$RFS_category, knnPred)
#cm

#fit <- train(RFS_category~., data = uromolV4, method = "pls", preProc = c("center", "scale"))
fit <- train(RFS_category~., data = uromolV4, method = "nbDiscrete")
fit <- train(RFS_category~., data = uromolV4, method = "awnb")
fit <- train(RFS_category~., data = uromolV4, method = "regLogistic")
fit <- train(RFS_category~., data = uromolV4, method = "logitBoost")
fit <- train(RFS_category~., data = uromolV4, method = "LMT")
fit <- train(RFS_category~., data = uromolV4, method = "regLogistic")
fit <- train(RFS_category~., data = uromolV4, method = "rpart")
fit <- train(RFS_category~., data = uromolV4, method = "rpart2")

