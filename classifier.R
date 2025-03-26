# Installs and Imports of packages
#install.packages("mltools")
#install.packages("tidyverse", dependencies =TRUE)
#install.packages("curl", dependencies =TRUE)
#install.packages("class")
#install.packages("ggfortify")
#install.packages("caTools")
#install.packages("caret")
#install.packages(c("survival", "survminer", "ggplot2"))

library(caret)
library(mltools)
library(tidyverse)
library(data.table)
library(class)
library(ggfortify)
library(caTools)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(broom)

# Importing data
wd <- getwd()
knowles <- readRDS(paste(wd, "/knowles_matched_TaLG_final.rds", sep=""))
uromol <- readRDS(paste(wd, "/UROMOL_TaLG.teachingcohort.rds", sep=""))

# Note any relationships of correspondence or parent-child between columns using gini impurities
#uromolDT <- data.table(uromol)
#knowlesDT <- data.table(knowles)
#uromol_gini <- gini_impurities(uromolDT, wide=TRUE)
#knowles_gini <- gini_impurities(knowlesDT, wide=TRUE)

#write_csv(uromol_gini, path = "uromol_gini.csv")
#write_csv(knowles_gini, path = "knowles_gini.csv")

# Drop rows with N/A for key columns
knowlesV1 <- drop_na(knowles, "Progression", "Recurrence")
uromolV1 <- drop_na(uromol, "Progression", "Recurrence")

# Drop unnecessary data
knowlesV1 <- subset(knowlesV1, select = -c(PFS_time., Tumor.stage, Tumor.grade))
uromolV1 <- subset(uromolV1, select = -c(PFS_time., Tumor.stage, Tumor.grade))

# Feature Selection
# remove genes with near zero variance. This removes ~4k genes 
allGenes <- colnames(uromolV1$exprs)
nzvGenes <- nearZeroVar(uromolV1$exprs, names=TRUE)
uromolV1$exprs <- uromolV1$exprs[, setdiff(allGenes, nzvGenes)]

# Drop columns not present in both datasets, as it makes testing worse. This drops ~10k genes
colsIntersect <- intersect(colnames(uromolV1), colnames(knowlesV1))
exprColsIntersect <- intersect(colnames(uromolV1$exprs), colnames(knowlesV1$exprs))

knowlesV1 <- knowlesV1[, colsIntersect]
knowlesV1$exprs <- knowlesV1$exprs[, exprColsIntersect]

uromolV1 <- uromolV1[, colsIntersect]
uromolV1$exprs <- uromolV1$exprs[, exprColsIntersect]

# Make a correlation matrix of expression data from training dataset
correlationMatrix <- cor(uromolV1$exprs)

# Sum the absolute values across columns to evaluate the overall correlation between genes
corrSum <- apply(correlationMatrix,1, function(x) sum(abs(x)))

# Select the 150 genes with the smallest absolute sum value as they are likely most informative for the classifier 
genes <- names(head(sort(corrSum),500))

knowlesV2 <- knowlesV1
uromolV2 <- uromolV1
knowlesV2$exprs <- knowlesV2$exprs[, genes]
uromolV2$exprs <- uromolV2$exprs[, genes]

# Visualize expression data as PCA
uromolV3 <- uromolV2
knowlesV3 <- knowlesV2

knowlesV3$ProgressionBoolean <- as.logical(knowlesV3$Progression)
knowlesV3$RecurrenceBoolean <- as.logical(knowlesV3$Recurrence)

uromolV3$ProgressionBoolean <- as.logical(uromolV3$Progression)
uromolV3$RecurrenceBoolean <- as.logical(uromolV3$Recurrence)

pcaUromol <- prcomp(uromolV3$exprs, scale. = TRUE)
pcaKnowles <- prcomp(knowlesV3$exprs, scale. = TRUE)

autoplot(pcaUromol, data = uromolV3, colour = 'ProgressionBoolean') + ggtitle("Uromol PCA labelled by Progression")
autoplot(pcaUromol, data = uromolV3, colour = 'RecurrenceBoolean') + ggtitle("Uromol PCA labelled by Recurrence")
autoplot(pcaUromol, data = uromolV3, colour = 'UROMOL2021.classification') + ggtitle("Uromol PCA labelled by UROMOL2021 classification")
autoplot(pcaKnowles, data = knowlesV3, colour = 'ProgressionBoolean') + ggtitle("Knowles PCA labelled by Progression")
autoplot(pcaKnowles, data = knowlesV3, colour = 'RecurrenceBoolean') + ggtitle("Knowles PCA labelled by Recurrence")
autoplot(pcaKnowles, data = knowlesV3, colour = 'UROMOL2021.classification') + ggtitle("Knowles PCA labelled by UROMOL2021 classification")

# Visualize RFS time
# remove rows with no recurrence, since the RFS in those cases is due to followup rate, not prognosis
uromolV3 <- subset(uromolV3, uromolV3$Recurrence==1)
knowlesV3 <- subset(knowlesV3, knowlesV3$Recurrence==1)

ggplot(uromolV3, aes(RFS_time)) + geom_density() + ggtitle("Density of uromol RFS time")
ggplot(knowlesV3, aes(RFS_time)) + geom_density() + ggtitle("Density of knowles RFS time")

# Prepare categories of RFS to be used as outcome variable for classification
#uromolV2 <- uromolV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))
#knowlesV2 <- knowlesV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))

# Make factor variables be read as factors
knowlesV2$UROMOL2021.classification <- gsub("_", " ", knowlesV2$UROMOL2021.classification)
uromolV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")] <- lapply(uromolV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")], as.factor)
knowlesV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")] <- lapply(knowlesV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")], as.factor)


# removing RFS_time since its definition is different for knowles dataset when recorrence == 0 and # replacing RFS_time with 0 whenever there is no recurrence, since it is N/A for cases with no recurrence in knowles dataset.
#max <- max(knowlesV2$RFS_time, na.rm = TRUE)
knowlesV2 <- knowlesV2 %>% mutate(RFS_time = ifelse(Recurrence == 0, 0, RFS_time))
knowlesV4 <- drop_na(knowlesV2)
uromolV4 <- drop_na(uromolV2)
knowlesV4RFS_time <- knowlesV4$RFS_time
uromolV4RFS_time <- uromolV4$RFS_time
knowlesV4 <- subset(knowlesV4, select = -RFS_time)
uromolV4 <- subset(uromolV4, select = -RFS_time)

# Train classifier
trainControl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)

fit <- train(Recurrence~., data = uromolV4, method = "rf", trControl = trainControl, preProcess=c("center", "scale"))
pred <- predict(fit, knowlesV4)
confusionMatrix(pred, knowlesV4$Recurrence)

model <- data.frame("pred" = pred, "true" = knowlesV4$Recurrence)
write_csv(model, file = "model_rf.csv")

# Kaplan-Meier survival curve 
model$time <- knowlesV4RFS_time
survTrue <- survfit(Surv(time = model$time, event = model$true)~ 1, data = model)
survPred <- survfit(Surv(time = model$time, event = model$pred)~ 1, data = model)
ggsurvplot_combine(list(survTrue, survPred), pval = TRUE, conf.int = TRUE, legend.labs = c("True", "Prediction"), palette = c("blue", "red"))
