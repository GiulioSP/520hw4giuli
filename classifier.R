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
uromolDT <- data.table(uromol)
knowlesDT <- data.table(knowles)
#uromol_gini <- gini_impurities(uromolDT, wide=TRUE)
#knowles_gini <- gini_impurities(knowlesDT, wide=TRUE)

#write_csv(uromol_gini, path = "uromol_gini.csv")
#write_csv(knowles_gini, path = "knowles_gini.csv")

# Dropping columns not present in both datasets, as it makes testing invalid
colsIntersect <- intersect(colnames(uromol), colnames(knowles))
exprColsIntersect <- intersect(colnames(uromol$exprs), colnames(knowles$exprs))

knowlesV1 <- knowles[, colsIntersect]
knowlesV1$exprs <- knowlesV1$exprs[, exprColsIntersect]

uromolV1 <- uromol[, colsIntersect]
uromolV1$exprs <- uromolV1$exprs[, exprColsIntersect]

# Drop rows with N/A for key columns
knowlesV1 <- drop_na(knowlesV1, "Progression", "Recurrence")
uromolV1 <- drop_na(uromolV1, "Progression", "Recurrence")

# Drop unnecessary data
knowlesV1 <- subset(knowlesV1, select = -c(PFS_time., Tumor.stage, Tumor.grade))
uromolV1 <- subset(uromolV1, select = -c(PFS_time., Tumor.stage, Tumor.grade))

# Feature Selection
# Make a correlation matrix of RNA expression data from training dataset
correlationMatrix <- cor(uromolV1$exprs)

# Sum the absolute values across columns to evaluate the overall correlation between genes
corrSum <- apply(correlationMatrix,1, function(x) sum(abs(x)))

# Select the 150 genes with the smallest absolute sum value as they are likely most informative for the classifier 
#genes <- names(head(sort(corrSum),500))
genes <- names(head(sort(corrSum),150))

knowlesV2 <- knowlesV1
uromolV2 <- uromolV1
knowlesV2$exprs <- knowlesV2$exprs[, genes]
uromolV2$exprs <- uromolV2$exprs[, genes]

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
#uromolV2 <- uromolV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))
#knowlesV2 <- knowlesV2 %>% mutate(RFS_category = ifelse(Recurrence == 0, "No Recurrence", ifelse(RFS_time >= 25 , "long RFS", "short RFS")))


# Make factor variables be read as factors
knowlesV2$UROMOL2021.classification <- gsub("_", " ", knowlesV2$UROMOL2021.classification)
uromolV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")] <- lapply(uromolV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")], as.factor)
knowlesV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")] <- lapply(knowlesV2[, c("Progression", "Recurrence", "Sex", "Concomitant.CIS", "BCG", "UROMOL2021.classification")], as.factor)


# removing RFS_time since its definition is different for knowles dataset when recorrence == 0 and # replacing RFS_time with 0 whenever there is no recurrence, since it is N/A for cases with no recurrence in knowles dataset.
knowlesV2 <- knowlesV2 %>% mutate(RFS_time = ifelse(Recurrence == 0, 0, RFS_time))
knowlesV4 <- drop_na(knowlesV2)
uromolV4 <- drop_na(uromolV2)
knowlesV4RFS_time <- knowlesV4$RFS_time
uromolV4RFS_time <- uromolV4$RFS_time
knowlesV4 <- subset(knowlesV4, select = -RFS_time)
uromolV4 <- subset(uromolV4, select = -RFS_time)



#uromolV4 %>% group_by(RFS_category) %>% count()
#knowlesV4 %>% group_by(RFS_category) %>% count()

# Train classifier
fit <- train(Recurrence~., data = uromolV4, method = "rf")
pred <- predict(fit, knowlesV4)
confusionMatrix(pred, knowlesV4$Recurrence)
model <- data.frame("pred" = pred, "true" = knowlesV4$Recurrence)
write_csv(model, file = "model2_rf.csv")
#prSummary(knowlesV4, lev = levels(knowlesV4$Recurrence))
# TODO ADD Kaplan-Meier curves
model$time <- knowlesV4RFS_time
model1 <- subset(model, select = c("time","pred"))
model1$source <- rep("pred", length(model$pred))
names(model1)[names(model1) == 'pred'] <- 'event'
model2 <- subset(model, select = c("time","true"))
model2$source <- rep("true", length(model$pred))
names(model2)[names(model2) == 'true'] <- 'event'

survModel <- rbind(model1, model2)

surObj <- Surv(time = survModel$time, event = survModel$event)
survivalfit <- survfit(surObj ~ source, data = survModel)
ggsurvplot(survivalfit)

survfitted <- survival::survfit(Surv(survModel$time, event = survModel$event) ~ source, data = survModel)
tidy_surv <- tidy(survfitted)

ggplot(tidy_surv, aes(time, true)) + geom_line() 


#ga_ctrl <- gafsControl(functions = caretGA)
#rf_ga <- gafs(x = uromolV4, y = uromolV4$RFS_category, iters = 20, gafsControl = ga_ctrl, method = "rpart")

#control <- trainControl(method="repeatedcv", number=10, repeats=3)

#correlationMatrix <- cor(uromolV4$exprs)

trellis.par.set(caretTheme())
cal_obj <- calibration(Class ~ Recurrence, data = pred, cuts = 13)
plot(cal_obj, type = "l", auto.key = list(columns = 3,
                                          lines = TRUE,
                                          points = FALSE))



library(survival)
library(survminer)
library(dplyr)
data(ovarian)
glimpse(ovarian)
surv_object <- Surv(time = ovarian$futime, event = ovarian$fustat)
surv_object 
