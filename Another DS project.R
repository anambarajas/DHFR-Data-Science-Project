library(datasets)
library(skimr)
library(Seurat)
library (caret)
library(ggplot2)
library(lattice)
library(kernlab)
library(e1071)

data("dhfr")


#head/tail of data

head(dhfr,4)
tail(dhfr, 4)

#checking a summary of data

summary(dhfr)

#checking total missing data

sum(is.na(dhfr))

#larger summary

skim(dhfr)

#see summary for group of species rather than dependent variables

dhfr %>%
  
  dplyr::group_by(Y) %>%
  skim()

#Panel plots
plot(dhfr, col = "purple")

#Scatterplot

plot(dhfr$moeGao_Abra_R, dhfr$moe2D_a_heavy, col = dhfr$Y, xlab = "Abra R", ylab =  "Heavy")

#Histogram

hist(dhfr$moeGao_Abra_R, col = "Purple")

#Plot with features

featurePlot(x = dhfr[,2:21], 
            y = dhfr$Y, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))

#Classification model

set.seed(100)

TrainingIndex <- createDataPartition(dhfr$Y, p=0.8, list = FALSE)
TrainingSet <- dhfr[TrainingIndex,] # Training Set
TestingSet <- dhfr[-TrainingIndex,] # Test Set

plot(TrainingSet, TestingSet)

#Build training model

Model <- train(Y ~ ., data = TrainingSet,
               method = "svmPoly",
               na.action = na.omit,
               preProcess=c("scale","center"),
               trControl= trainControl(method="none"),
               tuneGrid = data.frame(degree=1,scale=1,C=1)
)

# Build CV model, training leaving out one group of data and predicting it,
#this is done to every k group od data, in this example they use 10 folds,
#so divide 120/10

Model.cv <- train(Y ~ ., data = TrainingSet,
                  method = "svmPoly",
                  na.action = na.omit,
                  preProcess=c("scale","center"),
                  trControl= trainControl(method="cv", number=10),
                  tuneGrid = data.frame(degree=1,scale=1,C=1)
)


# Apply model for prediction

Model.training <-predict(Model, TrainingSet) # Apply model to make prediction on Training set
Model.testing <-predict(Model, TestingSet) # Apply model to make prediction on Testing set
Model.cv <-predict(Model.cv, TrainingSet) # Perform cross-validation

# Model performance (Displays confusion matrix and statistics)

Model.training.confusion <-confusionMatrix(Model.training, TrainingSet$Y)
Model.testing.confusion <-confusionMatrix(Model.testing, TestingSet$Y)
Model.cv.confusion <-confusionMatrix(Model.cv, TrainingSet$Y)

print(Model.training.confusion)
print(Model.testing.confusion)
print(Model.cv.confusion)

# Feature importance

Importance <- varImp(Model)
plot(Importance, top = 5)
plot(Importance, col = "red")



