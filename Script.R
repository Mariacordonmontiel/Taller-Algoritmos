library(readr)
library(e1071)
library(tidyverse)
library(caret)
library(dplyr)
library(MASS)
library(klaR)
library(mda)
library(earth)
library(e1071)
library(rpart)
library(rpart.plot)
library(randomForest)
library(mlbench)
library(randomForest)
library(ggplot2)
library(mlbench)
library(gbm)
library(ggplot2)

# Cargo y preparo la base de datos
classes <- read_delim("classes.csv", delim = ";",
                      escape_double = FALSE, col_names = FALSE, 
                      trim_ws = TRUE)
colnames(classes) <- c("Sample","Class")
genes <- read_delim("gene_expression.csv", 
                    delim = ";", escape_double = FALSE, col_names = FALSE,
                    trim_ws = TRUE)
nombres <- readLines("column_names.txt")
colnames(genes)<-nombres
data <- cbind(classes,genes)
data$Class <- factor(data$Class)
data <- data[,-1]
zero_columns <- sapply(data, function(col) all(col == 0))
data <- data[, !zero_columns]
str(data)


set.seed(123)

# Divido el conjunto de datos de entrenamiento y de prueba
trainIndex <- createDataPartition(data$Class, p=0.7,list=F)
trainData <- data[trainIndex,]
testData <- data[-trainIndex,]

########### K vecinos mas cercanos (k-NN) 
# Estandarizo los datos
train_scaled <- trainData %>%
  mutate(across(where(is.numeric),scale))
test_scaled <- testData %>%
  mutate(across(where(is.numeric),scale))

# Busco el valor máximo de k con validacion cruzada
control <- trainControl(method = "cv", number=10)
tuneGrid <- expand.grid(.k=1:10)
knn_model <- train(Class~., data=train_scaled, method="knn", trControl=control, tuneGrid=tuneGrid)
print(knn_model)
plot(knn_model)
# Obtengo el valor maximo de k
optimal_k <- knn_model$bestTune$k
print(paste("El número óptimo de K es:", optimal_k))

# En este caso el valor optiomo de k es 7

# Evaluo el modelo, primero hago las prediciones y despues creo la matriz para comprobar las predicciones de ese modelo
predictionsk <- predict(knn_model, test_scaled)
confusion_Matrixknn<-confusionMatrix(predictionsk, test_scaled$Class)
confusion_Matrixknn
accuracyk <- sum(predictionsk==test_scaled$Class)/nrow(test_scaled)
print(paste("Precision del modelo k-NN:", accuracyk))

# Guardo los resultados para crear una tabla posteriormente
accuracyknn<-confusion_Matrixknn$overall["Accuracy"]
accuracyknn
kappaknn <- confusion_Matrixknn$overall["Kappa"]
kappaknn


########### Arbol de decision
# Entreno el modelo
model<-rpart(Class~., data=train_scaled, method="class")
# Hago predicciones y evalulo esas predicciones
predictionsAr <- predict(model, newdata=test_scaled, type="class")
confusion_matrixAD <- confusionMatrix(predictionsAr, test_scaled$Class)
confusion_matrixAD
rpart.plot(model, main="Arbol de Decisión para el dataset")

# Guardo los resultados para crear una tabla posteriormente
accuracyAD <- confusion_matrixAD$overall["Accuracy"]
accuracyAD
kappaAD <- confusion_matrixAD$overall["Kappa"]
kappaAD


########### Naive Bayes
# Entreno el modelo
modelNB <- naiveBayes(Class~., train_scaled)
modelNB
# Hago predicciones y evalulo esas predicciones
predictionsNB <- predict(modelNB, test_scaled)
confusion_matrixNB <- confusionMatrix(predictionsNB, test_scaled$Class)
print(confusion_matrixNB)

# Guardo los resultados para crear una tabla posteriormente
accuracyNB <- confusion_matrixNB$overall["Accuracy"]
accuracyNB
kappaNB <- confusion_matrixNB$overall["Kappa"]
kappaNB

# Creo una tabla donde poder comparar los resultados obtenidos de todos los modelos utilizados
resultados <- data.frame(
  Metodo = c("K vecinos más cercanos (k-NN)",
             "Árbol de decisión",
             "Naive Bayes"),
  Accuracy = c(accuracyknn,
               accuracyAD,
               accuracyNB),
  Kappa = c(kappaknn,
            kappaAD,
            kappaNB)
)
print(resultados)





