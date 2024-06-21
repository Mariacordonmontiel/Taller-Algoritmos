#En primer lugar, se han de unificar los archivos para poder llevar a cabo la
#depuración de los datos.Se importan los diferentes archivos, de tal forma que el 
#txt de nombres de las columnas aparezca como tal en la base de datos 
#gene_expression. 

#Además, se importa classes.csv, un archivo en el que aparece cada sample 
#acompañado del tipo de tumor del que se trata.

library(readr)
library(dplyr)

data <- read_delim("gene_expression.csv",
                   delim = ";", escape_double = FALSE,
                   col_names = read_lines("column_names.txt"),
                   trim_ws = TRUE)

class <- read_delim("classes.csv", delim = ";", escape_double = FALSE, 
                    col_names = c("Sample","Class"),
                    trim_ws = TRUE)

#Se añaden las dos columnas, sample y class, al database original, que se ha
#nombrado como "data".

data$class <- class$Class
data$sample <- class$Sample



#Una vez tenemos el database completo, se procede a la depuración de los datos,
#para lo que se crea un nuevo database únicamente con las variables que han 
#de ser numéricas. Se comprueba que el numero de filas y columnas es correcto: 
#801 filas y 500 columnas.

data_to_clean = select(data, -sample, -class)
print(paste("El conjunto de datos tiene", nrow(data_to_clean), "filas y", 
            ncol(data_to_clean), "columnas."))



#A continuación, se comprueba que efectivamente todos los valores de la base de
#datos son numéricos. Según el resultado del print, no hay variables no 
#numéricas.

non_numeric_columns <- sapply(data_to_clean, function(x) !is.numeric(x))
non_numeric_names <- names(data_to_clean)[non_numeric_columns]
print(non_numeric_names)



#El siguiente paso es la identificación de valores nulos. Se observa que hay 3 
#variables cuyo valor es igual a 0, lo que nos imposibilita el análisis 
#posterior. Por ello, se filtran los datos de tal manera que se eliminen los 
#datos ausentes o nulos. Así, se eliminan las columnas cuyo sumatorio de valores
#sea 0, ya que esto significa que la columna está completamente vacía.

sumas <- colSums(data_to_clean)
columnascero <- names(sumas[sumas==0])
print(columnascero)
data_to_clean <- data_to_clean[, !names(data_to_clean) %in% columnascero]



#Al igual que se comprueba la existencia de valores nulos, también se comprueba
#si existen valores infinitos. En este caso, no hay valores infinitos.

inf_count_per_column <- sapply(data_to_clean, function(x) sum(is.infinite(x)))
total_inf <- sum(inf_count_per_column)
print(paste("Total de valores infinitos en el conjunto de datos:", total_inf))



#Esta base de datos ya depurada, se guarda en una nueva base de datos "data_clean",
#en la que también se añaden las columnas class y sample.

data_clean <- data_to_clean
data_clean$class <- class$Class
data_clean$sample <- class$Sample



#Por último, se procede a la normalización de los datos, es decir, se escala
#cada columna para que tenga media ~0 y desviación estándar ~1. Además, 
#se añaden de nuevo las columnas sample y class.

normalized_data <- as.data.frame(lapply(data_to_clean, scale))
normalized_data$sample <- data$sample
normalized_data$class <- data$class

#Así, se obtienen dos bases de datos depuradas, una normalizada y otra sin
#normalizar.


#Una vez los datos ya han sido depurados, se pueden implementar los métodos de 
#aprendizaje no supervisado y supervisado para el análisis de los tumores y la 
#expresión génica de diversos genes implicados.


############################################################
# APRENDIZAJE NO SUPERVISADO: REDUCCIÓN DE DIMENSIONALIDAD #
############################################################

#PCA (Análisis de componentes principales)

#Se importan las librerías necesarias para el análisis y la representación
#gráfica en 2D.

library(ggplot2)
library(stats)

#Se seleccionan solo las variables numéricas (todas menos las 2 últimas columnas)
#y se realiza el PCA. Posteriormente, se crea un nuevo dataframe en el que se 
#visualicen los dos primeros componentes principales y la variable Class.

data_pca <- normalized_data[, -((ncol(normalized_data)-1):ncol(normalized_data))]
pca_resultado <- prcomp(data_pca)
print(summary(pca_resultado))
data_PCA_2d <- data.frame(pca_resultado$x[,1:2])
data_PCA_2d$Class <- normalized_data$class

#Se grafican los resultados del PCA en función del tipo de tumor (Class)

ggplot(data_PCA_2d, aes(x = PC1, y = PC2, color = Class)) +
  geom_point() +
  ggtitle("PCA - 2 Componentes Principales") +
  xlab("Componente Principal 1") +
  ylab("Componente Principal 2") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "gray95"), 
        plot.title=element_text(hjust=0.5))



#LE (Laplacian Eigenmap)

install.packages("Rdimtools")
library(Rdimtools)

#Se define que los parámetros k-vecinos más cercanos sean el 20% de la muestra,
#y se guarda el resultado del análisis en un nuevo dataframe. Cabe destacar que
#se ha utilizado la base data_clean pero sin las variables no numéricas. Esto 
#se debe a que al hacer el análisis se detectan las variables no numéricas y 
#se imputan como NAs, lo cual no interesa que ocurra.

data_le <- data_clean[, -((ncol(data_clean)-1):ncol(data_clean))]
le.results <- do.lapeig(data_le, type=c("proportion",0.2), weighted=FALSE)
le.df <- data.frame(le.results$Y)


#Se representan gráficamente los resultados del LE según el tipo de tumor (Class).

ggplot(le.df, aes(x = X1, y = X2, color = data_clean$class)) +
  geom_point() +
  labs(title = "LE - 2 Dimensiones", x = 'dim1', y = 'dim2', color = "Class") +
  theme_classic() +
  theme(panel.grid.major = element_line(color = "gray90"), 
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "gray95"), 
        plot.title=element_text(hjust=0.5))
