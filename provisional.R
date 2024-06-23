
#### CARGA Y DEPURACIÓN DE DATOS ####

# Cargamos las librerías necesarias
library(readr)
library(dplyr)
library(factoextra)


# Cargamos los datos

data <- read_delim("gene_expression.csv",
                   delim = ";", escape_double = FALSE,
                   col_names = read_lines("column_names.txt"),
                   trim_ws = TRUE)

class <- read_delim("classes.csv", delim = ";", escape_double = FALSE, 
                    col_names = c("Sample","Class"),
                    trim_ws = TRUE)

#Se añaden las columnas 'Sample' y 'Class' al dataframe original

data$class <- class$Class
data$sample <- class$Sample

# Creamos un dataframe con las variables numéricas, exceptuando las columnas 'Sample' y 'Data'

data_to_clean = select(data, -sample, -class)

# Verificación de valores no numéricos

non_numeric_columns <- sapply(data_to_clean, function(x) !is.numeric(x))
non_numeric_names <- names(data_to_clean)[non_numeric_columns]
print(non_numeric_names)

#Eliminación de datos ausentes = 0

sumas <- colSums(data_to_clean)
columnascero <- names(sumas[sumas==0])
print(columnascero)
data_to_clean <- data_to_clean[, !names(data_to_clean) %in% columnascero]

data_clean <- data_to_clean
data_clean$class <- class$Class
data_clean$sample <- class$Sample

#Por último, se procede a la normalización de los datos (escala)

normalised_data <- as.data.frame(lapply(data_to_clean, scale))

# Añadimos de nuevo las columnas 'Sample' y 'Class'
normalised_data$sample <- data$sample
normalised_data$class <- data$class



###### APRENDIZAJE NO SUPERVISADO: MODELOS DE CLUSTERIZACIÓN ######


#### MODELO K-MEANS ####

# Aplicación de K-Means para encontrar el número óptimo de clusters
set.seed(123) # Establecer semilla para reproducibilidad

# Definir la función para calcular el total within sum of squares (WSS)
wss <- function(k) {
  #Ejecutar K-means y devolver el total withim sum of squares (WSS)
  kmeans(normalised_data[, -c(ncol(normalised_data)-1, ncol(normalised_data))], k, nstart = 25)$tot.withinss
}

# Calcular el total within sum of squares (WSS) para k = 1 a k = 10
k.values <- 1:10
wss_values <- sapply(k.values, wss)

# Graficar el "elbow curve" para encontrar el número óptimo de clusters
elbow_plot <- qplot(k.values, wss_values, geom = "line") +
  ggtitle("Elbow Curve para Encontrar el Número Óptimo de Clusters") +
  xlab("Número de Clusters K") +
  ylab("Total Within Sum of Squares (WSS)")

print(elbow_plot)

# Basado en la gráfica del "elbow curve", elegir el valor de K
optimal_k <- 10

# Aplicar K-Means con el valor óptimo de K
kmeans_result <- kmeans(normalised_data[, -c(ncol(normalised_data)-1, ncol(normalised_data))], centers = optimal_k, nstart = 25)

# Añadir el resultado de los clusters al dataframe original
data$Cluster <- as.factor(kmeans_result$cluster)

# Visualizar los clusters utilizando ggplot2
cluster_plot <- ggplot(data, aes(x = normalised_data[,1], y = normalised_data[,2], color = Cluster)) + 
  geom_point() + 
  ggtitle(paste("Clusters Identificados con K-Means (K =", optimal_k, ")")) + 
  xlab("Primera Dimensión") + 
  ylab("Segunda Dimensión")

print(cluster_plot)



#### CLUSTERIZACIÓN JERÁRQUICA AGLOMERATIVA #####

# Cargar las librerías necesarias
library(cluster)  # Para realizar análisis de clustering
library(factoextra) # Para visualización de análisis factorial y clustering

# Calcular la matriz de distancias
dist_matrix <- dist(normalised_data)

# Realizar el clustering jerárquico aglomerativo
hclust_model <- hclust(dist_matrix, method = 'ward.D')

# Convertir el modelo hclust a un objeto dendrogram
dendro <- as.dendrogram(hclust_model)

# Visualizar el dendrograma para explorar la estructura de clusters
plot(dendro, main = "Dendrograma de Clusterización Jerárquica Aglomerativa",
     xlab = "Índice de Observaciones", ylab = "Distancia")


# Definir el número de clusters deseados y obtener los colores para la visualización
num_clusters <- 7  # Número de clusters que deseamos visualizar
colors <- rainbow(num_clusters)

# Dibujar rectángulos alrededor de los clusters deseados
rect.hclust(hclust_model, k = num_clusters, border = colors)

# Añadir leyenda
legend("topright", legend = paste("Cluster", 1:num_clusters), col = colors, pch = 16)




