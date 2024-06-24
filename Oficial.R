
###### APRENDIZAJE NO SUPERVISADO: MODELOS DE CLUSTERIZACIÓN ######

#### MODELO K-MEANS ####

# Carga de librerías necesarias
library(readr)
library(dplyr)
library(factoextra)

# Aplicación de K-Means para encontrar el número óptimo de clusters
set.seed(123) # Establecer semilla para reproducibilidad

# Definir la función para calcular el total within sum of squares (WSS)
wss <- function(k) {
  #Ejecutar K-means y devolver el total withim sum of squares (WSS)
  kmeans(normalized_data[, -c(ncol(normalized_data)-1, ncol(normalized_data))], k, nstart = 25)$tot.withinss
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
kmeans_result <- kmeans(normalized_data[, -c(ncol(normalized_data)-1, ncol(normalized_data))], centers = optimal_k, nstart = 25)

# Añadir el resultado de los clusters al dataframe original
data$Cluster <- as.factor(kmeans_result$cluster)

# Visualizar los clusters utilizando ggplot2
cluster_plot <- ggplot(data, aes(x = normalized_data[,1], y = normalized_data[,2], color = Cluster)) + 
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
dist_matrix <- dist(normalized_data)

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




