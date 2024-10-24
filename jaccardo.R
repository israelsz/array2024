library("cluster") #Instalat
library("clusteval") #Instalar


grupo_1 <- c(1, 1, 2, 2, 1, 3, 3, 2, 1, 2) #10 objetos que fueron asignados a 3 grupos
grupo_2 <- c(1, 2, 2, 3, 1, 3, 2, 2, 3, 3) #10 objetos que fueron asignados a 3 grupos

indice_jaccard <- cluster_similarity(grupo_1,grupo_2, similarity = "jaccard")
indice_jaccard 

install.packages("phylotools",dependencies = TRUE)
install.packages("phangorn",dependencies = TRUE)
install.packages("ggtree",dependencies = TRUE)
install.packages("phytools",dependencies = TRUE)
install.packages("igraph",dependencies = TRUE)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")

install.packages("emoa",dependencies = TRUE)