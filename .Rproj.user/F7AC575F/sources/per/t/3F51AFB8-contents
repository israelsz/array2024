library("ape")
library("phangorn")
library("ggtree")       # For phylogenetic tree visualization
library("ggplot2")      # For plotting
library("gridExtra")    # For organizing plots in a grid
library("phytools")

datos_dna <- read.dna("C:\\Users\\Israel\\Documents\\Codes\\Taller de investigación\\array-code\\array-code-definitivo\\datasets\\primates_14.phylip")

###############################################################################

write.FASTA(datos_dna,"primates_14.fasta")
####################################################
# Convertir el objeto DNAbin a una matriz
matriz_dna <- as.matrix(datos_dna)

# Obtener los nombres de los taxones (nombres de las filas de la matriz)
nombres_taxones <- rownames(matriz_dna)


# Convertir la matriz DNAbin a una matriz de caracteres legibles
matriz_dna_caracter <- as.character(matriz_dna)

# Crear un data.frame con los nombres de los taxones y las secuencias legibles
df_dna_caracter <- data.frame(species = rownames(matriz_dna_caracter), matriz_dna_caracter, stringsAsFactors = FALSE)




# Crear un data.frame combinando los nombres de los taxones y las secuencias
df_dna <- data.frame(Taxon = nombres_taxones, matriz_dna, stringsAsFactors = FALSE)

# Verificar el resultado
head(df_dna)



############################################################################



datos_dna <- read.dna("C:\\Users\\Israel\\Documents\\Codes\\Taller de investigación\\array-code\\array-code-definitivo\\datasets\\primates_14.phylip")

# Convert the multiple sequence alignment result to a phyDat object for downstream analyses in phangorn
phyDat_msa_primate_sample = as.phyDat(datos_dna)

# Distance Calculation
# Calculate the Hamming distance matrix for the given aligned sequences
# This serves as a measure of pairwise sequence dissimilarity for tree construction
D_hamming = dist.hamming(phyDat_msa_primate_sample)

# Compute the Neighbor Joining tree
nj_tree = nj(D_hamming) #Aplica nj
nj_tree$edge.length[which(nj_tree$edge.length<0)]= 1 #Reemplaza los edges donde sea menor a cero por 1
nj_tree = midpoint(multi2di(nj_tree)) # Binariza los arboles


# Function to plot trees
plot_tree = function(tree_plot, title_plot, max_x) {
  g = ggtree(tree_plot, color = "#00A499", size = 1)
  
  # Customize the appearance of the plot and tip labels
  g = g + geom_tiplab(size = 4, color = "black", align = TRUE) +
    geom_nodepoint(size = 3, color = "#c7254e") +
    labs(title = title_plot, size = 6) +
    xlim(0, max_x) +
    theme(
      # Remove axis lines and text
      axis.line = element_blank(),
      axis.text = element_blank(),
      # Remove all grids
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Adjust margins and legend position
      plot.margin = margin(t = 0, r = -5, b = 0, l = 0, unit = "pt"),
      legend.position = 'top'
    )
  
  # Add labels for branch lengths with a smaller font size
  g = g + geom_text(aes(x = branch, label = round(branch, 2)), size = 3, color = "black", vjust = -0.5, hjust = 0)
  
  return(g)
}

# Create plots for the NJ tree
# Specify titles for each plot
g_NJ  = plot_tree(nj_tree, "Neighbor Joining (NJ)", 0.4)

# Organize the two plots in a one-row, two-column grid
grid.arrange(g_NJ, nrow = 1, ncol = 1)

arboles = rNNI(nj_tree, moves = 1)
g_nni  = plot_tree(arboles, "Arbol post NNI", 0.4)
grid.arrange(g_NJ,g_nni, nrow = 1, ncol = 2)

# Prueba mas de un arbol al mismo tiempo
multiarboles = rNNI(nj_tree, moves = 1, n = 100)


g1 = plot_tree(multiarboles[[1]], "arbol1", 0.4)
g2 = plot_tree(multiarboles[[2]], "arbol2", 0.4)
g3 = plot_tree(multiarboles[[3]], "arbol3", 0.4)
g4 = plot_tree(multiarboles[[4]], "arbol4", 0.4)
g5 = plot_tree(multiarboles[[5]], "arbol5", 0.4)
g6 = plot_tree(multiarboles[[6]], "arbol6", 0.4)
g7 = plot_tree(multiarboles[[7]], "arbol7", 0.4)
g8 = plot_tree(multiarboles[[8]], "arbol8", 0.4)
g9 = plot_tree(multiarboles[[9]], "arbol9", 0.4)
g10 = plot_tree(multiarboles[[10]], "arbol10", 0.4)

grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,g10,nrow = 2, ncol = 5)


# Saving to Newick
write.tree(multiarboles, file="multiarboles.nwk")


for (i in seq_along(multiarboles)) {
  # Verifica si el árbol en la posición i no es binario
  if (!is.binary(multiarboles[[i]])) {
    # Si no es binario, lo reemplaza con la versión binarizada
    multiarboles[[i]] = multi2di(multiarboles[[i]])
  }
}


# Arbol de consenso

# Compute the average tree by suppressing any console messages
average_tree = averageTree(multiarboles, method = "path.difference",quiet=T)

average_tree = midpoint(average_tree)
print(is.binary(average_tree))

#######################
# Arbol Medioide consenso
# Initialize an empty matrix for RF distances with dimensions 100x100
RF_distance_matrix = matrix(0, 10, 10)

# Calculate Robinson-Foulds (RF) distances between all pairs of trees 
for (a in 1:10) {
  for (b in 1:10) {
    # Calculate the RF distance between tree 'a' and tree 'b' and store it in the matrix
    RF_distance_matrix[a, b] = RF.dist(multiarboles[[a]], multiarboles[[b]], normalize = TRUE)
  }
}

# Find the index of the tree that minimizes the average RF distance
idx_med = which.min(colMeans(RF_distance_matrix))

# Medoid: Select the tree that corresponds to the medoid index
medoid_tree = multiarboles[[idx_med]]
#Reemplaza los edges donde sea menor a cero por 1
#medoid_tree$edge.length[which(medoid_tree$edge.length<=0)]= 1 
medoid_tree = midpoint(medoid_tree)
# Se verifica si es binario
print(is.binary(medoid_tree))


g_consenso  = plot_tree(medoid_tree, "Arbol Consenso", 0.4)
grid.arrange(g_consenso, nrow = 1, ncol = 1)

# Vamos a intentar cargar una red filo

red_filo = read.tree("./array-code-definitivo/phylonet.nwk")

red_filo$edge.length[which(red_filo$edge.length==NaN)] = 1 #Reemplaza los edges donde sea menor a cero por 1
red_filo$edge.length[is.na(red_filo$edge.length)] = 1

red_filo_binaria = midpoint(multi2di(red_filo))

# Graficar la red filogenética
plot(red_filo_binaria, type = "phylogram", edge.width = 2, no.margin = TRUE)

# Saving to Newick
write.tree(red_filo_binaria, file="red_filo_binaria.nwk")

library(igraph)

# Convertir la red filogenética a un objeto igraph (si es necesario)
red_igraph = as.igraph(red_filo)

# Graficar la red usando igraph
plot(red_igraph)


# Obtener las etiquetas de los nodos
etiquetas <- red_filo$tip.label

# Verificar si hay etiquetas duplicadas
duplicadas <- etiquetas[duplicated(etiquetas)]

# Imprimir las etiquetas duplicadas (si las hay)
if (length(duplicadas) > 0) {
  print(paste("Etiquetas duplicadas encontradas:", duplicadas))
} else {
  print("No hay etiquetas duplicadas")
}

###############
# Crear red filo directo en R
red_r = consensusNet(multiarboles, .8)

plot(red_r, type = "2D", show.edge.label=TRUE)

# Saving to Newick
write.tree(red_r, file="red_r.nwk")
red_igraph <- as.igraph(red_r)

# Graficar la red usando igraph
plot(red_igraph)