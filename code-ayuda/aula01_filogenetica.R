# Install the BiocManager package from CRAN
library("BiocManager") #install.packages("BiocManager")

# Use BiocManager to install the 'msa' package from Bioconductor
require("msa") #BiocManager::install("msa")

# Load the 'msa' package into the R session
library("msa")

# -----------------------------------------------------

# Load the 'Biostrings' package. 
# 'Biostrings' provides efficient containers to handle and analyze biological sequences.
library("Biostrings")

# Define a vector containing four DNA sequences
sequences = c("ACTGGCTG", "ACTGCTG", "AGTGACT", "TGTGACTGA")

# Convert the vector of sequences into a DNAStringSet object using the Biostrings package
biostrings_sequences = DNAStringSet(sequences)

# Perform multiple sequence alignment on the DNAStringSet object using the 'msa' package
msa_result_sample = msa(biostrings_sequences,method = "ClustalW")

# ---------------------------------------------------------------------------------------

# Define a vector containing amino acid sequences of varying lengths
amino_sequences =  c("ARNDCLTQ", "ARN", "ARMDCQK", "RRRTDGPSW")
amino_sequences = AAStringSet(amino_sequences) 
print(amino_sequences)

# Realiza la alineación múltiple de las secuencias de aminoácidos.
# El argumento 'type' se establece en "protein" para indicar que estamos trabajando con secuencias de proteínas.
amino_sequences = msa(amino_sequences, type = "protein")
print(amino_sequences)


# Carga la matriz de sustitución BLOSUM62, que es una matriz comúnmente utilizada para comparar secuencias de proteínas.
data(BLOSUM62)

# Calcula la puntuación de conservación para las secuencias de aminoácido
amino_sequences = msaConservationScore(amino_sequences, BLOSUM62)
print(amino_sequences)

# Now we need to open the salmon gene file
require(seqinr)
salmon_gene = read.fasta("ALKBH5.fasta", as.string = 'T')

# Convert the vector of sequences into a DNAStringSet object using the Biostrings package
salmon_sequences = DNAStringSet(unlist(salmon_gene))
# Perform multiple sequence alignment on the DNAStringSet object using the 'msa' package
msa_salmon_result_sample = msa(salmon_sequences,method = "ClustalW")
print(msa_salmon_result_sample)

# Parwise sequence distance

# Load the phangorn library: Provides tools for phylogenetic reconstruction and analysis
library("phangorn")
# Load the ape library: Essential for reading, writing, and manipulating phylogenetic trees
library("ape")
# Load the knitr library: Used for generating HTML tables
library("knitr")

# Define a vector containing four DNA sequences
# Convert the multiple sequence alignment result to a phyDat object for downstream analyses in phangorn
phyDat_msa_sample = as.phyDat(msa_result_sample)

# Assign custom names to the sequences
names(phyDat_msa_sample) = c("Seq 1", "Seq 2", "Seq 3", "Seq 4")

# Compute the Hamming distance matrix for the aligned sequences
D_hamming = dist.hamming(phyDat_msa_sample, ratio = FALSE)

# Round the values in the distance matrix to two decimal places
D_hamming = round(as.matrix(D_hamming), 2)
print(D_hamming)

# Print the Hamming distance matrix as an HTML table with a 20% width
# kable(D_hamming, format = "html", table.attr = "style='width:20%;'")

# Distance Position

# Compute the distance matrix using the 'dist.hamming' function (assuming it's from a specific package)
D_p = dist.hamming(phyDat_msa_sample,ratio = T)
# Convert the distance matrix to a matrix and round the values
D_p = round(as.matrix(D_p),2)
print(D_p)


# Edition Distance

# Load the "stringdist" library, which is used to calculate distances between character strings.
library("stringdist")

# Convert the sequences into character strings (if they are not already).
seq_chars = as.character(sequences)

# Calculate the Levenshtein distance matrix (edit distance) between the character sequences.
D_Edit = stringdistmatrix(seq_chars, seq_chars, method = "lv")

# Assign custom names to the sequences
colnames(D_Edit) = c("Seq 1", "Seq 2", "Seq 3", "Seq 4")
rownames(D_Edit) = c("Seq 1", "Seq 2", "Seq 3", "Seq 4")

# Round the values in the distance matrix to two decimal places.
D_Edit = round(D_Edit, 2)
print(D_Edit)

# Jukes-Cantor

# Compute the p distance matrix for the aligned sequences
D_JK69 = dist.dna(as.DNAbin(phyDat_msa_sample), model = "JC69")
# Convert the distance matrix to a matrix and round the values
D_JK69 = round(as.matrix(D_JK69),2)
print(D_JK69)

# Kimura 2 parameter

# Calculate the distance matrix with the K80 model and store it in 'D_K80'
D_K80 = dist.dna(as.DNAbin(phyDat_msa_sample), model = "K80")
# Convert the distance matrix to a matrix and round the values to two decimal places
D_K80 = round(as.matrix(D_K80), 2)

# Replace any NA (Not Available) values in the matrix with 0
D_K80[which(is.na(D_K80))] = 0
print(D_K80)


# Now we going to generate tree

# Set a seed for reproducibility of random processes
set.seed(1)
# Generate a random tree with 5 tips using the 'rtree' function
random_tree = rtree(5)
# Plot the generated random tree
plot(random_tree)

# Get the tip labels from the 'random_tree' object and store them in 'labels'
labels = random_tree$tip.label
# Print the tip labels
print(labels)

# Assign new tip labels to 'random_tree'
random_tree$tip.label = c("Species 1", "Species 2", "Species 3", "Species 4", "Species 5")
# Plot the tree with the updated tip labels
plot(random_tree)


# Adding / remove branches

# Remove a tip labeled "Species 3" from the 'random_tree' object
random_tree1 = drop.tip(random_tree, "Species 3")

# Add a new tip labeled "Species 3" connected to "Species 1" with an edge length of 1
random_tree2 = add.tips(random_tree1, "Species 3", "Species 1", edge.length = 1)

# Plot the trees side by side
# Set the plot layout to display two plots in a single row
par(mfrow = c(1, 2))
plot(random_tree1)
plot(random_tree2)


# Define all edge lengths in 'random_tree' equal to 1
random_tree2$edge.length = random_tree2$edge.length*0 + 1
# Plot the tree with the modified edge lengths
plot(random_tree2)


# Root the 'random_tree' at a specified node labeled as "Species 1" (Unroot() can be applied to remove the root)
rooted_tree1 = root(random_tree, "Species 1")
# Root the 'random_tree' at its midpoint
rooted_tree2 = midpoint(random_tree)
# Plot the trees
par(mfrow = c(1, 2))
plot(rooted_tree1)
plot(rooted_tree2)


# Create a multifurcating tree
# Remove tip labeled "Species 3" from 'random_tree' and then add it back as a multifurcation with 7 descendants, each with an edge length of 1
random_tree1 = drop.tip(random_tree, "Species 3")
random_tree1 = add.tips(random_tree1, "Species 3", 7, edge.length = 1)

# Solve polytomy by converting the multifurcating tree to a bifurcating tree
random_tree2 = multi2di(random_tree1)

# Plot the trees side by side
# Set the plot layout to display two plots in a single row
par(mfrow = c(1, 2))

# Plot the trees
plot(random_tree1)
plot(random_tree2)


# Load the necessary libraries
library("ggtree")       # For phylogenetic tree visualization
library("ggplot2")      # For plotting
library("gridExtra")    # For organizing plots in a grid

# Create the first dendrogram-style plot
g1 = ggtree(rooted_tree2, layout = "dendrogram", ladderize = TRUE, color = "#00A499")

# Customize the appearance of the first plot and tip labels
g1 = g1 + geom_tiplab(size = 4, color = "black", hjust = 1, angle = 90)
g1 = g1 + theme(plot.margin = margin(t = 0, r = 0, b = 50, l = 0, unit = "pt"), legend.position = "top")
g1 = g1 + geom_nodepoint(size = 3, color = "#c7254e") + theme(legend.position = 'top')
g1 = g1 + labs(title = "Midpoint rooted tree (dendrogram)")

# Create the second circular-style plot
g2 = ggtree(rooted_tree2, layout = "circular", ladderize = TRUE, color = "#00A499")

# Customize the appearance of the second plot and tip labels
g2 = g2 + geom_tiplab(size = 4, color = "black", hjust = 1, angle = 90)
g2 = g2 + theme(plot.margin = margin(t = 0, r = 0, b = 50, l = 0, unit = "pt"), legend.position = "top")
g2 = g2 + geom_nodepoint(size = 3, color = "#c7254e") + theme(legend.position = 'top')
g2 = g2 + labs(title = "Midpoint rooted tree (circular)")

# Organize the two plots in a two-column grid
grid.arrange(g1, g2, ncol = 2)

# READ/WRITE Trees

# Saving to Newick
write.tree(rooted_tree2, file="output_newick_file.nwk")
# Saving to Nexus (using 'ape' package)
write.nexus(rooted_tree2, file="output_nexus_file.nex")

# For Newick
tree_newick = read.tree("output_newick_file.nwk")
# For Nexus
tree_nexus = read.nexus("output_nexus_file.nex")

plot(tree_newick)
plot(tree_nexus)

print(is.rooted(tree_newick))
print(is.binary(tree_newick))
print(is.rooted(tree_nexus))
print(is.binary(tree_nexus))

# Defining a multiPhylo structure

list_tree = c(tree_newick, tree_nexus, random_tree1, random_tree2)

print(list_tree)

for(x in 1:4){
  plot(list_tree[x])
}

# Now we going to make the first inference of filogenetics tree

# Define a vector containing four DNA sequences
# Convert the multiple sequence alignment result to a phyDat object for downstream analyses in phangorn
phyDat_msa_salmon_sample = as.phyDat(msa_salmon_result_sample)


# Distance Calculation
# Calculate the Hamming distance matrix for the given aligned sequences
# This serves as a measure of pairwise sequence dissimilarity for tree construction
D_hamming = dist.hamming(phyDat_msa_salmon_sample)

# Nearest Neighbor Clustering (NNC)
# Construct a phylogenetic tree using the Nearest Neighbor Clustering method
# This method groups sequences based on the nearest (smallest) pairwise distance
tree_NNC = upgma(D_hamming, "single")
tree_NNC = midpoint(tree_NNC)

# Furthest Neighbor (FN)
# Construct a phylogenetic tree using the Furthest Neighbor method
# This method groups sequences based on the furthest (largest) pairwise distance
tree_FN = upgma(D_hamming, "complete")
tree_FN = midpoint(tree_FN)

# Weighted Pair Group Method with Arithmetic Mean (WPGMA)
# Construct a phylogenetic tree using the WPGMA method
# This method considers all pairwise distances for clustering and calculates average distances
tree_WPGMA = wpgma(D_hamming)
tree_WPGMA = midpoint(tree_WPGMA)

# Unweighted Pair-Group Centroid Method (UPGMC)
# Construct a phylogenetic tree using the UPGMC method
# This method clusters sequences based on the centroid distance without considering the number of sequences in each cluster
tree_UPGMC = upgma(D_hamming, "centroid")
tree_UPGMC = midpoint(tree_UPGMC)

# Weighted Pair-Group Centroid Method (WPGGMC)
# Construct a phylogenetic tree using the WPGGMC method
# This method clusters sequences based on the centroid distance and considers the number of sequences in each cluster
tree_WPGGMC = wpgma(D_hamming, "centroid")
tree_WPGGMC = midpoint(tree_WPGGMC)

# Unweighted Pair Group Method with Arithmetic Mean (UPGMA)
# Construct a phylogenetic tree using the UPGMA method
# This method is similar to WPGMA but does not consider the number of sequences in each cluster for centroid calculation
tree_UPGMA = upgma(D_hamming)
tree_UPGMA = midpoint(tree_UPGMA)

# Ward's Minimum Variance Method
# Construct a phylogenetic tree using Ward's method
# This method minimizes the total within-cluster variance
tree_WARD = upgma(D_hamming, "ward")
tree_WARD = midpoint(tree_WARD)

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

# Generate plots for each tree
g_NNC    = plot_tree(tree_NNC, "Nearest Neighbor Clustering (NNC)", 0.3)
g_FN     = plot_tree(tree_FN, "Furthest Neighbor (FN)", 0.3)
g_WPGMA  = plot_tree(tree_WPGMA, "Weighted Pair Group Method with Arithmetic Mean (WPGMA)", 0.3)
g_UPGMC  = plot_tree(tree_UPGMC, "Unweighted Pair-Group Centroid Method (UPGMC)", 0.3)
g_WPGGMC = plot_tree(tree_WPGGMC, "Weighted Pair-Group Centroid Method (WPGGMC)", 0.3)
g_UPGMA  = plot_tree(tree_UPGMA, "Unweighted Pair Group Method with Arithmetic Mean (UPGMA)", 0.3)
g_ward   = plot_tree(tree_WARD, "Ward's Minimum Variance Method", 0.3)

# Arrange the plots in a grid
grid.arrange(g_NNC, g_FN, g_WPGMA, g_UPGMC, g_WPGGMC, g_UPGMA, g_ward, nrow = 4, ncol = 2)


# Challenge 04
# Construct UPGMA and WPGMA

primates_data = read.dna("primates_14.phylip")

# Convert the multiple sequence alignment result to a phyDat object for downstream analyses in phangorn
phyDat_msa_primate_sample = as.phyDat(primates_data)


# Distance Calculation
# Calculate the Hamming distance matrix for the given aligned sequences
# This serves as a measure of pairwise sequence dissimilarity for tree construction
D_hamming = dist.hamming(phyDat_msa_primate_sample)

# Unweighted Pair Group Method with Arithmetic Mean (UPGMA)
# Construct a phylogenetic tree using the UPGMA method
# This method is similar to WPGMA but does not consider the number of sequences in each cluster for centroid calculation
tree_UPGMA = upgma(D_hamming)
tree_UPGMA = midpoint(tree_UPGMA)

# Weighted Pair Group Method with Arithmetic Mean (WPGMA)
# Construct a phylogenetic tree using the WPGMA method
# This method considers all pairwise distances for clustering and calculates average distances
tree_WPGMA = wpgma(D_hamming)
tree_WPGMA = midpoint(tree_WPGMA)

g_WPGMA  = plot_tree(tree_WPGMA, "Weighted Pair Group Method with Arithmetic Mean (WPGMA)", 0.3)
g_WPGMA  = plot_tree(tree_WPGMA, "Weighted Pair Group Method with Arithmetic Mean (WPGMA)", 0.3)


# Arrange the plots in a grid
grid.arrange(g_WPGMA, g_WPGMA, ncol = 2)

# Additive methods
# This trees was constructed on primates_14 data

# Compute the Neighbor Joining tree
nj_tree = nj(D_hamming)
nj_tree$edge.length[which(nj_tree$edge.length<0)]=0
nj_tree = midpoint(multi2di(nj_tree))

# Compute the Neighbor Joining tree
bionj_tree = bionj(D_hamming)
bionj_tree$edge.length[which(bionj_tree$edge.length<0)]=0
bionj_tree = midpoint(multi2di(bionj_tree))

# Create plots for the UPGMA and WPGMA trees
# Specify titles for each plot
g_NJ    = plot_tree(nj_tree,    "Neighbor Joining (NJ)",0.3)
g_BioNJ = plot_tree(bionj_tree, "Bio-Neighbor Joining (BioNJ)",0.3)

# Organize the two plots in a one-row, two-column grid
grid.arrange(g_NJ, g_BioNJ, nrow = 1, ncol = 2)



# Compute the Minimum Evolution tree
me_tree = fastme.bal(D_hamming)
me_tree = multi2di(me_tree)

# Compute the Minimum Evolution tree
ols_tree = fastme.ols(D_hamming)
ols_tree = multi2di(ols_tree)

# Organize the two plots in a one-row, two-column grid
g_me_tree  = plot_tree(me_tree,  "Balanced ME",0.55)
g_ols_tree = plot_tree(ols_tree, "Ordinary Least Squares",0.55)
grid.arrange(g_me_tree, g_ols_tree, nrow = 1, ncol = 2)


# The function takes a phylogenetic tree as input and returns the ME score
ME_score = function(tree_input) {
  # The ME score is the sum of the lengths of all the edges in the tree
  # Extract the edge lengths from the input tree and calculate the ME score
  res_score = round(sum(tree_input$edge.length),2)
  
  # Return the calculated ME score
  return(res_score)
}

# Calculate and print the ME scores for the given trees (me_tree and ols_tree)
print(paste("Minimum Evolution Score for me_tree:", ME_score(me_tree)))

print(paste("Minimum Evolution Score for ols_tree:", ME_score(ols_tree)))



# Define a function to calculate the Minimum Evolution (ME) score using Least Squares (LS)
# The function takes a tree and a distance matrix as input and returns the LS score
LS_score = function(tree_input, distance_input) {
  # Calculate the cophenetic matrix from the input tree
  # The cophenetic matrix represents the patristic (tree-derived) distances between each pair of taxa
  D = cophenetic(tree_input)
  
  # Extract the input distance matrix
  d = distance_input
  
  # Calculate the LS score by summing the squared differences between the observed distances (d)
  # and the cophenetic distances (D)
  res_score = round(sum((D - d)^2),2)
  
  # Return the calculated LS score
  return(res_score)
}

# Calculate and print the ME scores for the given trees (me_tree and ols_tree)
# using the provided Hamming distance matrix (D_hamming)
print(paste("Minimum Evolution Score for me_tree:", LS_score(me_tree, D_hamming)))

print(paste("Minimum Evolution Score for ols_tree:", LS_score(ols_tree, D_hamming)))



# Parsimony Tree

suppressWarnings({
  suppressMessages({
    options(verbose = FALSE)
    
    # Distance Calculation
    # Compute the Hamming distance matrix for the given aligned sequences
    # This will serve as a measure of pairwise sequence dissimilarity
    D_hamming = dist.hamming(phyDat_msa_primate_sample)
    
    # Base Tree Construction
    # Compute the Neighbor Joining (NJ) tree as a base tree
    # NJ is a distance-based method for constructing phylogenetic trees
    nj_tree = nj(D_hamming)
    
    # Correct Negative Branch Lengths
    # Negative branch lengths are not biologically meaningful,
    # hence set any negative branch lengths to 0
    nj_tree$edge.length[which(nj_tree$edge.length < 0)] = 0
    
    # Rooting the Tree
    # Root the NJ tree at its midpoint to ensure ultrametricity,
    # and convert any multifurcating nodes to bifurcating nodes
    nj_tree = midpoint(multi2di(nj_tree))
    
    # Parsimony Optimization
    # Optimize the parsimony of the NJ tree
    # The optim.parsimony function refines the tree topology to minimize the parsimony score
    par_tree = optim.parsimony(nj_tree, phyDat_msa_primate_sample)
    
    
    #Edge estimation
    par_tree=acctran(par_tree, phyDat_msa_primate_sample)
    par_tree=midpoint(par_tree)
    
    # Calculating and Printing Parsimony
    # Calculate and print the parsimony score of the optimized tree
    # Lower parsimony scores indicate more parsimonious (simpler) trees
    print(paste("Parsimony Score for the optimized tree (NJ):", parsimony(par_tree, phyDat_msa_primate_sample)))
    
    # Parsimony For BioNJ
    
    # Base Tree Construction
    # Compute the Neighbor Joining (NJ) tree as a base tree
    # NJ is a distance-based method for constructing phylogenetic trees
    bionj_tree = bionj(D_hamming)
    
    # Correct Negative Branch Lengths
    # Negative branch lengths are not biologically meaningful,
    # hence set any negative branch lengths to 0
    bionj_tree$edge.length[which(bionj_tree$edge.length < 0)] = 0
    
    # Rooting the Tree
    # Root the NJ tree at its midpoint to ensure ultrametricity,
    # and convert any multifurcating nodes to bifurcating nodes
    bionj_tree = midpoint(multi2di(bionj_tree))
    
    # Parsimony Optimization
    # Optimize the parsimony of the NJ tree
    # The optim.parsimony function refines the tree topology to minimize the parsimony score
    par_tree = optim.parsimony(bionj_tree, phyDat_msa_primate_sample)
    
    
    #Edge estimation
    par_tree=acctran(par_tree, phyDat_msa_primate_sample)
    par_tree=midpoint(par_tree)
    
    # Calculating and Printing Parsimony
    # Calculate and print the parsimony score of the optimized tree
    # Lower parsimony scores indicate more parsimonious (simpler) trees
    print(paste("Parsimony Score for the optimized tree (BioNJ):", parsimony(par_tree, phyDat_msa_primate_sample)))
  })
})

g_par_tree = plot_tree(par_tree, "Parsimony tree",150)
grid.arrange(g_par_tree, nrow = 1, ncol = 1)

print(par_tree[1])



# Likehood method

# Suppress warnings and messages
suppressWarnings({
  suppressMessages({
    options(verbose = FALSE)
    
    # Distance Calculation
    # Compute the Hamming distance matrix for the given aligned sequences
    # This will serve as a measure of pairwise sequence dissimilarity
    D_hamming = dist.hamming(phyDat_msa_primate_sample)
    
    # Base Tree Construction
    # Compute the Neighbor Joining (NJ) tree as a base tree
    # NJ is a distance-based method for constructing phylogenetic trees
    nj_tree = nj(D_hamming)
    
    # Correct Negative Branch Lengths
    # Negative branch lengths are not biologically meaningful,
    # hence set any negative branch lengths to 0
    nj_tree$edge.length[which(nj_tree$edge.length < 0)] = 0
    
    # Model Testing
    # Conduct a model test to determine the best-fitting evolutionary model based on AIC
    mt = modelTest(phyDat_msa_primate_sample, nj_tree, model = c("JC", "K80", "HKY", "GTR"))
    var_mt = attr(mt, "env")
    evolutionary_model = eval(get(mt$Model[which.min(mt$AIC)], var_mt), var_mt)
    
    # Initialize PML Tree
    # Create an initial PML (Phylogenetic Maximum Likelihood) tree using the NJ tree and selected evolutionary model
    initial_pml_tree = pml(tree = nj_tree, data = phyDat_msa_primate_sample, model = evolutionary_model$model,
                           bf = evolutionary_model$bf, Q = evolutionary_model$Q, inv = evolutionary_model$inv,
                           k = evolutionary_model$k, shape = evolutionary_model$shape)
    
    # Optimize PML Tree
    # Optimize the PML tree through various parameters and tree rearrangement methods
    pml_tree = optim.pml(initial_pml_tree, model = evolutionary_model$model, optInv = T, optGamma = F,
                         rearrangement = "stochastic", optBf = F, optQ = F,
                         optEdge = F, optNni = TRUE, control = pml.control(trace = 0))
  })
})

# Visualize the Optimized PML Tree
# Create a visualization of the optimized PML tree with a specified title
pml_tree$tree = multi2di(pml_tree$tree)
g_pml_tree = plot_tree(pml_tree$tree, "Likelihood tree", 13)
grid.arrange(g_pml_tree, nrow = 1, ncol = 1)

print(pml_tree)

install.packages("kmer")

library("kmer")
kmer_distance = round(kdistance(phyDat_msa_primate_sample, k = 5),2)
print(kmer_distance)

# Compute the Neighbor Joining tree
nj_km_tree = nj(kmer_distance)
nj_km_tree$edge.length[which(nj_km_tree$edge.length<0)]=0
nj_km_tree = midpoint(multi2di(nj_km_tree))

# Visualize the Optimized PML Tree
# Create a visualization of the optimized PML tree with a specified title
g_nj_km_tree = plot_tree(nj_km_tree, "NJ kmer tree", 0.5)
grid.arrange(g_nj_km_tree, nrow = 1, ncol = 1)