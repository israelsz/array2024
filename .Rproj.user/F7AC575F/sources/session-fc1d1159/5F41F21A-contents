library("nsga2R") 
library("emoa")

set.seed(1234)
# randomly generate a polulation of fifty chromosomes, each with two objectives
y = matrix(runif(20, -5, 5), nrow=10, ncol=2)
rankIdxList <- fastNonDominatedSorting(y)
rankIdxList

#######################################################

library(mco)
popSize <- 50
lowerBounds <- rep(0,30)
upperBounds <- rep(1,30)
varNo <- length(lowerBounds)
objDim <- 2
set.seed(1234)
population <- t(sapply(1:popSize, function(u) array(runif(length(lowerBounds),
                                                          lowerBounds,upperBounds))))
population <- cbind(population, t(apply(population,1,zdt2)))
ranking <- fastNonDominatedSorting(population[,(varNo+1):(varNo+objDim)])
rnkIndex <- integer(popSize)
i <- 1
while (i <= length(ranking)) {
  rnkIndex[ranking[[i]]] <- i
  i <- i + 1
}
population <- cbind(population,rnkIndex)
objRange <- apply(population[,(varNo+1):(varNo+objDim)], 2, max) -
  apply(population[,(varNo+1):(varNo+objDim)], 2, min)
cd <- crowdingDist4frnt(population,ranking,objRange)
cd

###########################################################
popSize <- 50
lowerBounds <- rep(0,30)
upperBounds <- rep(1,30)
varNo <- length(lowerBounds)
objDim <- 2
set.seed(1234)

# Generamos una población aleatoria de 50 soluciones con 30 variables de decisión
population <- t(sapply(1:popSize, function(u) array(runif(length(lowerBounds), lowerBounds, upperBounds))))

# Añadimos las funciones objetivo ZDT2 para cada solución en la población
population <- cbind(population, t(apply(population, 1, zdt2)))

# Clasificamos las soluciones en frentes de no dominación usando fastNonDominatedSorting
ranking <- fastNonDominatedSorting(population[,(varNo+1):(varNo+objDim)])

# Creamos un índice de los frentes
rnkIndex <- integer(popSize)
i <- 1
while (i <= length(ranking)) {
  rnkIndex[ranking[[i]]] <- i
  i <- i + 1
}

# Agregamos el rango de no dominación a la población
population <- cbind(population, rnkIndex)

# Calculamos el rango de cada función objetivo (máximo - mínimo)
objRange <- apply(population[,(varNo+1):(varNo+objDim)], 2, max) -
  apply(population[,(varNo+1):(varNo+objDim)], 2, min)

# Calculamos las distancias de hacinamiento usando crowdingDist4frnt
cd <- crowdingDist4frnt(population, ranking, objRange)

# Mostramos las distancias de hacinamiento
cd

#############################################
# Supongamos que ya tienes los puntajes de parsimonia y verosimilitud de cada red
# Aquí solo simulo estos valores, pero tú ya los tendrías calculados

# Criterios: Parsimonia y Verosimilitud (como las funciones objetivo)
parsimonia <- c(50, 48, 55, 53, 47, 60, 52, 54, 49, 51)
verosimilitud <- c(-300, -305, -290, -295, -310, -280, -302, -298, -307, -299)

# Crear la población que contiene estos dos criterios para cada red
# Aquí cada fila representa una red y las dos columnas son sus puntajes
population <- cbind(parsimonia, verosimilitud)

# Paso 1: Aplicar fastNonDominatedSorting para clasificar las redes en frentes
ranking <- fastNonDominatedSorting(population)

# Paso 2: Calcular el rango de cada función objetivo (máximo - mínimo)
# Este rango es necesario para calcular la distancia de hacinamiento
objRange <- apply(population, 2, max) - apply(population, 2, min)

# Paso 3: Aplicar crowdingDist4frnt para calcular la distancia de hacinamiento
crowdingDistances <- crowdingDist4frnt(population, ranking, objRange)

# Mostrar las distancias de hacinamiento para cada red
print("Distancias de hacinamiento por red:")
print(crowdingDistances)

# Si necesitas recortar las soluciones (basado en la diversidad o algún otro criterio),
# podrías usar estas distancias de hacinamiento para tomar decisiones.

###############################################################################
###############################################################################
# Cargar la librería emoa
library(emoa)

# Simular 100 valores de parsimonia y verosimilitud
set.seed(123)
parsimonia <- runif(100, min=40, max=70)  # 100 valores de parsimonia aleatorios
verosimilitud <- runif(100, min=-350, max=-250)  # 100 valores de verosimilitud aleatorios

# Crear la población que contiene estos dos criterios para cada red
population <- cbind(parsimonia, verosimilitud)

# Paso 1: Aplicar fastNonDominatedSorting para clasificar las redes en frentes
ranking <- fastNonDominatedSorting(population)

# Paso 2: Función personalizada para calcular la distancia de hacinamiento
calculate_crowding_distance <- function(pop, ranking) {
  # Inicializar un vector para almacenar las distancias
  distances <- rep(0, nrow(pop))
  
  # Iterar sobre cada frente en ranking
  for (i in 1:length(ranking)) {
    front <- ranking[[i]]  # Obtener las soluciones del frente actual
    
    if (length(front) < 2) {
      # Si el frente tiene menos de 2 soluciones, asignar Inf
      distances[front] <- Inf
    } else {
      # Calcular la distancia de hacinamiento con la función crowding_distance de emoa
      front_values <- pop[front, ]  # Valores de las soluciones en el frente
      distances[front] <- crowding_distance(front_values)  # Aplicar crowding_distance
    }
  }
  
  return(distances)
}

# Paso 3: Calcular las distancias de hacinamiento
crowdingDistances <- calculate_crowding_distance(population, ranking)

# Paso 4: Agregar las distancias de hacinamiento a la población y ordenarla
population_with_dist <- data.frame(population, crowding_dist = crowdingDistances)
population_with_dist$ranking <- unlist(lapply(1:length(ranking), function(i) rep(i, length(ranking[[i]]))))

# Ordenar por ranking y por crowding_dist (de mayor a menor)
population_with_dist <- population_with_dist[order(population_with_dist$ranking, -population_with_dist$crowding_dist),]

# Mostrar las primeras 10 soluciones
head(population_with_dist, 10)

###############################################################################
# Cargar la librería emoa
library(emoa)

# Implementación ajustada de crow_distance para parsimonia y verosimilitud
crow_distance <- function(table_fitness) {
  
  # Paso 1: Extraer las funciones objetivo (parsimonia y verosimilitud)
  fitness_values <- as.matrix(table_fitness[, 1:2])
  
  # Paso 2: Calcular la distancia de hacinamiento con crowding_distance
  # Transponemos la matriz de funciones objetivo y pasamos a crowding_distance
  distances <- crowding_distance(t(fitness_values))
  
  # Paso 3: Agregar las distancias de hacinamiento a la tabla de fitness
  table_fitness <- data.frame(table_fitness, crow_dist = distances)
  
  # Paso 4: Ordenar la tabla por ranking y luego por distancia de hacinamiento en orden descendente
  # Asumimos que la columna "ranking" ya está calculada y presente en table_fitness
  table_fitness <- table_fitness[with(table_fitness, order(ranking, -crow_dist)),]
  
  # Devolver la tabla con las distancias de hacinamiento incluidas
  return(table_fitness)
}

# Simular 100 valores de parsimonia y verosimilitud
set.seed(123)
parsimonia <- runif(100, min=40, max=70)  # 100 valores de parsimonia aleatorios
verosimilitud <- runif(100, min=-350, max=-250)  # 100 valores de verosimilitud aleatorios
ranking <- sample(1:5, 100, replace = TRUE)  # Asumimos que ya tienes un ranking de frentes

# Crear la tabla de fitness que contiene parsimonia, verosimilitud y ranking
table_fitness <- data.frame(parsimonia, verosimilitud, ranking)

# Aplicar la función para calcular la distancia de hacinamiento
result <- crow_distance(table_fitness)

# Mostrar las primeras filas del resultado
head(result, 10)

###############################################################################
###########################DEFINITIVO GOGOGOGO###################################################
# Cargar la librería emoa
library(emoa)

# Implementación ajustada de crow_distance para parsimonia y verosimilitud
crow_distance <- function(table_fitness) {
  
  # Paso 1: Extraer las funciones objetivo (parsimonia y verosimilitud)
  fitness_values <- as.matrix(table_fitness[, 1:2])
  
  # Paso 2: Calcular la distancia de hacinamiento con crowding_distance
  # Transponemos la matriz de funciones objetivo y pasamos a crowding_distance
  distances <- crowding_distance(t(fitness_values))
  
  # Paso 3: Agregar las distancias de hacinamiento a la tabla de fitness
  table_fitness <- data.frame(table_fitness, crow_dist = distances)
  
  # Paso 4: Ordenar la tabla por ranking y luego por distancia de hacinamiento en orden descendente
  # Asumimos que la columna "ranking" ya está calculada y presente en table_fitness
  table_fitness <- table_fitness[with(table_fitness, order(ranking, -crow_dist)),]
  
  # Devolver la tabla con las distancias de hacinamiento incluidas
  return(table_fitness)
}

# Simular 100 valores de parsimonia y verosimilitud
set.seed(123)
parsimonia <- runif(100, min=40, max=70)  # 100 valores de parsimonia aleatorios
verosimilitud <- runif(100, min=-350, max=-250)  # 100 valores de verosimilitud aleatorios

# Crear la tabla de fitness que contiene parsimonia y verosimilitud
population <- data.frame(parsimonia, verosimilitud)

# Calcular el ranking usando fastNonDominatedSorting
ranking <- fastNonDominatedSorting(as.matrix(population))

# Crear una columna de ranking para cada solución basada en fastNonDominatedSorting
# Las soluciones en el mismo frente reciben el mismo ranking
ranking_vector <- integer(nrow(population))
for (i in seq_along(ranking)) {
  ranking_vector[ranking[[i]]] <- i
}

# Agregar el ranking a la tabla de fitness
population$ranking <- ranking_vector

# Aplicar la función para calcular la distancia de hacinamiento
result <- crow_distance(population)

# Mostrar las primeras filas del resultado
head(result, 10)


