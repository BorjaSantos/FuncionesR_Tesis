# FUNCION PARA OBTENER LA MATRIZ DE CORRELACIONES MEDIANTE LA UTILIZACION DE LAS REGLAS DE RUBIN.

# INICIO: 5 / SEPTIEMBRE / 2014
# FIN: 8 / SEPTIEMBRE / 2014

# INPUT: 
# lista.matrices: Muestra de posibles matrices de correlacion.
# n: numero de variables.

# OUTPUT:
# matriz imputada: Matriz de correlaciones para el metaanalisis obtenida mediante las reglas de Rubin.
# sd.matrix: Desviacion estandar asociadada a la imputacion multiple.

# QUE HACE LA FUNCION:
# A partir del numero de variables y una lista con posibles matrices de correlacion, estima mediante las reglas de Rubin, 
# la matriz de correlacion a emplear en el metaanalisis multivariado y la matriz con las desviaciones estandar asociadas.

matriz.imputada <- function(lista.matrices,n){
  matrix.aux <- matrix(0,nrow=length(lista.matrices),ncol=n*n)
  mean.matrix <- matrix(0,nrow=n,ncol=n)
  sd.matrix <- matrix(0,nrow=n,ncol=n)
  sd.vector <- vector("numeric",length=n*n)
  for(i in 1:length(lista.matrices)){
    matrix.aux[i,] <- as.vector(lista.matrices[[i]])
  }
  vector.imputaciones <- colMeans(matrix.aux)
  mean.matrix <- matrix(vector.imputaciones,ncol=n,nrow=n,byrow=TRUE)
  for (j in 1:dim(matrix.aux)[2]){
    sd.vector[j] <- sd(matrix.aux[,j])
  }
  sd.matrix <- matrix(sd.vector,ncol=n,nrow=n,byrow=TRUE)
  output <- list(matriz.imputada = mean.matrix, sd.imputacion = sd.matrix)
  return(output)
}