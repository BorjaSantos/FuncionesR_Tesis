# FUNCION PARA GENERAR LA LISTA DE MATRICES DE VARIANZAS-COVARIANZAS A EMPLEAR EN EL METAANALISIS
# MULTIVARIADO

# INCIO: 8 / SEPTIEMBRE / 2014
# FIN: 8 / SEPTIEMBRE / 2014

# INPUTS:
# matriz.correlacion: matriz de correlacion obtenida mediante imputaciones multiples y la regla de Rubin.
# matrices.sd: matrices diagonales con las desviaciones estandar de los estimadores del efecto de cada estudio.

# OUTPUT:
# matrices.covarianzas: Lista con las pmatrices de covarianzas para cada uno de los estudios incluidos en el metaanalisis.

matrices.varcovar <- function(matriz.correlacion,matrices.sd){
  matrices.covarianzas <- vector("list",length(matrices.sd))
  for( i in 1:length(matrices.sd)){
    matrices.covarianzas[[i]] <- matrices.sd[[i]]%*%matriz.correlacion%*%matrices.sd[[i]] 
  }
  return(matrices.covarianzas)
}