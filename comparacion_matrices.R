######################################################
######################################################
#####        COMPARACION ENTRE MATRICES         ######
######################################################
######################################################

# FECHA INICIO: 18 / NOVIEMBRE / 2014.
# FECHA FIN: 18 / NOVIEMBRE / 2014.

# ARGUMENTOS DE ENTRADA (INPUTS):
# i) lista.matrices1: Lista de matrices candidatas a ser matrices de correlacion.
# ii) lista.matrices2: Transformación de las matrices de la lista anterior a definidas positivas.

# ARGUMENTOS DE SALIDA (OUTPUT):
# i) Norma de la diferencia entre las matrices corresponidentes de ambas listas, para saber si estan proximas o no.

# QÚE HACE LA FUNCION:
# Paso 1: Calcula la diferencia entre las matrices correspondientes.
# Paso 2: Calcula la norma de Fröbenius de la matriz diferencia.

comparacion.matrices <- function(lista.matricesA,lista.matricesB){
  # Control de los argumentos de la funcion
  if (length(lista.matricesA)!=length(lista.matricesB)){
    stop("El tamaño de las listas de matrices a comparar no es el mismo")
  }
  # Lista donde se van a almacenar las diferencias entre las matrices de ambas listas
  dif.matrix.list <- vector("list",length(lista.matricesA))
  # Vector donde se va a almacenar la norma de la matriz diferencia
  frobenius.norma <- vector("numeric",length(lista.matricesA))
  for (i in 1:length(lista.matricesA)){
      dif.matrix.list[[i]] <- lista.matricesA[[i]] - lista.matricesB[[i]]
      frobenius.norma[i] <- norm(dif.matrix.list[[i]],type="F")
  }
  # Output de la funcion
  output <- list(matriz.dif = dif.matrix.list, norma.diferencia = frobenius.norma)
  return(output)
}