######################################################
######################################################
########## TRANSFORMACION MATRICES SEGUN #############
########### REBONATO & JACKEL 1999  ##################
######################################################
######################################################

# FUNCIONES PARA MUESTREAR MATRICES DE CORRELACION

# FECHA INICIO: 4 / NOVIEMBRE / 2014.
# FECHA FIN: 6 / NOVIEMBRE / 2014.


# INPUT: 
# Lista de matrices simetricas

# OUTPUT
# Lista de matrices de correlacion.

# QUE HACE ESTA FUNCION
# Paso 1: Comprobar que son semidefinidas positivas, si no lo son, transformarlas para que si lo sean.
# Paso 2: Devolver una lista de tamaño N de matrices de correlacion.


rebonato.jackel.espectral <- function(lista.matrices){
  matrices.trasformadas <- vector("list",length(lista.matrices))
  matriz.t <- diag(0,dim(matrices.correlacion[[1]]))
  for(i in 1:length(lista.matrices)){
    determinante <- det(lista.matrices[[i]])
    if (determinante < 0){
      descomposicion <- eigen(lista.matrices[[i]])
      for (j in 1:length(descomposicion$values)){
        if(descomposicion$values[j] < 0){
          descomposicion$values[j] <- 0
        }
        delta <- diag(descomposicion$values) #corregir
        sqrt.delta <- diag(sqrt(descomposicion$values))
        b.prima <- eigen(lista.matrices[[i]])$vectors%*%sqrt.delta
        for(k in 1:dim(delta)[1]) {
          matriz.t[k,k] <- 1/(sum(((descomposicion$vectors[k,])^2*diag(delta))))
        }
        sqrt.matriz.t <-diag(sqrt(diag(matriz.t)))
        b <- sqrt.matriz.t%*%b.prima
        m.aux <- b%*%t(b)
      }
      matrices.trasformadas[[i]] <- m.aux 
    }else{
      matrices.trasformadas[[i]] <- lista.matrices[[i]]  
    }
  }
  matrices.trasformadas <- suppressWarnings(matrices.trasformadas)
  return(matrices.trasformadas)
}

