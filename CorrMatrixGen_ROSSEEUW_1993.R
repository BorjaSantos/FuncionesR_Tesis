######################################################
######################################################
#####TRANSFORMACION MATRICES SEGUN ROSSEEUW_1993######
######################################################
######################################################


#######################################################
#### NOTA: CAMBIAR 0.00001 por epsilon en las funciones,
#### así queda más flexible la elección de la constante
#### de la transformación. 
#######################################################
# FUNCIONES PARA MUESTREAR MATRICES DE CORRELACION

# FECHA INICIO: 1 / SEPTIEMBRE / 2014.
# FECHA FIN: 5/ SEPTIRMBRE / 2014.

# CASO METAANALISIS MULTIVARIADO CON 2 VARIABLES DE RESULTADO

# ARGUMENTOS DE ENTRADA (INPUTS):
# i) Tamaño de la simulacion (N)
# ii) Distribucion de probabilidad asociada a cada coeficiente de correlacion ((n)(n+1))/2)

# ARGUMENTOS DE SALIDA (OUTPUT):
# i) Lista de matrices de correlacion.

# QÚE HACE LA FUNCION:
# Paso 1: Indicar el tamaño de la matriz de correlacion y determinar la distribucion de parobabilidad de cada elemento.
# Paso 2: Generar una muestra aleatoria de tamaño N con cada distribucion de probabilidad.
# Paso 3: Construir N matrices simetricas con cada elemento de las muestras generadas en el paso anterior.

# Libreria necesaria para poder construir matrices simetricas a partir de un vector con sus elementos.
library(miscTools) 

corr.matrix.gen.2vars <- function(N,coefficients,semilla){
  set.seed(semilla)
  # CRITERIOS DE CONTROL DE LOS PARMETROS DE LA FUNCION
  if(N <=0){
    stop("El tamano de la muestra aleatoria de matrices de correlacion no es correcta")  
  }
  if(length(coefficients) != N){
    stop("El numero de coeficientes de correlacion no es correcto")
  }
  # CREACION DE LA LISTA DE MATRICES
  matrix.list <- vector("list",N)
  for (i in 1:N){
    matrix.list[[i]] <- symMatrix(c(1,coefficients[i],1),byrow = FALSE)
  }
  return(matrix.list)
}

# INPUT: 
# Lista de matrices simetricas

# OUTPUT
# Lista de matrices de correlacion.

# QUE HACE ESTA FUNCION
# Paso 1: Comprobar que son semidefinidas positivas, si no lo son, transformarlas para que si lo sean.
# Paso 2: Devolver una lista de tamaño N de matrices de correlacion.

transformar.matrices.2vars <- function(lista.matrices){
  lista.transformadas <- vector("list",length(lista.matrices))
  for(i in 1:length(lista.matrices)){
    determinante <- det(lista.matrices[[i]])
    if (determinante < 0){
      descomposicion <- eigen(lista.matrices[[i]])
      for (j in 1:length(descomposicion$values)){
        if(descomposicion$values[j] < 0){
          descomposicion$values[j] <- 0.00001
        }
        m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
        m.aux2 <- diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[1,1])))%*%m.aux%*%diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2])))
      }
      lista.transformadas[[i]] <- m.aux2  
    }else{
      lista.transformadas[[i]] <- lista.matrices[[i]]  
    }
  }
  return(lista.transformadas)
}

# CASO METAANALISIS MULTIVARIADO CON 3 VARIABLES DE RESULTADO

# ARGUMENTOS DE ENTRADA (INPUTS):
# i) Tamaño de la simulacion (N)
# ii) Distribucion de probabilidad asociada a cada coeficiente de correlacion ((n)(n+1))/2)

# ARGUMENTOS DE SALIDA (OUTPUT):
# i) Lista de matrices de correlacion.

# QÚE HACE LA FUNCION:
# Paso 1: Indicar el tamaño de la matriz de correlacion y determinar la distribucion de parobabilidad de cada elemento.
# Paso 2: Generar una muestra aleatoria de tamaño N con cada distribucion de probabilidad.
# Paso 3: Construir N matrices simetricas con cada elemento de las muestras generadas en el paso anterior.

# Libreria necesaria para poder construir matrices simetricas a partir de un vector con sus elementos.
library(miscTools) 

corr.matrix.gen.3vars <- function(N,coefficients,semilla){
    set.seed(semilla)
    # CRITERIOS DE CONTROL DE LOS PARMETROS DE LA FUNCION
    if(N <=0){
      stop("El tamano de la muestra aleatoria de matrices de correlacion no es correcta")  
    }
    if(length(coefficients) != N*3){
      stop("El numero de coeficientes de correlacion no es correcto")
    }
    # CREACION DE LA LISTA DE MATRICES
    matrix.list <- vector("list",N)
    for (i in 1:N){
      matrix.list[[i]] <- symMatrix(c(1,coefficients[i],coefficients[i+N],1,coefficients[i+(2*N)],1),byrow = FALSE)
    }
    return(matrix.list)
}

# INPUT: 
# Lista de matrices simetricas

# OUTPUT
# Lista de matrices de correlacion.

# QUE HACE ESTA FUNCION
# Paso 1: Comprobar que son semidefinidas positivas, si no lo son, transformarlas para que si lo sean.
# Paso 2: Devolver una lista de tamaño N de matrices de correlacion.

transformar.matrices.3vars <- function(lista.matrices){
  lista.transformadas <- vector("list",length(lista.matrices))
  for(i in 1:length(lista.matrices)){
    determinante <- det(lista.matrices[[i]])
    if (determinante < 0){
      descomposicion <- eigen(lista.matrices[[i]])
      for (j in 1:length(descomposicion$values)){
        if(descomposicion$values[j] < 0){
          descomposicion$values[j] <- 0.00001
        }
        m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
        m.aux2 <- diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3])))%*%m.aux%*%diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3])))
      }
      lista.transformadas[[i]] <- m.aux2  
    }else{
      lista.transformadas[[i]] <- lista.matrices[[i]]  
    }
  }
  return(lista.transformadas)
}

# CASO METAANALISIS MULTIVARIADO CON 4 VARIABLES DE RESULTADO

# ARGUMENTOS DE ENTRADA (INPUTS):
# i) Tamaño de la simulacion (N)
# ii) Distribucion de probabilidad asociada a cada coeficiente de correlacion ((n)(n+1))/2)

# ARGUMENTOS DE SALIDA (OUTPUT):
# i) Lista de matrices de correlacion.

# QÚE HACE LA FUNCION:
# Paso 1: Indicar el tamaño de la matriz de correlacion y determinar la distribucion de parobabilidad de cada elemento.
# Paso 2: Generar una muestra aleatoria de tamaño N con cada distribucion de probabilidad.
# Paso 3: Construir N matrices simetricas con cada elemento de las muestras generadas en el paso anterior.

# Libreria necesaria para poder construir matrices simetricas a partir de un vector con sus elementos.
library(miscTools) 

corr.matrix.gen.4vars <- function(N,coefficients,semilla){
  set.seed(semilla)
  # CRITERIOS DE CONTROL DE LOS PARMETROS DE LA FUNCION
  if(N <=0){
    stop("El tamano de la muestra aleatoria de matrices de correlacion no es correcta")  
  }
  if(length(coefficients) != N*6){
    stop("El numero de coeficientes de correlacion no es correcto")
  }
  # CREACION DE LA LISTA DE MATRICES
  matrix.list <- vector("list",N)
  for (i in 1:N){
    matrix.list[[i]] <- symMatrix(c(1,coefficients[i],coefficients[i+N],coefficients[i+2*N],
                                    1,coefficients[i+3*N],coefficients[i+4*N],
                                    1,coefficients[i+5*N],
                                    1),byrow = FALSE)
  }
  return(matrix.list)
}

# INPUT: 
# Lista de matrices simetricas

# OUTPUT
# Lista de matrices de correlacion.

# QUE HACE ESTA FUNCION
# Paso 1: Comprobar que son semidefinidas positivas, si no lo son, transformarlas para que si lo sean.
# Paso 2: Devolver una lista de tamaño N de matrices de correlacion.

transformar.matrices.4vars <- function(lista.matrices){
  lista.transformadas <- vector("list",length(lista.matrices))
  for(i in 1:length(lista.matrices)){
    determinante <- det(lista.matrices[[i]])
    if (determinante < 0){
      descomposicion <- eigen(lista.matrices[[i]])
      for (j in 1:length(descomposicion$values)){
        if(descomposicion$values[j] < 0){
          descomposicion$values[j] <- 0.00001
        }
        m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
        m.aux2 <- diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3]),1/sqrt(m.aux[4,4])))%*%m.aux%*%diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3]),1/sqrt(m.aux[4,4])))
      }
      lista.transformadas[[i]] <- m.aux2  
    }else{
      lista.transformadas[[i]] <- lista.matrices[[i]]  
    }
  }
  return(lista.transformadas)
}


# CASO METAANALISIS MULTIVARIADO CON 5 VARIABLES DE RESULTADO

# ARGUMENTOS DE ENTRADA (INPUTS):
# i) Tamaño de la simulacion (N)
# ii) Distribucion de probabilidad asociada a cada coeficiente de correlacion ((n)(n+1))/2)

# ARGUMENTOS DE SALIDA (OUTPUT):
# i) Lista de matrices de correlacion.

# QÚE HACE LA FUNCION:
# Paso 1: Indicar el tamaño de la matriz de correlacion y determinar la distribucion de parobabilidad de cada elemento.
# Paso 2: Generar una muestra aleatoria de tamaño N con cada distribucion de probabilidad.
# Paso 3: Construir N matrices simetricas con cada elemento de las muestras generadas en el paso anterior.

# Libreria necesaria para poder construir matrices simetricas a partir de un vector con sus elementos.
library(miscTools) 

corr.matrix.gen.5vars <- function(N,coefficients,semilla){
  set.seed(semilla)
  # CRITERIOS DE CONTROL DE LOS PARMETROS DE LA FUNCION
  if(N <=0){
    stop("El tamano de la muestra aleatoria de matrices de correlacion no es correcta")  
  }
  if(length(coefficients) != N*10){
    stop("El numero de coeficientes de correlacion no es correcto")
  }
  # CREACION DE LA LISTA DE MATRICES
  matrix.list <- vector("list",N)
  for (i in 1:N){
    matrix.list[[i]] <- symMatrix(c(1,coefficients[i],coefficients[i+N],coefficients[i+2*N],coefficients[i+3*N],
                                    1,coefficients[i+4*N],coefficients[i+5*N],coefficients[i+6*N],
                                    1,coefficients[i+7*N],coefficients[i+8*N],
                                    1,coefficients[i+9*N],
                                    1),byrow = FALSE)
  }
  return(matrix.list)
}

# INPUT: 
# Lista de matrices simetricas

# OUTPUT
# Lista de matrices de correlacion.

# QUE HACE ESTA FUNCION
# Paso 1: Comprobar que son semidefinidas positivas, si no lo son, transformarlas para que si lo sean.
# Paso 2: Devolver una lista de tamaño N de matrices de correlacion.

transformar.matrices.5vars <- function(lista.matrices){
  lista.transformadas <- vector("list",length(lista.matrices))
  for(i in 1:length(lista.matrices)){
    determinante <- det(lista.matrices[[i]])
    if (determinante < 0){
      descomposicion <- eigen(lista.matrices[[i]])
      for (j in 1:length(descomposicion$values)){
        if(descomposicion$values[j] < 0){
          descomposicion$values[j] <- 0.00001
        }
        m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
        m.aux2 <- diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3]),1/sqrt(m.aux[4,4]),1/sqrt(m.aux[5,5])))%*%m.aux%*%diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3]),1/sqrt(m.aux[4,4]),1/sqrt(m.aux[5,5])))
      }
      lista.transformadas[[i]] <- m.aux2  
    }else{
      lista.transformadas[[i]] <- lista.matrices[[i]]  
    }
  }
  return(lista.transformadas)
}

# CASO METAANALISIS MULTIVARIADO CON 6 VARIABLES DE RESULTADO

# ARGUMENTOS DE ENTRADA (INPUTS):
# i) Tamaño de la simulacion (N)
# ii) Distribucion de probabilidad asociada a cada coeficiente de correlacion ((n)(n+1))/2)

# ARGUMENTOS DE SALIDA (OUTPUT):
# i) Lista de matrices de correlacion.

# QÚE HACE LA FUNCION:
# Paso 1: Indicar el tamaño de la matriz de correlacion y determinar la distribucion de parobabilidad de cada elemento.
# Paso 2: Generar una muestra aleatoria de tamaño N con cada distribucion de probabilidad.
# Paso 3: Construir N matrices simetricas con cada elemento de las muestras generadas en el paso anterior.

# Libreria necesaria para poder construir matrices simetricas a partir de un vector con sus elementos.
library(miscTools) 

corr.matrix.gen.6vars <- function(N,coefficients,semilla){
  set.seed(semilla)
  # CRITERIOS DE CONTROL DE LOS PARMETROS DE LA FUNCION
  if(N <=0){
    stop("El tamano de la muestra aleatoria de matrices de correlacion no es correcta")  
  }
  if(length(coefficients) != N*15){
    stop("El numero de coeficientes de correlacion no es correcto")
  }
  # CREACION DE LA LISTA DE MATRICES
  matrix.list <- vector("list",N)
  for (i in 1:N){
    matrix.list[[i]] <- symMatrix(c(1,coefficients[i],coefficients[i+N],coefficients[i+2*N],coefficients[i+3*N],coefficients[i+4*N],
                                    1,coefficients[i+5*N],coefficients[i+6*N],coefficients[i+7*N],coefficients[i+8*N],
                                    1,coefficients[i+9*N],coefficients[i+10*N],coefficients[i+11*N],
                                    1,coefficients[i+12*N],coefficients[i+13*N],
                                    1,coefficients[i+14*N],
                                    1),byrow = FALSE)
  }
  return(matrix.list)
}

# INPUT: 
# Lista de matrices simetricas

# OUTPUT
# Lista de matrices de correlacion.

# QUE HACE ESTA FUNCION
# Paso 1: Comprobar que son semidefinidas positivas, si no lo son, transformarlas para que si lo sean.
# Paso 2: Devolver una lista de tamaño N de matrices de correlacion.

transformar.matrices.6vars <- function(lista.matrices){
  lista.transformadas <- vector("list",length(lista.matrices))
  for(i in 1:length(lista.matrices)){
    determinante <- det(lista.matrices[[i]])
    if (determinante < 0){
      descomposicion <- eigen(lista.matrices[[i]])
      for (j in 1:length(descomposicion$values)){
        if(descomposicion$values[j] < 0){
          descomposicion$values[j] <- 0.00001
        }
        m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
        m.aux2 <- diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3]),1/sqrt(m.aux[4,4]),1/sqrt(m.aux[5,5]),1/sqrt(m.aux[6,6])))%*%m.aux%*%diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3]),1/sqrt(m.aux[4,4]),1/sqrt(m.aux[5,5]),1/sqrt(m.aux[6,6])))
      }
      lista.transformadas[[i]] <- m.aux2  
    }else{
      lista.transformadas[[i]] <- lista.matrices[[i]]  
    }
  }
  return(lista.transformadas)
}
























# CHEQUEO DE QUE LAS MATRICES GENERADAS SEAN DE CORRELACION (SEMIDEFINIDAS POSITIVAS)
# SI matrix.list[[i]] TIENE DETERMINANTE >=0 PERFECTO.
# SI NO, DEBE SER MODIFICADA Y ENCONTRAR LA MATRIZ MAS PROXIMA A ELLA QUE LO SEA.
convertir.matriz.correlacion <- function(matriz.lista){
  matrix.list <- matriz.lista 
  lista.modificadas <- vector("list",length(matriz.lista))
  print("Longitud lista matrices")
  print(length(matrix.list))
  for (i in 1:length(matrix.list)){
    print("Matriz nuemro...")
    print(i)
    print(matrix.list[[i]])
    d <- det(matrix.list[[i]])
    print("Determinante de la matriz")
    print(d)
    if (d < 0){
      print("guay")
      desc <- eigen(matrix.list[[i]])
      print(desc$values)
      for (i in 1:length(desc$values)){
        if(desc$values[i] < 0){
          print("guay2")
          desc$values[i] <- 0.00001
        }  
      }
      print(desc$values)
      m <- desc$vectors%*%diag(desc$values)%*%t(desc$vectors)
      m2 <- diag(c(1/sqrt(m[1,1]),1/sqrt(m[2,2]),1/sqrt(m[3,3])))%*%m%*%diag(c(1/sqrt(m[1,1]),1/sqrt(m[2,2]),1/sqrt(m[3,3])))         
      print("Matriz transformada porque tocaba")
      print(m2)
      print("Comprobacion de determinante mayor que cero")
      print(det(m2))
      lista.modificadas[[i]] <- m2
    }
    if(d>=0){
      lista.modificadas[[i]] <- matriz.lista[[i]]
    }
      
  }
  print("Lista de matrices original")
  print(matriz.lista)
  print("Lista de matrices modificadas")
  print(lista.modificadas)
  #return(matrix.list)
}
# OUTPUT DE LA FUNCION

# comando pruebas:
prueba <- corr.matrix.gen.3vars(5,c(rbeta(5,0.5,0.6),runif(5,0.2,1),runif(5,-1,1)),2014)

transformar.matrices <- function(lista.matrices){
  lista.transformadas <- vector("list",length(lista.matrices))
  for(i in 1:length(lista.matrices)){
    determinante <- det(lista.matrices[[i]])
    if (determinante < 0){
      descomposicion <- eigen(lista.matrices[[i]])
      for (j in 1:length(descomposicion$values)){
        if(descomposicion$values[j] < 0){
          descomposicion$values[j] <- 0.00001
        }
      m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
      m.aux2 <- diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3])))%*%m.aux%*%diag(c(1/sqrt(m.aux[1,1]),1/sqrt(m.aux[2,2]),1/sqrt(m.aux[3,3])))
      }
    lista.transformadas[[i]] <- m.aux2  
    }else{
    lista.transformadas[[i]] <- lista.matrices[[i]]  
    }
  }
  return(lista.transformadas)
}

