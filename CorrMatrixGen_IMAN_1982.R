######################################################
######################################################
#####  TRANSFORMACION MATRICES SEGUN IMAN_1982  ######
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


trasformar.matrices.iman2 <- function(lista.matrices,epsilon,opcion){
  if(opcion=="A"){
  output <- lapply(lista.matrices,opcion.1)
  }
  if(opcion=="B"){
    output <- lapply(lista.matrices,opcion.2)
  }
  if(opcion=="C"){
    output <- lapply(lista.matrices,opcion.3)
  }
  if(opcion=="D"){
    output <- lapply(lista.matrices,opcion.4)
  }
  if(opcion=="E"){
    output <- lapply(lista.matrices,opcion.5)
  }
  
  return(output)
}



############################### FUNCIONES AUXILIARES##############################################

# METODO DE CORRECCION A (iGUALAR A 1 LA DIAGONAL DE LAS MATRICES)

opcion.1 <- function(a){
  while(det(a)<=0){
   descomposicion <- eigen(a)
   for (j in 1:length(descomposicion$values)){
         if(descomposicion$values[j] < epsilon){
               descomposicion$values[j] <- epsilon
           }
         m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
         diag(m.aux) <- rep(1,dim(m.aux)[1])
  }
  a <- m.aux
}
return(a)
}

##################################################################################################

##################################################################################################

# METODO DE CORRECCION B (iGUALAR A 1 LA DIAGONAL DE LAS MATRICES 
# Y EL RESTO DE COMPONENTES A +/- 0.9999)

opcion.2 <- function(a){
  while(det(a)<=0){
    descomposicion <- eigen(a)
    for (j in 1:length(descomposicion$values)){
      if(descomposicion$values[j] < epsilon){
        descomposicion$values[j] <- epsilon
      }
      m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
      for(k in 1:dim(m.aux)[1]){
        for(l in 1:dim(m.aux)[2]){
          if(m.aux[k,l] <= -1){
            m.aux[k,l] <- -0.999
          }
          if(m.aux[k,l] >= 1){
            m.aux[k,l] <- 0.999
          }
        }
      }
      diag(m.aux) <- rep(1,dim(m.aux)[1])
      #print(m.aux)
    }
    a <- m.aux
  }
  return(a)
}

################################################################################################

################################################################################################

# METODO DE CORRECCION C (a'[i,j]=a[i,j]/(sqrt(r[i,i]*r[j,j])))

opcion.3 <- function(a){
  while(det(a)<=0){
  descomposicion <- eigen(a)
  for (j in 1:length(descomposicion$values)){
    if(descomposicion$values[j] < epsilon){
      descomposicion$values[j] <- epsilon
    }
    m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
    m.aux2 <- m.aux
    for(k in 1:dim(m.aux)[1]){
      for(l in 1:dim(m.aux)[2]){
        m.aux[k,l] <- m.aux2[k,l]/(sqrt(m.aux2[k,k]*m.aux2[l,l]))
      }
    }
    #print(m.aux)
  }
  a <- m.aux
}
return(a)
}

###############################################################################################

###############################################################################################
# METODO DE CORRECCION D (igual que C pero forzando la diagonal a 1)

opcion.4 <- function(a){
  while(det(a)<=0){
    descomposicion <- eigen(a)
    for (j in 1:length(descomposicion$values)){
      if(descomposicion$values[j] < epsilon){
        descomposicion$values[j] <- epsilon
      }
      m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
      m.aux2 <- m.aux
      for(k in 1:dim(m.aux)[1]){
        for(l in 1:dim(m.aux)[2]){
          m.aux[k,l] <- m.aux2[k,l]/(sqrt(m.aux2[k,k]*m.aux2[l,l]))
        }
      }
    }
    a <- m.aux
  }
  diag(a) <- rep(1,dim(a)[1])
  return(a)
}

###############################################################################################
# METODO DE CORRECCION E (a'[i,j]=a[i,j]/(sqrt((r[i,i]-0.05)*(r[j,j]-0.05))))

opcion.5 <- function(a){
  while(det(a)<=0){
    descomposicion <- eigen(a)
    for (j in 1:length(descomposicion$values)){
      if(descomposicion$values[j] < epsilon){
        descomposicion$values[j] <- epsilon
      }
      m.aux <- descomposicion$vectors%*%diag(descomposicion$values)%*%t(descomposicion$vectors)
      m.aux2 <- m.aux
      for(k in 1:dim(m.aux)[1]){
        for(l in 1:dim(m.aux)[2]){
          m.aux[k,l] <- m.aux2[k,l]/(sqrt((diag(m.aux2)[k]-0.05))*sqrt((diag(m.aux2)[l]-0.05)))
        }
      }
    }
    a <- m.aux
  }
  # Cambiar diagonal a a 1.
  diag(a) <- rep(1,dim(a)[1])
  return(a)
}




















# EJEMPLO:
# ESTIMADOR DEL EFECTO VARIABLE 1
smd1 <- rnorm(6,0.1) 
var1 <- runif(6,5,20)

# ESTIMADOR DEL EFECTOS VARIABLE 2
smd2 <- rbeta(6,0.1,4) 
var2 <- runif(6,5,20)

# ESTIMADOR DEL EFECTO VARIABLE 3
smd3 <- runif(6,-2,2) 
var3 <- runif(6,5,20)

# BASE DE DATOS:
estudio <- c(1:6)

datos <- cbind(estudio,smd1,var1,smd2,var2,smd3,var3)
datos

# Calculo de las matrices de varianzas-covarianzas VCOV = S*R*S

# Matrices de desviaciones estandar
sd1 <- diag(c(sqrt(datos[1,3]),sqrt(datos[1,5]),sqrt(datos[1,7])))
sd2 <- diag(c(sqrt(datos[2,3]),sqrt(datos[2,5]),sqrt(datos[2,7])))
sd3 <- diag(c(sqrt(datos[3,3]),sqrt(datos[3,5]),sqrt(datos[3,7])))
sd4 <- diag(c(sqrt(datos[4,3]),sqrt(datos[4,5]),sqrt(datos[4,7])))
sd5 <- diag(c(sqrt(datos[5,3]),sqrt(datos[5,5]),sqrt(datos[5,7])))
sd6 <- diag(c(sqrt(datos[6,3]),sqrt(datos[6,5]),sqrt(datos[6,7])))
sd <- list(sd1,sd2,sd3,sd4,sd5,sd6)

# Matrices de correlacion:
matrices.correlacion <- corr.matrix.gen.3vars(5,c(rbeta(5,0.5,0.6),runif(5,0.2,1),runif(5,-1,1)),2014)
lapply(matrices.correlacion,det)

prueba <- transformar.matrices.iman(matrices.correlacion,0.001,"A")
det(prueba[[1]])

