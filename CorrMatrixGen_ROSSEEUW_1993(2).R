                               ######################################################
                               ######################################################
                               ########  TRANSFORMACION MATRICES SEGUN   ############
                               ############ ROUSSEEUW & MOLENBERGHS #################
                               ##################### 1993 ###########################
                               ######################################################
                               ######################################################

# FECHA INICIO: 6 / NOVIEMBRE / 2014.
# FECHA FIN:  / NOVIEMBRE / 2014.


# INPUT: 
# Lista de matrices simetricas

# OUTPUT
# Lista de matrices de correlacion.

# QUE HACE ESTA FUNCION
# Paso 1: Comprobar que son semidefinidas positivas, si no lo son, transformarlas para que si lo sean.
# Paso 2: Devolver una lista de tamaño N de matrices de correlacion.

# SHRINKING METHODS

###############################################################################################################################
shrink.lineal <- function(lista.matrices){
  list.output <- lapply(lista.matrices,shrink.lineal.aux)
  return(list.output) # NOTA: si aparece la matriz identidad es porque no habia valores propios en [1,0]
}
  
# Funciones auxiliares:

shrink.lineal.aux <- function(a){
  if(det(a)<0){
  # Descomposicion de la matriz en R = T*D*T^{t}
  valores.propios <- eigen(a)$values
  max.val.propio <- aux(valores.propios)
  a.new <- max.val.propio*a+(1-max.val.propio)*diag(dim(a)[1])
  return(a.new)
  }else{
    return(a)
  }  
}   
aux <- function(v){
  a <- which(v<=1)[1]
  val <- v[a]
  if(val>0){
    return(val)
  }else{
    val <- 0
    return(val)
  }
}
###############################################################################################################################  
  
###############################################################################################################################

shrink.no.lineal <- function(lista.matrices,epsilon,link){
  epsilon <- epsilon
  if(link=="A"){
  output <- lapply(lista.matrices,opcion.1)  
  }
  if(link=="B"){
    output <- lapply(lista.matrices,opcion.2)  
  }
  return(output)
}
                               
# Funciones auxiliares:

opcion.1 <-  function(a){
  while(det(a)<=0){
  for(i in 1:dim(a)[1]){
    for(j in 1:dim(a)[2]){
      if(a[i,j]< (-atanh(epsilon))){
        a[i,j] <- atanh(tanh(a[i,j])+epsilon)
      }
      if(abs(a[i,j])< atanh(epsilon)){
        a[i,j] <- 0
      }
      if(a[i,j] > atanh(epsilon)){
        a[i,j] <- atanh(tanh(a[i,j])-epsilon)
      }
    }
  }
  diag(a) <- rep(1,dim(a)[1])
  }
  return(a)
}  

opcion.2 <- function(a){
  while(det(a)<=0){
    for(i in 1:dim(a)[1]){
      for(j in 1:dim(a)[2]){
        if(a[i,j]< (-1)*atanh(epsilon)){
          a[i,j] <- atanh(tanh(a[i,j])+epsilon)
        }
        if(abs(a[i,j])< atanh(epsilon)){
          a[i,j] <- 0
        }
        if(a[i,j] > atanh(epsilon)){
          a[i,j] <- atanh(tanh(a[i,j])-epsilon)
        }
      }
    }
    diag(a) <- rep(1,dim(a)[1])
  }
  return(a)  
}

##############################################################################################################################
##############################################################################################################################
# FALTA SHRINKING METHOD