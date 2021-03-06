################################################################################################################################
#################################FUNCION 1: SIMULACION DE ESTUDIOS CON DATOS INDIVIDUALES#######################################
################################################################################################################################

# PARAMTEROS IMPORTANTES EN LA SIMULACION:
# 1. Numero de estudios
# 2. Numero de variables
# 3. Tamano de la muestra
# 4. Media grupo tratamiento
# 5. SD grupo tratamiento
# 6. Media grupo control
# 7. SD grupo control.
# 8. Correlacion entre outcomes.

# VALORES DE LOS PARAMETROS
# Numero de estudios: 5, 10, 15, 25, 50.
# Numero de variables: 2, 4, 6.
# Tamano de la muestra: 20, 40, 60, 80, 100 (la mitad en cada rama) 
# m1.trat / m1.trat: media grupo tratamiento / media grupo control para la variable 1.
# sd1.trat / sd1.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 1.
# a1 / b1: limites superior e inferior del rango de valores de la variable 1.
# m2.trat / m2.trat: media grupo tratamiento / media grupo control para la variable 2.
# sd2.trat / sd2.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 2.
# a2 / b2: limites superior e inferior del rango de valores de la variable 2.

# Librerias necesarias para la simulacion
library(truncdist) # De esta manera se truncaran los valores de las diferentes puntuaciones, para que no se simulen puntuaciones
# fuera de rango.
# Semilla para poder replicar los datos
set.seed(18052013)

# Funcion que va a simular los datos en el escenario que queramos

simulacion_datos_2vars <- function(n.estudios,m1.trat,sd1.trat,m1.ctrl,sd1.ctrl,a1,b1,m2.trat,sd2.trat,m2.ctrl,sd2.ctrl,a2,b2,
                                   tamano.muestra,semilla,replicaciones){
  # Control de los paramtero de la funcion
  if(n.estudios > 50 || n.estudios < 5){
    stop("Numero de estudios incorrecto")
  }
  if(tamano.muestra > 100 || tamano.muestra < 20){
    stop("Tamano de muestra incorrecto")
  }
  if(replicaciones < 5 || replicaciones > 150){
    stop("El numero de replicaciones es incorrecto")
  }
  # Libreria necesaria para truncar las distribuciones
  library(truncdist)
  # Funcion propiamente dicha
  database <- vector("list",replicaciones) # Lista donde se almacenaran las replicaciones del escenario deseado
  for(i in 1:replicaciones){
    database[[i]] <- vector("list",n.estudios) # Lista donde se almacenaran el numero de estudios determinado
    set.seed(semilla)
    for(j in 1:n.estudios){
      Var1.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.trat,sd=sd1.trat),0) # Variable 1 tratamiento.
      Var1.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.ctrl,sd=sd1.ctrl),0) # Variable 1 control.
      Var2.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.trat,sd=sd2.trat),0) # Variable 2 tratamiento.
      Var2.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.ctrl,sd=sd2.ctrl),0) # Variable 2 control.
      database[[i]][[j]] <- as.data.frame(cbind(Var1.Trat,Var1.Ctrl,Var2.Trat,Var2.Ctrl))
      semilla <- semilla + 1
    }
  }
  
  return(database)
}

###############################################################################################################################
###################FUNCION 2: SIMULACION DE ESTUDIOS CON DATOS INDIVIDUALES (SIMPLIFICADA)#####################################
###############################################################################################################################

# COMPARACION DE LOS MÉTODOS DE METAANALISIS DESARROLLADOS FRENTE A LOS PROPUESTOS EN LA TESIS.
# LA COMPARACION SE HARA MEDIANTE SIMULACION BAJO DIFERENTES CONDICIONES.

# FECHA INICIO: 4 / MARZO / 2014
# FECHA FIN: 10/ MARZO / 2014

# PARAMTEROS IMPORTANTES EN LA SIMULACION:
# 1. Numero de estudios
# 2. Numero de variables
# 3. Tamano de la muestra
# 4. Media grupo tratamiento
# 5. SD grupo tratamiento
# 6. Media grupo control
# 7. SD grupo control.
# 8. Correlacion entre outcomes.

# VALORES DE LOS PARAMETROS
# Numero de estudios: 5, 10, 15, 25, 50.
# Numero de variables: 2, 4, 6.
# Tamano de la muestra: 20, 40, 60, 80, 100 (la mitad en cada rama) 
# Numero de estudios: 5, 10, 15, 25, 50.
# Numero de variables: 2, 4, 6.
# Tamano de la muestra: 20, 40, 60, 80, 100 (la mitad en cada rama) 
# m1.trat / m1.trat: media grupo tratamiento / media grupo control para la variable 1.
# sd1.trat / sd1.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 1.
# a1 / b1: limites superior e inferior del rango de valores de la variable 1.
# m2.trat / m2.trat: media grupo tratamiento / media grupo control para la variable 2.
# sd2.trat / sd2.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 2.
# a2 / b2: limites superior e inferior del rango de valores de la variable 2.

# Librerias necesarias para la simulacion
library(truncdist) # De esta manera se truncaran los valores de las diferentes puntuaciones, para que no se simulen puntuaciones
# fuera de rango.
# Semilla para poder replicar los datos
set.seed(18052013)

# Funcion que va a simular los datos en el escenario que queramos

simulacion_datos2_2vars <- function(n.estudios,m1.trat,sd1.trat,m1.ctrl,sd1.ctrl,a1,b1,m2.trat,sd2.trat,m2.ctrl,sd2.ctrl,a2,b2,
                                    tamano.muestra,semilla,replicaciones){
  # Control de los paramtero de la funcion
  if(n.estudios > 50 || n.estudios < 5){
    stop("Numero de estudios incorrecto")
  }
  if(tamano.muestra > 100 || tamano.muestra < 20){
    stop("Tamano de muestra incorrecto")
  }
  if(replicaciones < 5 || replicaciones > 150){
    stop("El numero de replicaciones es incorrecto")
  }
  # Funcion propiamente dicha
  database <- vector("list",replicaciones) # Lista donde se almacenaran las replicaciones del escenario deseado
  for(i in 1:replicaciones){
    database[[i]] <- vector("list",n.estudios) # Lista donde se almacenaran el número de estudios determinado
    set.seed(semilla)
    Var1.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.trat,sd=sd1.trat),0) # Variable 1 tratamiento.
    Var1.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.ctrl,sd=sd1.ctrl),0) # Variable 1 control.
    Var2.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.trat,sd=sd2.trat),0) # Variable 2 tratamiento.
    Var2.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.ctrl,sd=sd2.ctrl),0) # Variable 2 control.
    database[[i]] <- as.data.frame(cbind(Var1.Trat,Var1.Ctrl,Var2.Trat,Var2.Ctrl))
    semilla <- semilla + (n.estudios)
  }
  return(database)
}

##############################################################################################################################
####################FUNCION 3: g DE HEDGES' DE UNA SIMULACION CON REPLICACIONES (SIMPLIFICADA)################################
##############################################################################################################################

g_hedges_simulacion2_2vars <- function(database,correlacion=0){
  library(compute.es)
  datos_meta <- diag(0,length(database),5)
  for(i in 1:length(database)){
    g1 <- mes(mean(database[[i]][,1]),sd(database[[i]][,1]),dim(database[[i]])[1],mean(database[[i]][,2]),sd(database[[i]][,2]),dim(database[[i]])[1])[12][,1]
    g2 <- mes(mean(database[[i]][,3]),sd(database[[i]][,3]),dim(database[[i]])[1],mean(database[[i]][,4]),sd(database[[i]][,4]),dim(database[[i]])[1])[12][,1]
    var.g1 <- mes(mean(database[[i]][,1]),sd(database[[i]][,1]),dim(database[[i]])[1],mean(database[[i]][,2]),sd(database[[i]][,2]),dim(database[[i]])[1])[13][,1]
    var.g2 <- mes(mean(database[[i]][,3]),sd(database[[i]][,3]),dim(database[[i]])[1],mean(database[[i]][,4]),sd(database[[i]][,4]),dim(database[[i]])[1])[13][,1]
    covar.g1g2 <- correlacion*sqrt(var.g1)*sqrt(var.g2)
    input <- c(g1,g2,var.g1,covar.g1g2,var.g2)
    datos_meta[i,] <- input
  }
  colnames(datos_meta) <- c("g1","g2","var.g1","covar.g1g2","var.g2")
  return(as.data.frame(datos_meta))
}

##############################################################################################################################
##########################FUNCION 4: g DE HEDGES' DE UNA SIMULACION CON REPLICACIONES#########################################
##############################################################################################################################

g_hedges_simulacion_2vars <- function(database,correlacion=0){
  library(compute.es)
  datos_meta <- vector("list",length(database))
  g1 <- vector(mode="numeric",length=length(database[[1]]))
  g2 <- vector(mode="numeric",length=length(database[[1]]))
  var.g1 <- vector(mode="numeric",length=length(database[[1]]))
  var.g2 <- vector(mode="numeric",length=length(database[[1]]))
  covar.g1g2 <- vector(mode="numeric",length=length(database[[1]])) 
  for(i in 1:length(database)){
    for(j in 1:length(database[[i]])){
      g1[j] <- mes(mean(database[[i]][[j]][,1]),sd(database[[i]][[j]][,1]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,2]),sd(database[[i]][[j]][,2]),dim(database[[i]][[j]])[1])[12][,1]
      g2[j] <- mes(mean(database[[i]][[j]][,3]),sd(database[[i]][[j]][,3]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,4]),sd(database[[i]][[j]][,4]),dim(database[[i]][[j]])[1])[12][,1]
      var.g1[j] <- mes(mean(database[[i]][[j]][,1]),sd(database[[i]][[j]][,1]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,2]),sd(database[[i]][[j]][,2]),dim(database[[i]][[j]])[1])[13][,1]
      var.g2[j] <- mes(mean(database[[i]][[j]][,3]),sd(database[[i]][[j]][,3]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,4]),sd(database[[i]][[j]][,4]),dim(database[[i]][[j]])[1])[13][,1]
      covar.g1g2[j] <- correlacion*sqrt(var.g1[j])*sqrt(var.g2[j])
      input <- cbind(g1,g2,var.g1,covar.g1g2,var.g2)
    }
    datos_meta[[i]] <- as.data.frame(input)
  }
  return(datos_meta)
}

##############################################################################################################################
##########################FUNCION 5: METAANALISIS MULTIVARIADO DE SIMULACION (MODIFICADA)#####################################
##############################################################################################################################

###############
# FUNCION 5.1 #
###############

#datos <- g_hedges_simulacion2(database,correlacion=0)

meta_multi_simulacion2 <- function(data,metodo="reml"){
  library(mvmeta)
  meta_multi <- mvmeta(cbind(g1,g2),as.data.frame(data)[3:5],data=data,method=metodo)
  coeficientes <- meta_multi$coefficients
  coeficientes_inferencia <- summary(meta_multi)$coefficients
  coeficientes_var_cov <- summary(meta_multi)$corRandom
  output <- list(coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  return(output)
}

###############
# FUNCION 5.2 #
###############

meta_multi_simulacion2_bis <- function(n.estudios=5,n.vars=4,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion=0,
                                       metodo="reml"){
  library(mvmeta)
  datos <- g_hedges_simulacion2(n.estudios=5,n.vars=4,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion=0)
  data <- as.data.frame(datos)
  meta_multi <- mvmeta(cbind(g1,g2),as.data.frame(data)[3:5],data=data,method=metodo)
  coeficientes <- meta_multi$coefficients
  coeficientes_inferencia <- summary(meta_multi)$coefficients
  coeficientes_var_cov <- summary(meta_multi)$corRandom
  output <- list(coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  return(output)  
}
##############################################################################################################################
##########################FUNCION 6: METAANALISIS MULTIVARIADO DE SIMULACION (UNO POR REPLICACION)############################
##############################################################################################################################

###############
# FUNCION 6.1 #
###############

#datos <- g_hedges_simulacion_2vars(database,correlacion=0)

meta_multi_simulacion <- function(data,metodo="reml"){
  library(mvmeta)
  meta_multi <- vector(mode="list",length=length(data))
  for(i in 1:length(data)){
    datos <- data[[i]]
    meta_resultados <- mvmeta(cbind(g1,g2),S=as.data.frame(data)[3:5],data=datos,method=metodo)
    coeficientes <- meta_resultados$coefficients
    coeficientes_inferencia <- summary(meta_resultados)$coefficients
    coeficientes_var_cov <- summary(meta_resultados)$corRandom
    meta_multi[[i]] <- list(dat=data[[i]],coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  }
  return(meta_multi)  
}

###############
# FUNCION 6.2 #
###############

meta_multi_simulacion_bis <- function(n.estudios=5,n.vars=4,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion=0,
                                         metodo="reml"){
  library(mvmeta)
  datos <- g_hedges_simulacion(n.estudios,n.vars,tamano.muestra,semilla,replicaciones,correlacion)
  meta_multi <- vector(mode="list",length=replicaciones)
  for(i in 1:replicaciones){
    data <- as.data.frame(datos[[i]])
    meta_resultados <- mvmeta(cbind(g1,g2),S=as.data.frame(data)[3:5],data=data,method=metodo)
    coeficientes <- meta_resultados$coefficients
    coeficientes_inferencia <- summary(meta_resultados)$coefficients
    coeficientes_var_cov <- summary(meta_resultados)$corRandom
    meta_multi[[i]] <- list(dat=datos[[i]],coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  }
  return(meta_multi)  
  
}