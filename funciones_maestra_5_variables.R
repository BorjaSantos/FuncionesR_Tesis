################################################################################################################################
#################################FUNCION 1: SIMULACION DE ESTUDIOS CON DATOS INDIVIDUALES#######################################
################################################################################################################################

# PARAMTEROS IMPORTANTES EN LA SIMULACION:
# 1. Numero de estudios
# 2. Numero de variables
# 3. Tamaño de la muestra
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
# m3.trat / m3.trat: media grupo tratamiento / media grupo control para la variable 3.
# sd3.trat / sd3.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 3.
# a3 / b3: limites superior e inferior del rango de valores de la variable 3. 
# m4.trat / m4.trat: media grupo tratamiento / media grupo control para la variable 4.
# sd4.trat / sd4.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 4.
# a4 / b4: limites superior e inferior del rango de valores de la variable 4. 
# m5.trat / m5.trat: media grupo tratamiento / media grupo control para la variable 5.
# sd5.trat / sd5.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 5.
# a5 / b5: limites superior e inferior del rango de valores de la variable 5. 

# Libreri???as necesarias para la simulacion
library(truncdist) # De esta manera se truncaran los valores de las diferentes puntuaciones, para que no se simulen puntuaciones
# fuera de rango.
# Semilla para poder replicar los datos
set.seed(18052013)

# Funcion que va a simular los datos en el escenario que queramos

simulacion_datos <- function(n.estudios,m1.trat,sd1.trat,m1.ctrl,sd1.ctrl,a1,b1,m2.trat,sd2.trat,m2.ctrl,sd2.ctrl,a2,b2,
                             m3.trat,sd3.trat,a3,b3,m4,trat,sd4.trat,m4.ctrl,sd4.ctrl,a4,b4,m5.trat,sd5.trat,a5,b5,
                             tamano.muestra,semilla,
                             replicaciones){
  # Control de los paramtero de la funcion
  if(n.estudios > 50 || n.estudios < 5){
    stop("Numero de estudios incorrecto")
  }
  if((n.vars > 12 || n.vars < 4) && is.integer(n.vars/2)!= TRUE){
    stop("Numero de variables por estudio incorrecto")
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
    for(j in 1:n.estudios){
      Var1.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.trat,sd=sd1.trat),0) # Puntuaciones variable 1 en tratamiento
      Var1.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.ctrl,sd=sd1.ctrl),0) # Puntuaciones variable 1 en control
      Var2.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.trat,sd=sd2.trat),0) # Puntuaciones variable 2 tratamiento
      Var2.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.ctrl,sd=sd2.ctrl),0) # Puntuaciones variable 2 control
      Var3.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a3,b=b3,mean=m3.trat,sd=sd3.trat),0) # Puntuaciones variable 3 tratamiento
      Var3.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a3,b=b3,mean=m3.ctrl,sd=sd3.ctrl),0) # Puntuaciones variable 3 control
      Var4.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a4,b=b4,mean=m4.trat,sd=sd4.trat),0) # Puntuaciones variable 4 tratamiento
      Var4.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a4,b=b4,mean=m4.ctrl,sd=sd4.ctrl),0) # Puntuaciones variable 4 control
      Var5.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a6,b=b5,mean=m5.trat,sd=sd5.trat),0) # Puntuaciones variable 5 tratamiento
      Var5.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a5,b=b5,mean=m5.ctrl,sd=sd5.ctrl),0) # Puntuaciones variable 5 control
      database[[i]][[j]] <- as.data.frame(cbind(Var1.Trat,Var1.Ctrl,Var2.Trat,Var2.Ctrl,Var3.Trat,Var3.Ctrl,Var4.Trat,
                                                Var4.Ctrl,Var5.Trat,Var5.Ctrl))
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
# m1.trat / m1.trat: media grupo tratamiento / media grupo control para la variable 1.
# sd1.trat / sd1.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 1.
# a1 / b1: limites superior e inferior del rango de valores de la variable 1.
# m2.trat / m2.trat: media grupo tratamiento / media grupo control para la variable 2.
# sd2.trat / sd2.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 2.
# a2 / b2: limites superior e inferior del rango de valores de la variable 2.
# m3.trat / m3.trat: media grupo tratamiento / media grupo control para la variable 3.
# sd3.trat / sd3.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 3.
# a3 / b3: limites superior e inferior del rango de valores de la variable 3. 
# m4.trat / m4.trat: media grupo tratamiento / media grupo control para la variable 4.
# sd4.trat / sd4.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 4.
# a4 / b4: limites superior e inferior del rango de valores de la variable 4. 
# m5.trat / m5.trat: media grupo tratamiento / media grupo control para la variable 5.
# sd5.trat / sd5.trat: desviacion estandar grupo tratamiento / desviacion estandar grupo control para la variable 5.
# a5 / b5: limites superior e inferior del rango de valores de la variable 5. 

# Librerias necesarias para la simulacion
library(truncdist) # De esta manera se truncaran los valores de las diferentes puntuaciones, para que no se simulen puntuaciones
# fuera de rango.
# Semilla para poder replicar los datos
set.seed(18052013)

# Funcion que va a simular los datos en el escenario que queramos

simulacion_datos2 <- function(n.estudios,m1.trat,sd1.trat,m1.ctrl,sd1.ctrl,a1,b1,m2.trat,sd2.trat,m2.ctrl,sd2.ctrl,a2,b2,
                              m3.trat,sd3.trat,a3,b3,m4,trat,sd4.trat,m4.ctrl,sd4.ctrl,a4,b4,m5.trat,sd5.trat,a5,b5,
                              tamano.muestra,semilla,
                              replicaciones){
  # Control de los paramtero de la funcion
  if(n.estudios > 50 || n.estudios < 5){
    stop("Numero de estudios incorrecto")
  }
  if((n.vars > 12 || n.vars < 4) && is.integer(n.vars/2)!= TRUE){
    stop("Numero de variables por estudio incorrecto")
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
    Var1.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.trat,sd=sd1.trat),0) # Puntuaciones variable 1 en tratamiento
    Var1.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a1,b=b1,mean=m1.ctrl,sd=sd1.ctrl),0) # Puntuaciones variable 1 en control
    Var2.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.trat,sd=sd2.trat),0) # Puntuaciones variable 2 tratamiento
    Var2.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a2,b=b2,mean=m2.ctrl,sd=sd2.ctrl),0) # Puntuaciones variable 2 control
    Var3.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a3,b=b3,mean=m3.trat,sd=sd3.trat),0) # Puntuaciones variable 3 tratamiento
    Var3.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a3,b=b3,mean=m3.ctrl,sd=sd3.ctrl),0) # Puntuaciones variable 3 control
    Var4.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a4,b=b4,mean=m4.trat,sd=sd4.trat),0) # Puntuaciones variable 4 tratamiento
    Var4.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a4,b=b4,mean=m4.ctrl,sd=sd4.ctrl),0) # Puntuaciones variable 4 control
    Var5.Trat <- round(rtrunc(tamano.muestra,spec="norm",a=a6,b=b5,mean=m5.trat,sd=sd5.trat),0) # Puntuaciones variable 5 tratamiento
    Var5.Ctrl <- round(rtrunc(tamano.muestra,spec="norm",a=a5,b=b5,mean=m5.ctrl,sd=sd5.ctrl),0) # Puntuaciones variable 5 control
    database[[i]] <- as.data.frame(cbind(Var1.Trat,Var1.Ctrl,Var2.Trat,Var2.Ctrl,Var3.Trat,Var3.Ctrl,Var4.Trat,
                                         Var4.Ctrl,Var5.Trat,Var5.Ctrl))
    semilla <- semilla + (n.estudios) 
  }
 return(database)
}

##############################################################################################################################
####################FUNCION 3: g DE HEDGES' DE UNA SIMULACION CON REPLICACIONES (SIMPLIFICADA)################################
##############################################################################################################################

g_hedges_simulacion2_5vars <- function(database,correlacion12=0,
                                       correlacion13=0,correlacion14=0,correlacion15=0,correlacion23=0,correlacion24=0,
                                       correlacion25=0,correlacion34=0,correlacion35=0,correlacion45=0){
  library(compute.es)
  datos_meta <- diag(0,length(database),20)
  for(i in 1:replicaciones){
    g1 <- mes(mean(database[[i]][,1]),sd(database[[i]][,1]),dim(database[[i]])[1],mean(database[[i]][,2]),sd(database[[i]][,2]),dim(database[[i]])[1])[12][,1]
    g2 <- mes(mean(database[[i]][,3]),sd(database[[i]][,3]),dim(database[[i]])[1],mean(database[[i]][,4]),sd(database[[i]][,4]),dim(database[[i]])[1])[12][,1]
    g3 <- mes(mean(database[[i]][,5]),sd(database[[i]][,5]),dim(database[[i]])[1],mean(database[[i]][,6]),sd(database[[i]][,6]),dim(database[[i]])[1])[12][,1]
    g4 <- mes(mean(database[[i]][,7]),sd(database[[i]][,7]),dim(database[[i]])[1],mean(database[[i]][,8]),sd(database[[i]][,8]),dim(database[[i]])[1])[12][,1]
    g5 <- mes(mean(database[[i]][,9]),sd(database[[i]][,9]),dim(database[[i]])[1],mean(database[[i]][,10]),sd(database[[i]][,10]),dim(database[[i]])[1])[12][,1]
    var.g1 <- mes(mean(database[[i]][,1]),sd(database[[i]][,1]),dim(database[[i]])[1],mean(database[[i]][,2]),sd(database[[i]][,2]),dim(database[[i]])[1])[13][,1]
    var.g2 <- mes(mean(database[[i]][,3]),sd(database[[i]][,3]),dim(database[[i]])[1],mean(database[[i]][,4]),sd(database[[i]][,4]),dim(database[[i]])[1])[13][,1]
    var.g3 <- mes(mean(database[[i]][,5]),sd(database[[i]][,5]),dim(database[[i]])[1],mean(database[[i]][,6]),sd(database[[i]][,6]),dim(database[[i]])[1])[13][,1]
    var.g4 <- mes(mean(database[[i]][,7]),sd(database[[i]][,7]),dim(database[[i]])[1],mean(database[[i]][,8]),sd(database[[i]][,8]),dim(database[[i]])[1])[13][,1]
    var.g5 <- mes(mean(database[[i]][,9]),sd(database[[i]][,9]),dim(database[[i]])[1],mean(database[[i]][,10]),sd(database[[i]][,10]),dim(database[[i]])[1])[13][,1]
    covar.g1g2 <- correlacion12*sqrt(var.g1)*sqrt(var.g2)
    covar.g1g3 <- correlacion13*sqrt(var.g1)*sqrt(var.g3)
    covar.g1g4 <- correlacion14*sqrt(var.g1)*sqrt(var.g4)
    covar.g1g5 <- correlacion15*sqrt(var.g1)*sqrt(var.g5)
    covar.g2g3 <- correlacion23*sqrt(var.g2)*sqrt(var.g3)
    covar.g2g4 <- correlacion24*sqrt(var.g2)*sqrt(var.g4)
    covar.g2g5 <- correlacion25*sqrt(var.g2)*sqrt(var.g5)
    covar.g3g4 <- correlacion34*sqrt(var.g3)*sqrt(var.g4)
    covar.g3g5 <- correlacion35*sqrt(var.g3)*sqrt(var.g5)
    covar.g4g5 <- correlacion45*sqrt(var.g4)*sqrt(var.g5)
    input <- c(g1,g2,g3,g4,g5,var.g1,covar.g1g2,covar.g1g3,covar.g1g4,covar.g1g5,var.g2,covar.g2g3,covar.g2g4,covar.g2g5,
               var.g3,covar.g3g4,covar.g3g5,var.g4,covar.g4g5,var.g5)
    datos_meta[i,] <- input
  }
  colnames(datos_meta) <- c("g1","g2","g3","g4","g5","var.g1","covar.g1g2","covar.g1g3","covar.g1g4","covar.g1g5","var.g2",
                            "covar.g2g3","covar.g2g4","covar.g2g5","var.g3","var.g3g4","var.g3g5","var.g4","covar.g4g5",
                            "var.g5")
  return(as.data.frame(datos_meta))
}

##############################################################################################################################
##########################FUNCION 4: g DE HEDGES' DE UNA SIMULACION CON REPLICACIONES#########################################
##############################################################################################################################

g_hedges_simulacion_5vars <- function(database,correlacion12=0,
                                      correlacion13=0,correlacion14=0,correlacion15=0,correlacion23=0,correlacion24=0,
                                      correlacion25=0,correlacion34=0,correlacion35=0,correlacion45=0){
  library(compute.es)
  datos_meta <- vector("list",length(database))
  g1 <- vector(mode="numeric",length=length(database[[1]]))
  g2 <- vector(mode="numeric",length=length(database[[1]]))
  g3 <- vector(mode="numeric",length=length(database[[1]]))
  g4 <- vector(mode="numeric",length=length(database[[1]]))
  g5 <- vector(mode="numeric",length=length(database[[1]]))
  var.g1 <- vector(mode="numeric",length=length(database[[1]]))
  var.g2 <- vector(mode="numeric",length=length(database[[1]]))
  var.g3 <- vector(mode="numeric",length=length(database[[1]]))
  var.g4 <- vector(mode="numeric",length=length(database[[1]]))
  var.g5 <- vector(mode="numeric",length=length(database[[1]]))
  covar.g1g2 <- vector(mode="numeric",length=length(database[[1]])) 
  covar.g1g3 <- vector(mode="numeric",length=length(database[[1]]))
  covar.g1g4 <- vector(mode="numeric",length=length(database[[1]])) 
  covar.g1g5 <- vector(mode="numeric",length=length(database[[1]])) 
  covar.g2g3 <- vector(mode="numeric",length=length(database[[1]]))
  covar.g2g4 <- vector(mode="numeric",length=length(database[[1]])) 
  covar.g2g5 <- vector(mode="numeric",length=length(database[[1]])) 
  covar.g3g4 <- vector(mode="numeric",length=length(database[[1]])) 
  covar.g3g5 <- vector(mode="numeric",length=length(database[[1]])) 
  covar.g4g5 <- vector(mode="numeric",length=length(database[[1]])) 
  for(i in 1:length(database)){
    for(j in 1:length(database[[i]])){
      g1[j] <- mes(mean(database[[i]][[j]][,1]),sd(database[[i]][[j]][,1]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,2]),sd(database[[i]][[j]][,2]),dim(database[[i]][[j]])[1])[12][,1]
      g2[j] <- mes(mean(database[[i]][[j]][,3]),sd(database[[i]][[j]][,3]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,4]),sd(database[[i]][[j]][,4]),dim(database[[i]][[j]])[1])[12][,1]
      g3[j] <- mes(mean(database[[i]][[j]][,5]),sd(database[[i]][[j]][,5]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,6]),sd(database[[i]][[j]][,6]),dim(database[[i]][[j]])[1])[12][,1]
      g4[j] <- mes(mean(database[[i]][[j]][,7]),sd(database[[i]][[j]][,7]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,8]),sd(database[[i]][[j]][,8]),dim(database[[i]][[j]])[1])[12][,1]
      g5[j] <- mes(mean(database[[i]][[j]][,9]),sd(database[[i]][[j]][,9]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,10]),sd(database[[i]][[j]][,10]),dim(database[[i]][[j]])[1])[12][,1]
      var.g1[j] <- mes(mean(database[[i]][[j]][,1]),sd(database[[i]][[j]][,1]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,2]),sd(database[[i]][[j]][,2]),dim(database[[i]][[j]])[1])[13][,1]
      var.g2[j] <- mes(mean(database[[i]][[j]][,3]),sd(database[[i]][[j]][,3]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,4]),sd(database[[i]][[j]][,4]),dim(database[[i]][[j]])[1])[13][,1]
      var.g3[j] <- mes(mean(database[[i]][[j]][,5]),sd(database[[i]][[j]][,5]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,6]),sd(database[[i]][[j]][,6]),dim(database[[i]][[j]])[1])[13][,1]
      var.g4[j] <- mes(mean(database[[i]][[j]][,7]),sd(database[[i]][[j]][,7]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,8]),sd(database[[i]][[j]][,8]),dim(database[[i]][[j]])[1])[13][,1]
      var.g5[j] <- mes(mean(database[[i]][[j]][,9]),sd(database[[i]][[j]][,9]),dim(database[[i]][[j]])[1],mean(database[[i]][[j]][,10]),sd(database[[i]][[j]][,10]),dim(database[[i]][[j]])[1])[13][,1]
      covar.g1g2[j] <- correlacion12*sqrt(var.g1[j])*sqrt(var.g2[j])
      covar.g1g3[j] <- correlacion13*sqrt(var.g1[j])*sqrt(var.g3[j])
      covar.g1g4[j] <- correlacion14*sqrt(var.g1[j])*sqrt(var.g4[j])
      covar.g1g5[j] <- correlacion15*sqrt(var.g1[j])*sqrt(var.g5[j])
      covar.g2g3[j] <- correlacion23*sqrt(var.g2[j])*sqrt(var.g3[j])
      covar.g2g4[j] <- correlacion24*sqrt(var.g2[j])*sqrt(var.g4[j])
      covar.g2g5[j] <- correlacion25*sqrt(var.g2[j])*sqrt(var.g5[j])
      covar.g3g4[j] <- correlacion34*sqrt(var.g3[j])*sqrt(var.g4[j])
      covar.g3g5[j] <- correlacion35*sqrt(var.g3[j])*sqrt(var.g5[j])
      covar.g4g5[j] <- correlacion45*sqrt(var.g4[j])*sqrt(var.g5[j])
      input <- cbind(g1,g2,g3,g4,g5,var.g1,covar.g1g2,covar.g1g3,covar.g1g4,covar.g1g5,var.g2,covar.g2g3,covar.g2g4,covar.g2g5,
                     var.g3,covar.g3g4,covar.g3g5,var.g4,covar.g4g5,var.g5)
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

# datos <- g_hedges_simulacion2_5vars(database,correlacion12=0,correlacion13=0,correlacion14=0,correlacion15=0,correlacion23=0,
#                                     correlacion24=0,correlacion25=0,correlacion34=0,correlacion35=0,correlacion45=0)

meta_multi_simulacion2_5vars <- function(data,metodo="reml"){
  library(mvmeta)
  meta_multi <- mvmeta(cbind(g1,g2,g3,g4,g5),as.data.frame(data)[6:20],data=data,method=metodo)
  coeficientes <- meta_multi$coefficients
  coeficientes_inferencia <- summary(meta_multi)$coefficients
  coeficientes_var_cov <- summary(meta_multi)$corRandom
  output <- list(coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  return(output)
}

###############
# FUNCION 5.2 #
###############

meta_multi_simulacion2_5vars_bis <- function(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
                                             correlacion13=0,correlacion14=0,correlacion15=0,correlacion23=0,correlacion24=0,
                                             correlacion25=0,correlacion34=0,correlacion35=0,correlacion45=0,metodo="reml"){
  library(mvmeta)
  datos <- g_hedges_simulacion2_5vars(n.estudios,tamano.muestra,semilla,replicaciones,correlacion12,
                                      correlacion13,correlacion14,correlacion23,correlacion24,correlacion34)
  data <- as.data.frame(datos)
  meta_multi <- mvmeta(cbind(g1,g2,g3,g4,g5),as.data.frame(data)[6:20],data=data,method=metodo)
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

# datos <- g_hedges_simulacion_5vars(database,correlacion12=0,correlacion13=0,correlacion14=0,correlacion15=0,correlacion23=0,
#                                    correlacion24=0,correlacion25=0,correlacion34=0,correlacion35=0,correlacion45=0)

meta_multi_simulacion_5vars <- function(data,metodo="reml",replicaciones=5){
  library(mvmeta)
  meta_multi <- vector(mode="list",length=replicaciones)
  for(i in 1:length(data)){
    datos <- data[[i]]
    meta_resultados <- mvmeta(cbind(g1,g2,g3,g4,g5),S=as.data.frame(data)[6:20],data=datos,method=metodo)
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

meta_multi_simulacion_5_vars_bis <- function(n.estudios=5,tamano.muestra=20,semilla=18052013,replicaciones=5,correlacion12=0,
                                             correlacion13=0,correlacion14=0,correlacion15=0,correlacion23=0,correlacion24=0,
                                             correlacion25=0,correlacion34=0,correlacion35=0,correlacion45=0,metodo="reml"){
  library(mvmeta)
  datos <- g_hedges_simulacion_5vars(n.estudios,tamano.muestra,semilla,replicaciones,correlacion12,correlacion13,
                                     correlacion14,correlacion15,correlacion23,correlacion24,correlacion25,correlacion34,
                                     correlacion35,correlacion45)
  meta_multi <- vector(mode="list",length=replicaciones)
  for(i in 1:replicaciones){
    data <- as.data.frame(datos[[i]])
    meta_resultados <- mvmeta(cbind(g1,g2,g3,g4,g5),S=as.data.frame(data)[6:20],data=data,method=metodo)
    coeficientes <- meta_resultados$coefficients
    coeficientes_inferencia <- summary(meta_resultados)$coefficients
    coeficientes_var_cov <- summary(meta_resultados)$corRandom
    meta_multi[[i]] <- list(dat=datos[[i]],coefs=coeficientes,inferencia=coeficientes_inferencia,var_cov=coeficientes_var_cov)
  }
  return(meta_multi)  
}