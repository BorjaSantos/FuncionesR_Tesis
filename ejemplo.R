# EJEMPLO 1 DE METAANALISIS MULTIVARIANTE CON TRES ESTIMADORES DEL EFECTO.

# ESTIMADOR DEL EFECTO VARIABLE 1
smd1 <- rnorm(6,0.1) 
var1 <- runif(6,5,20)

# ESTIMADOR DEL EFECTOS VARIABLE 2
smd2 <- rnorm(6,0.1) 
var2 <- runif(6,5,20)

# ESTIMADOR DEL EFECTO VARIABLE 3
smd3 <- rnorm(6,0.1) 
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
library(miscTools)
matrices.correlacion <- corr.matrix.gen.3vars(5,c(rbeta(5,0.5,0.6),runif(5,0.2,1),runif(5,-1,1)),2014)
lapply(matrices.correlacion,det)
# Matrices de correlacion corregidas:
matrices.correlacion2 <- transformar.matrices.3vars(matrices.correlacion) 
lapply(matrices.correlacion2,det)

# Matriz de correlacion imputada:
matriz.correlacion <- matriz.imputada(matrices.correlacion2,3)$matriz.imputada

# Matrices de varianzas-covarianzas:
matrices.metaanalisis <- matrices.varcovar(matriz.correlacion,sd)


# metaanalisis multivariante:
library(mvmeta)
meta_prueba <- mvmeta(cbind(var1,var2,var3),matrices.metaanalisis)

# EJEMPLO 2 DE METAANALISIS MULTIVARIANTE CON TRES ESTIMADORES DEL EFECTO.

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
# Matrices de correlacion corregidas:
matrices.correlacion2 <- transformar.matrices.3vars(matrices.correlacion) 
lapply(matrices.correlacion2,det)

# Matriz de correlacion imputada:
matriz.correlacion <- matriz.imputada(matrices.correlacion2,3)$matriz.imputada

# Matrices de varianzas-covarianzas:
matrices.metaanalisis <- matrices.varcovar(matriz.correlacion,sd)


# metaanalisis multivariante:
library(mvmeta)
meta_prueba <- mvmeta(cbind(var1,var2,var3),matrices.metaanalisis)

# COMPROBACION DE LAS FUNCIONES PROGRAMADAS MEDIANTE SU APLICACION A UN METAANALISIS BIVARIADO PUBLICADO.

# Base de datos:
datos <- berkey98

# Calculo de las matrices de varianzas-covarianzas VCOV = S*R*S

# Matrices de desviaciones estandar
sd1 <- diag(c(sqrt(datos[1,5]),sqrt(datos[1,7])))
sd2 <- diag(c(sqrt(datos[2,5]),sqrt(datos[2,7])))
sd3 <- diag(c(sqrt(datos[3,5]),sqrt(datos[3,7])))
sd4 <- diag(c(sqrt(datos[4,5]),sqrt(datos[4,7])))
sd5 <- diag(c(sqrt(datos[5,5]),sqrt(datos[5,7])))
sd <- list(sd1,sd2,sd3,sd4,sd5)

# Matrices de correlacion:
matrices.correlacion <- corr.matrix.gen.2vars(5,c(runif(5,0.2,1)),2014)
lapply(matrices.correlacion,det)
# Matrices de correlacion corregidas:
matrices.correlacion2 <- transformar.matrices.2vars(matrices.correlacion) 
lapply(matrices.correlacion2,det)

# Matriz de correlacion imputada:
matriz.correlacion <- matriz.imputada(matrices.correlacion2,2)$matriz.imputada

# Matrices de varianzas-covarianzas:
matrices.metaanalisis <- matrices.varcovar(matriz.correlacion,sd)
meta_prueba <- mvmeta(cbind(PD,AL),matrices.metaanalisis,data=datos)
meta_prueba