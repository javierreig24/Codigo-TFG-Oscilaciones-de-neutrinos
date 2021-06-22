#En este codigo se presentara el método MCMC para encontrar el mínimo principal
#de una doble gaussiana, así como de una gaussiana simple. Se estudiaran ademas
#los diferentes parámetros de los que depende el algoritmo

rm(list=ls()) #Inicialmente borramos todas las variables que esten en la memoria

#Cargamos algunas librerias y fuentes que nos seran utiles
library(ggplot2)
library(latex2exp)
library(extrafont) 
font_import(pattern = "lmroman*") 
#En este algoritmo resumiré la mayoría de aspectos vistos en este algoritmo de cadenas 
loadfonts(device = "win")
par(family = "LM Roman 12")
windowsFonts("LM Roman 12"= windowsFont("LM Roman 12"))

testf <- function(xx) #Test para una funcion de dos variables
{x1 <- xx[1] #Cogemos la primera variable
x2 <- xx[2] #Cogemos la segunda variable

#Mi función es la suma de dos distribuciones normales
make.mvn <- function(mean, vcv) {
  logdet <- as.numeric(determinant(vcv, TRUE)$modulus)
  tmp <- length(mean) * log(2 * pi) + logdet
  vcv.i <- solve(vcv)
  
  function(x) {
    dx <- x - mean
    exp(-(tmp + rowSums((dx %*% vcv.i) * dx))/2)
  }
}
mu1 = c(-1, 1) #media de una gaussiana
mu2 = c(2, -2) #media de la otra gaussiana
vcv1 = matrix(c(1, .25, .25, 1.5), 2, 2) #matriz de distribucion normal en 2D
vcv2 = matrix(c(2, -.5, -.5, 2), 2, 2)
f1 = make.mvn(mu1, vcv1)
f2 = make.mvn(mu2, vcv2)
y= -(f1(xx) + f2(xx)) #Suma  de las dos gaussianas negativas
return(y)
}

#A continuacion el algoritmo  de MCMC simulated annealing, como vemos los 
#parametros de entrada son la funcion que se quiere minimizar (func), el valor
#inicial (s0), el numero de iteraciones (niter) y el paso (step) y la anchura
#de la gaussiana para generar el siguiente punto, sigma
simulated_annealing <- function(func, s0, niter , step, sigma ) {
  
  # Inicializamos las distintas variables del algoritmo
  ## s_b y f_b son el punto y el valor de la funcion mejores
  ## s_c y f_c son el punto y el valor de la funcion actuales
  ## s_n y f_n son el punto y el valor y la funcion nuevos
  s_b <- s_c <- s_n <- s0 #Inicializamos todos al valor inicial
  f_b <- f_c <- f_n <- func(s_n) 
  #vector para acumular resultados
  track1<-vector() 
  track2<-vector()
  track1<-c(track1,s0[1])#c función que combina los argumentos
  track2<-c(track2,s0[2])
  
  #Empieza el cuerpo del algoritmo que se repite hasta las iteraciones, niter
  for (k in 1:niter) {     
    Temp <- (1 - step)^k #Temperatura se va haciendo cada vez mas sensible, temperatura geometrica
    #Temp <- 1/log(1 + k) Temperatura logaritmica
    # Consideramos unos vecinos a nuestro punto aleatorios
    #El nuevo punto es un valor gaussiano aleatorio  centado en el actual, s_c 
    #con una anchura de 0.15 para ambas variables
    s_n <- rnorm(2, s_c, sigma) #El vecino nuevo es aleatorio, 2 valores por las dos variables, la media es el valor de antes y desviación estándar0.15 
    f_n <- func(s_n) #Calcula el nuevo valor
    # actualizamos el estado actua
    if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
      #Para que el nuevo valor pase a ser el actual, este ha de ser menor que 
      #el valor guardado como actual o que un numero aleatorio entre 0 y 1
      #sea menor que la exponencial exp(-(f_n - f_c) / Temp)
      s_c <- s_n
      f_c <- f_n
      track1<-c(track1,s_n[1]) #añade el de antes mas el nuevo
      track2<-c(track2,s_n[2])
    }
    #actualizamos el mejor estado, si el nuevo valor de la funcion es menos que
    #el que esta guardado como actual
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n 
    }
  }
  return(list(best_value = f_b, best_state = s_b, values1 = track1, values2 = track2))
}
#El algoritmo devuelve los diferentes valores de la cadena y los mejores valores

#A continuacion algunos de los analisis con este algoritmo. En primer lugar,
#representaremos el recorrido de las cadenas sobre las curvas de nivel de nuestra
#funcion, empezamos con parametros que dejan a la cadena recorrer poco espacio de 
#fases y son cada vez mas precisos. Funcion geometrica, anchura de la 
#gaussiana para el siguiente punto 0.15, 1200 iteraciones por cadena y un paso de 0.1

x = seq(-5, 6, length=71)
y = seq(-7, 6, length=71)
xy = expand.grid(x=x, y=y)
z = matrix(apply(as.matrix(xy), 1, testf), length(x), length(y))
par(mfrow=c(1,3))
image(x, y, z, las=1,xlab=expression(italic('x')),ylab=expression(italic('y')), main=expression('(a)'),cex.main=2, cex.lab=2,  cex.axis=1.2 )
contour(x, y, z, add=TRUE)

#Funcion para repetir las cadenas
repetir<-function(ss){
  sol<-simulated_annealing(testf, ss ,1200,0.1,0.15);
  y<-rbind(sol$values1,sol$values2)
  ybest = sol$best_state
  return (list(values=y, best=ybest))}

#Representamos 10 cadenas
replicate(10,{y<-repetir(runif(2,-5,5))$values;lines(y[1,],y[2,], col="#00000088",lwd=3)})

#A continuacion contamos las cadenas que se van de los dos minimos, las que van al
#minimo principal, las que van al secundario y de entre todas las cadenas,
#sacamos el valor de la función minimo que sera el mejor
cadenasquesevan=0;
minprincipal=0;
minsecundario=0;
bestvalue = 0;
ybest = c(0,0);

#Sabiendo donde estan los minimos podemos contar las cadenas que se van y las que no
for(k in 1:1000){
  y<-repetir(runif(2,-5,5))$best;
  if(( -1.1 < y[1]) && (y[1] < -.9)){
    minprincipal<-minprincipal+1 #Aqui contamos las cadenas que van al principal
    if(testf(y)<bestvalue){ #Esto nos devolvera el mejor mínimo de todas las cadenas
      ybest<-y
     bestvalue<-testf(y)
    }
  }
  else{
    if(( 1.9 < y[1]) && (y[1] < 2.1)){
      minsecundario <-minsecundario+1 #Las cadenas que van al secundario
    }else{
      cadenasquesevan = cadenasquesevan+1 #Las cadenas que se pierden
    }
  }
}

#Obtenemos los diferentes valores para analizar el resultado
print(cadenasquesevan)
print(minprincipal)
print(minsecundario)
print(ybest)
print(bestvalue)


#Vemos de estos resultados que muy pocas cadenas se van, y la mayoria van al 
#minimo principal, aunque un gran numero de las cadenas van tambien al secundario
#esto es asi porque no hemos dejado a la cadena recorrer todo el espacio de 
#fases y es muy sensible al valor inicial, pero es muy precisa

#Hacemos ahora dejando recorrer más espacio a las cadenas
#Aumentamos las iteraciones, con la temperatura logaritmica, aumentamos anchura

#Nuevo algoritmo con la temperatura logaritmica
simulated_annealing <- function(func, s0, niter , step, sigma ) {
  
  s_b <- s_c <- s_n <- s0 
  f_b <- f_c <- f_n <- func(s_n) 

  track1<-vector() 
  track2<-vector()
  track1<-c(track1,s0[1])
  track2<-c(track2,s0[2])
  
  for (k in 1:niter) {     
    Temp <- 1/log(1 + k) #Temperatura logaritmica
    s_n <- rnorm(2, s_c, sigma) 
    f_n <- func(s_n) 
    
    if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
      s_c <- s_n
      f_c <- f_n
      track1<-c(track1,s_n[1]) 
      track2<-c(track2,s_n[2])
    }
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n 
    }
  }
  return(list(best_value = f_b, best_state = s_b, values1 = track1, values2 = track2))
}


x = seq(-5, 6, length=71)
y = seq(-7, 6, length=71)
xy = expand.grid(x=x, y=y)
z = matrix(apply(as.matrix(xy), 1, testf), length(x), length(y))

image(x, y, z, las=1, xlab=expression(italic('x')),ylab=expression(italic('y')),main=expression('(b)'),cex.main=2, cex.lab=2,  cex.axis=1.2)
contour(x, y, z, add=TRUE)

repetir<-function(ss){
  sol<-simulated_annealing(testf, ss ,2000,1.2,0.8);
  y<-rbind(sol$values1,sol$values2)
  ybest = sol$best_state
  return (list(values=y, best=ybest))}

replicate(10,{y<-repetir(runif(2,-5,5))$values;lines(y[1,],y[2,], col="#00000088",lwd=1)})

cadenasquesevan=0;
minprincipal=0;
minsecundario=0;
bestvalue = 0;
ybest = c(0,0);

for(k in 1:1000){
  y<-repetir(runif(2,-5,5))$best;
  if(( -1.1 < y[1]) && (y[1] < -.9)){
    minprincipal<-minprincipal+1 #Aquí contamos las cadenas que van al principal
    if(testf(y)<bestvalue){ #Esto nos devolverá el mejor mínimo de todas las cadenas
      ybest<-y
      bestvalue<-testf(y)
    }
  }
  else{
    if(( 1.9 < y[1]) && (y[1] < 2.1)){
      minsecundario <-minsecundario+1 #Las cadenas que van al secundario
    }else{
      cadenasquesevan = cadenasquesevan+1 #Las cadenas que se pierden
    }
  }
}

print(cadenasquesevan)
print(minprincipal)
print(minsecundario)
print(ybest)
print(bestvalue)
#Si hacemos esto muchísimas cadenas se nos van y no convergen bien, aunque haya
#más proporción de las que se van al mínimo principal, la mayoria no van a 
#ninguno de los dos minimos

#Aplicamos ahora este algoritmo a una sola gaussiana
make.mvn <- function(mean, vcv) {
  logdet <- as.numeric(determinant(vcv, TRUE)$modulus)
  tmp <- length(mean) * log(2 * pi) + logdet
  vcv.i <- solve(vcv)
  
  function(x) {
    dx <- x - mean
    exp(-(tmp + rowSums((dx %*% vcv.i) * dx))/2)
  }
}

testf <- function(xx) #Test para una función de dos variables
{x1 <- xx[1] #Cogemos la primera variable
x2 <- xx[2] #Cogemos la segunda variable
mu1 = c(-1, 1) 
vcv1 = matrix(c(1, .25, .25, 1.5), 2, 2) 
f1 = make.mvn(mu1, vcv1)
y= -f1(xx)
return(y)
}


simulated_annealing <- function(func, s0, niter , step, sigma ) {

  s_b <- s_c <- s_n <- s0 
  f_b <- f_c <- f_n <- func(s_n) 
  track1<-vector() 
  track2<-vector()
  track1<-c(track1,s0[1])
  track2<-c(track2,s0[2])
  for (k in 1:niter) {     
    Temp <- (1 - step)^k
    s_n <- rnorm(2, s_c, sigma)
    f_n <- func(s_n) 
    if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
      s_c <- s_n
      f_c <- f_n
      track1<-c(track1,s_n[1]) 
      track2<-c(track2,s_n[2])
    }
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n 
    }
  }
  return(list(best_value = f_b, best_state = s_b, values1 = track1, values2 = track2))
}


x = seq(-5, 6, length=71)
y = seq(-7, 6, length=71)
xy = expand.grid(x=x, y=y)
z = matrix(apply(as.matrix(xy), 1, testf), length(x), length(y))
image(x, y, z, las=1,xlab=expression(italic('x')),ylab=expression(italic('y')),main=expression('(c)'),cex.main=2, cex.lab=2,  cex.axis=1.2)
contour(x, y, z, add=TRUE)

repetir<-function(ss){
  sol<-simulated_annealing(testf, ss ,1200,0.1, 0.15);
  y<-rbind(sol$values1,sol$values2)
  ybest = sol$best_state
  return (list(values=y, best=ybest))}

replicate(10,{y<-repetir(runif(2,-5,5))$values;lines(y[1,],y[2,], col="#00000088",lwd=3)})

cadenasquesevan=0;
minprincipal=0;
minsecundario=0;
bestvalue = 0;
ybest = c(0,0);

for(k in 1:1000){
  y<-repetir(runif(2,-5,5))$best;
  if(( -1.1 < y[1]) && (y[1] < -.9)){
    minprincipal<-minprincipal+1 
    if(testf(y)<bestvalue){ 
      ybest<-y
      bestvalue<-testf(y)
    }
  }
  else{
    cadenasquesevan = cadenasquesevan+1
    
  }
}

print(cadenasquesevan)
print(minprincipal)
print(minsecundario)
print(ybest)
print(bestvalue)

#Concluimos que al poner una sola gaussiana la inmensa mayoría convergen al minimo
#lo que corrobora la potencia del metodo


