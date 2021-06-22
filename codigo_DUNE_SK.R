
#Con este codigo ajustaremos los datos de las oscilaciones que generaremos para
#los experimentos Super Kamiokande y DUNE

rm(list=ls()) #Eliminamos las variables que tengamos
#Cargamos una serie de librerias y fuentes que nos seran de utilidad
library(latex2exp)
library(extrafont) 
font_import(pattern = "lmroman*") 
#En este algoritmo resumiré la mayoría de aspectos vistos en este algoritmo de cadenas 
loadfonts(device = "win")
par(family = "LM Roman 12")
windowsFonts("LM Roman 12"= windowsFont("LM Roman 12"))

#Ajuste de datos de Super-Kamiokande para neutrinos atmosfericos

#Datos del experimento 
theta = pi/4
E = 10 #GeV
deltam = 2.5*10^(-3) #Electronvolts al cuadrado
R = 6371 #Radio de la Tierra (km)
Ratm <- R+15 #Radio de los neutrinos atmosféricos (km)
Rdet <- R-1 #Profundidad del detector (km)

#Funcion de las oscilaciones
oscil <- function(L, angulo, masacuad){
  y = 1- (sin(2*angulo))^2*(sin(1.27*(masacuad*L)/(E)))^2
  return(y)
}

#Longitud recorrida en funcion del coseno del angulo de mezcla
L <-function(c){ #Para generar los puntos como argumento el coseno del ángulo
  y = sqrt(Ratm^2-Rdet^2*(1-c^2)) - Rdet*c
  return(y)
}

#Generamos los valores de la funcion como si fueran nuestros datos experimentales
#Se generan valores aleatorios que siguen una distribucion gaussiana
#centrada en el valor de la funcion calculado con los datos del experimento SK
#y de anchura 0.05, generamos 10 datos
valoresfun <-vector()
c2 = seq(-1, 1, length.out=10)
for(k in 1:10){
  valoresfun<- c(valoresfun,rnorm(1,oscil(L(c2[k]),theta,deltam),0.05))
}

testf <- function(xx){ #Vamos a coger nuestra función de oscilaciones
  #x1 angulo y x2 incremento de masa al cuadrado
  x1 <- xx[1] 
  x2 <- xx[2]

#Creamos un cector con los valores de L
valoresL<-vector()
y = 0;
c2 = seq(-1, 1, length.out=10)#puntos equiespaciados ente -1 y 1 del coseno
for (k in 1:10){
  c = c2[k] 
valoresL = c(valoresL, L(c)) #Almacena los valores de L
y <- y +(valoresfun[k] - oscil(L(c),x1,x2))^2 #Ya hemos generado nuestra función
#que tendremos que minimizar es la suma cuadratica de las distancias entre nuestros
#puntos generados y la funcion dependiendo de los parametros a ajustar que en este
#caso son la masa y el angulo
}
return(y)
}

#A continuacion el algoritmo  de MCMC simulated annealing, como vemos los 
#parametros de entrada son la funcion que se quiere minimizar (func), el valor
#inicial (s0), el numero de iteraciones (niter) y el paso (step). En este caso
#la anchura de la gaussiana para generar el siguiente punto esta dentro de esta
#funcion pero tambien es un parametro de entrada (que fijamos) del algoritmo
#La explicacion del algoritmo es la misma que en el caso del ajuste al minimo
#de las gaussianas

simulated_annealing <- function(func, s0, niter , step) {
  

  s_b <- s_c <- s_n <- s0 
  f_b <- f_c <- f_n <- func(s_n) 
 
  track1<-vector() 
  track2<-vector()
  track1<-c(track1,s0[1])
  track2<-c(track2,s0[2])

 
  for (k in 1:niter) {     
    Temp <- (1 - step)^k #Temperatura geometrica
    #Para generar los nuevos valores en este caso, se hacen de los dos parametros
    #por separado pues tienen valores muy distintos
    s_n[1] <- rnorm(1, s_c[1], .004) 
    #Sabemos (y en los experimentos reales tambien) que nuestro valor esta dentro de un
    #intervalo, impondremos entonces que no se generen valores fuera de ese intervalo
    while (s_n[1] < 0.6 || s_n[1] > 0.9){
      s_n[1]<-rnorm(1,s_c[1], .004)
    }
    #Lo mismo para el segundo parametro
    s_n[2] <- rnorm(1, s_c[2], 0.00004) 
    while (s_n[2]<0.0015 || s_n[2]>0.0035){
      s_n[2]<-rnorm(1,s_c[2], 0.00004)
    }
    f_n <- testf(s_n) #Calcula el nuevo valor
   
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

repetir<-function(ss){
  sol<-simulated_annealing(testf, ss ,30000,0.07);
  y<-rbind(sol$values1,sol$values2)
  ybest = sol$best_state
  return (list(values=y, best=ybest))}
 
#Inicializamos el mejor valor a 5 pues estamos buscando un minimo que es 0
bestvalue = 5;

#De un gran número de cadenas que generamos seleccionamos el mejor valor, es decir el
#valor minimo, lo hacemos asi para tener una poblacion estadistica significativa
#pues el valor que se obtiene depende mucho del punto inicial

for(k in 1:300){
  y<-repetir(c(runif(1,0.6,0.9),runif(1,0.0015,0.0035)))$best;
    if(testf(y) < bestvalue ){ #Esto nos devolverá el mejor mínimo de todas
      ybest<-y
      bestvalue<-testf(y)
  }
}

print(ybest)
print(bestvalue)
#Obtenemos una buena estimación de nuestros parámetros

#Representamos ahora los 10 datos y la función a la que se han ajustado.
c1 = seq(-1, 1, length.out=1000);
c2 = seq(-1, 1, length.out=10);
par(mfrow=c(1,2))
plot(c1,oscil(L(c1),ybest[1],ybest[2]),type='l', xlab=expression(italic('x')),ylab=expression(italic('y')), cex.axis = 1)
points(c2, valoresfun)


#Hacemos lo mismo que para SK pero ahora con DUNE
rm(list=ls())

#En este caso tenemos 5 parametros a ajustar, tres angulos y dos diferencias de 
#masa al cuadrado
th12 = 0.5903
th23 = 0.866
th13=0.15
m21=7.39*10^(-5)
m31=2.451*10^(-3)
dens=2.848 #La densidad es un dato del experimento DUNE
a = (1/3500)*(dens/3)
L = 1285 #Distancia recorrida en DUNE

#Funcion de las oscilaciones a primer orden
oscil <- function(E, the12,the23,the13,mas21,mas31){
  y = ((sin(the23))^2*(sin(2*the13))^2*(sin(1.27*mas31*L/E-a*L))^2*(1.27*mas31*L/E)^2)/(1.27*mas31*L/E-a*L)^2 + sin(2*the23)*sin(2*the13)*sin(2*the12)*sin(1.27*mas31*L/E-a*L)*(1.27*mas31*L/E)*sin(a*L)*(1.27*mas21*L/E)*cos(1.27*mas31*L/E)/(a*L*(1.27*mas31*L/E-a*L))+(cos(the23))^2*(sin(2*the12))^2*(sin(a*L))^2*(1.27*mas21*L/E)^2/(a*L)^2
  return(y)
}

#Generamos los valores igual que como lo hemos hecho en SK. En este caso generamos 20
E2 <- seq(0.5,8,length.out=20);
valoresfun <-vector()
for (k in 1:20){
  valoresfun = c(valoresfun, rnorm(1, oscil(E2[k],th12,th23,th13,m21,m31),0.005)) #Almacena valores funcion
}


testf <- function(xx){ #Vamos a crear la funcion a minimizar al igual que con SK
  x1 <- xx[1] 
  x2 <- xx[2]
  x3 <- xx[3]
  x4 <- xx[4]
  x5 <- xx[5]
  
  y = 0;
  E2 <- seq(0.5,8,length.out=20);
  for (k in 1:20){
    E =  E2[k]
    y <- y +(valoresfun[k] - oscil(E,x1,x2,x3,x4,x5))^2 #Ya hemos generado nuestra función
  }
  return(y)
}

simulated_annealing <- function(func, s0, niter , step) {
 
  s_b <- s_c <- s_n <- s0 
  f_b <- f_c <- f_n <- func(s_n) 
  
  track1<-vector() 
  track2<-vector()
  track1<-c(track1,s0[1])
  track2<-c(track2,s0[2])
  
  
  for (k in 1:niter) {     
    Temp <- (1 - step)^k 
    #Hacemos lo mismo que en el caso de SK, acotamos los valores donde puede estar
    s_n[1] <- rnorm(1, s_c[1], .004) 
    while (s_n[1] < 0.50 || s_n[1] > 0.65){
      s_n[1]<-rnorm(1,s_c[1], .004)
    }
    
    s_n[2] <- rnorm(1, s_c[2], 0.004) 
    while (s_n[2]<0.80 || s_n[2]>0.95){
      s_n[2]<-rnorm(1,s_c[2], 0.004)
    }
    s_n[3] <- rnorm(1, s_c[3], 0.004) 
    while (s_n[3]<0.10 || s_n[3]>0.20){
      s_n[3]<-rnorm(1,s_c[3], 0.004)
    }
    s_n[4] <- rnorm(1, s_c[4], 0.0000004) 
    while (s_n[4]<6*10^(-5) || s_n[4]>8*10^(-5)){
      s_n[4]<-rnorm(1,s_c[4], 0.0000004)
    }
    s_n[5] <- rnorm(1, s_c[5], 0.00004) 
    while (s_n[5]<0.0020 || s_n[5]>0.0030){
      s_n[5]<-rnorm(1,s_c[5], 0.00004)
    }
    
    f_n <- testf(s_n)
    
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


repetir<-function(ss){
  sol<-simulated_annealing(testf, ss ,50000,0.07);
  y<-rbind(sol$values1,sol$values2)
  ybest = sol$best_state
  return (list(values=y, best=ybest))}

bestvalue = 5;

#De un gran número de cadenas que generamos seleccionamos el mejor valor

for(k in 1:5){
  y<-repetir(c(runif(1,0.50,0.65),runif(1,0.80,0.95),runif(1,0.10,0.20),runif(1,6*10^(-5),8*10^(-5)),runif(1,0.0020,0.0030)))$best;
  if(testf(y) < bestvalue ){ #Esto nos devolverá el mejor mínimo de todas
    ybest<-y
    bestvalue<-testf(y)
  }
}

E1 <- seq(0.5, 8, length.out=1000);
E2 <- seq(0.5,8,length.out=20);
plot(E1,oscil(E1,ybest[1],ybest[2],ybest[3],ybest[4],ybest[5]),type='l',  xlab=expression(italic('x')),ylab=expression(italic('y')), cex.axis=1 )
points(E2, valoresfun)



print(ybest)
print(bestvalue)




