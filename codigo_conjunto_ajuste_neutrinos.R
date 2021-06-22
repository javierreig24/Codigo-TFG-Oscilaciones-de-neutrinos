
#Codigo SK con MCMC
rm(list=ls())
library(latex2exp)
library(extrafont) 
font_import(pattern = "lmroman*") 
#En este algoritmo resumiré la mayoría de aspectos vistos en este algoritmo de cadenas 
loadfonts(device = "win")
par(family = "LM Roman 12")
windowsFonts("LM Roman 12"= windowsFont("LM Roman 12"))


theta = pi/4
E = 10 #GeV
deltam = 2.5*10^(-3) #Electronvolts al cuadrado
R = 6371 #Radio de la Tierra (km)
Ratm <- R+15 #Radio de los neutrinos atmosféricos (km)
Rdet <- R-1 #Profundidad del detector (km)


oscil <- function(L, angulo, masacuad){
  y = 1- (sin(2*angulo))^2*(sin(1.27*(masacuad*L)/(E)))^2
  return(y)
}

L <-function(c){ #Para generar los puntos como argumento el coseno del ángulo
  y = sqrt(Ratm^2-Rdet^2*(1-c^2)) - Rdet*c
  return(y)
}

#Generamos los valores de la función a parte para que no sean distintos cada vez
#que se llama a la función
valoresfun <-vector()
c2 = seq(-1, 1, length.out=10)
for(k in 1:10){
  valoresfun<- c(valoresfun,rnorm(1,oscil(L(c2[k]),theta,deltam),0.05))
}

testf <- function(xx){ #Vamos a coger nuestra función de oscilaciones
  #x1 angulo y x2 deltamasaalcuadrado
  x1 <- xx[1] 
  x2 <- xx[2]

valoresL<-vector()
y = 0;
c2 = seq(-1, 1, length.out=10)
for (k in 1:10){
  c = c2[k] #puntos equiespaciados ente -1 y 1 del coseno
valoresL = c(valoresL, L(c)) #Almacena los valores de L
y <- y +(valoresfun[k] - oscil(L(c),x1,x2))^2 #Ya hemos generado nuestra función
}
return(y)
}

simulated_annealing <- function(func, s0, niter , step) {
  
  # Initialize
  ## s stands for state
  ## f stands for function value
  ## b stands for best
  ## c stands for current
  ## n stands for neighbor
  s_b <- s_c <- s_n <- s0 #Todos iguales al principio que es lo que damos al valor inicial
  f_b <- f_c <- f_n <- func(s_n) 
  #vector to accumulate results
  track1<-vector() 
  track2<-vector()
  track1<-c(track1,s0[1])#c función que combina los argumentos
  track2<-c(track2,s0[2])
  # message("It\tBest\tCurrent\tNeighbour\tTemp")
  #  message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", 0L, f_b, f_c, f_n, 1))
 
  for (k in 1:niter) {     
    Temp <- (1 - step)^k #Temperatura se va haciendo cada vez más sensible
    #Temp <- 1/log(1 + k) #Temperatura se va haciendo cada vez m??s sensible
    # consider a random neighbor
    s_n[1] <- rnorm(1, s_c[1], .004) 
    
    while (s_n[1] < 0.6 || s_n[1] > 0.9){
      s_n[1]<-rnorm(1,s_c[1], .004)
    }
    s_n[2] <- rnorm(1, s_c[2], 0.00004) 
    while (s_n[2]<0.0015 || s_n[2]>0.0035){
      s_n[2]<-rnorm(1,s_c[2], 0.00004)
    }
    f_n <- testf(s_n) #Calcula el nuevo valor
    # update current state
    if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
      #Se impone que el vecino sea menor que el original y además que un número aleatorio constante sea menor que esa exponencial con la temperatura
      s_c <- s_n
      f_c <- f_n
      track1<-c(track1,s_n[1]) #añade el de antes más el nuevo
      track2<-c(track2,s_n[2])
    }
    # update best state
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n 
    }
    # message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", k, f_b, f_c, f_n, Temp))
  }
  return(list(iterations = niter, best_value = f_b, best_state = s_b, values1 = track1, values2 = track2))
}

repetir<-function(ss){
  sol<-simulated_annealing(testf, ss ,30000,0.07);
  y<-rbind(sol$values1,sol$values2)
  ybest = sol$best_state
  return (list(values=y, best=ybest))}



bestvalue = 5;

#De un gran número de cadenas que generamos seleccionamos el mejor valor

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

#Representamos ahora los 100 datos y la función a la que se han ajustado.
c1 = seq(-1, 1, length.out=1000);
c2 = seq(-1, 1, length.out=10);
par(mfrow=c(1,2))
plot(c1,oscil(L(c1),ybest[1],ybest[2]),type='l', xlab=expression(italic('x')),ylab=expression(italic('y')), cex.axis = 1)
points(c2, oscil(L(c2),theta,deltam))

print(valoresfun)

#Codigo DUNE MCMC
rm(list=ls())
th12 = 0.5903
th23 = 0.866
th13=0.15
m21=7.39*10^(-5)
m31=2.451*10^(-3)
dens=2.848
a = (1/3500)*(dens/3)
L = 1285

oscil <- function(E, the12,the23,the13,mas21,mas31){
  y = ((sin(the23))^2*(sin(2*the13))^2*(sin(1.27*mas31*L/E-a*L))^2*(1.27*mas31*L/E)^2)/(1.27*mas31*L/E-a*L)^2 + sin(2*the23)*sin(2*the13)*sin(2*the12)*sin(1.27*mas31*L/E-a*L)*(1.27*mas31*L/E)*sin(a*L)*(1.27*mas21*L/E)*cos(1.27*mas31*L/E)/(a*L*(1.27*mas31*L/E-a*L))+(cos(the23))^2*(sin(2*the12))^2*(sin(a*L))^2*(1.27*mas21*L/E)^2/(a*L)^2
  return(y)
}

E2 <- seq(0.5,8,length.out=20);
valoresfun <-vector()
for (k in 1:20){
  valoresfun = c(valoresfun, rnorm(1, oscil(E2[k],th12,th23,th13,m21,m31),0.005)) #Almacena valores funcion
}

print(valoresfun)

testf <- function(xx){ #Vamos a coger nuestra función de oscilaciones
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
  
  # Initialize
  ## s stands for state
  ## f stands for function value
  ## b stands for best
  ## c stands for current
  ## n stands for neighbor
  s_b <- s_c <- s_n <- s0 #Todos iguales al principio que es lo que damos al valor inicial
  f_b <- f_c <- f_n <- func(s_n) 
  #vector to accumulate results
  track1<-vector() 
  track2<-vector()
  track1<-c(track1,s0[1])#c función que combina los argumentos
  track2<-c(track2,s0[2])
  # message("It\tBest\tCurrent\tNeighbour\tTemp")
  #  message(sprintf("%i\t%.4f\t%.4f\t%.4f\t%.4f", 0L, f_b, f_c, f_n, 1))
  
  for (k in 1:niter) {     
    Temp <- (1 - step)^k #Temperatura se va haciendo cada vez más sensible
    #Temp <- 1/log(1 + k) #Temperatura se va haciendo cada vez más sensible
    # consider a random neighbor
    s_n[1] <- rnorm(1, s_c[1], .004) 
    while (s_n[1] < 0.55 || s_n[1] > 0.61){
      s_n[1]<-rnorm(1,s_c[1], .004)
    }
    
    s_n[2] <- rnorm(1, s_c[2], 0.004) 
    while (s_n[2]<0.83 || s_n[2]>0.89){
      s_n[2]<-rnorm(1,s_c[2], 0.004)
    }
    s_n[3] <- rnorm(1, s_c[3], 0.004) 
    while (s_n[3]<0.11 || s_n[3]>0.17){
      s_n[3]<-rnorm(1,s_c[3], 0.004)
    }
    s_n[4] <- rnorm(1, s_c[4], 0.0000004) 
    while (s_n[4]<6.7*10^(-5) || s_n[4]>7.5*10^(-5)){
      s_n[4]<-rnorm(1,s_c[4], 0.0000004)
    }
    s_n[5] <- rnorm(1, s_c[5], 0.00004) 
    while (s_n[5]<0.0022 || s_n[5]>0.0028){
      s_n[5]<-rnorm(1,s_c[5], 0.00004)
    }
    
    f_n <- testf(s_n) #Calcula el nuevo valor
    # update current state
    if (f_n < f_c || runif(1, 0, 1) < exp(-(f_n - f_c) / Temp)) {
      #Se impone que el vecino sea menor que el original y además que un número aleatorio constante sea menor que esa exponencial con la temperatura
      s_c <- s_n
      f_c <- f_n
      track1<-c(track1,s_n[1]) #añade el de antes más el nuevo
      track2<-c(track2,s_n[2])
    }
    # update best state
    if (f_n < f_b) {
      s_b <- s_n
      f_b <- f_n 
    }
  }
  return(list(iterations = niter, best_value = f_b, best_state = s_b, values1 = track1, values2 = track2))
}

s0 = vector()
s0[1] = runif(1,0.55,0.61)
s0[2] = runif(1,0.83,0.89)
s0[3] = runif(1,0.11,0.17)
s0[4] = runif(1,6.7*10^(-5),7.5*10^(-5))
s0[5] = runif(1,0.0022,0.0028)

sol = simulated_annealing (testf,s0,1,0.07)


repetir<-function(ss){
  sol<-simulated_annealing(testf, ss ,50000,0.07);
  y<-rbind(sol$values1,sol$values2)
  ybest = sol$best_state
  return (list(values=y, best=ybest))}

bestvalue = 5;

#De un gran número de cadenas que generamos seleccionamos el mejor valor

for(k in 1:300){
  y<-repetir(c(runif(1,0.55,0.61),runif(1,0.83,0.89),runif(1,0.11,0.17),runif(1,6.7*10^(-5),7.5*10^(-5)),runif(1,0.0022,0.0028)))$best;
  if(testf(y) < bestvalue ){ #Esto nos devolverá el mejor mínimo de todas
    ybest<-y
    bestvalue<-testf(y)
  }
}

E1 <- seq(0.5, 8, length.out=1000);
E2 <- seq(0.5,8,length.out=20);
plot(E1,oscil(E1,ybest[1],ybest[2],ybest[3],ybest[4],ybest[5]),type='l',  xlab=expression(italic('x')),ylab=expression(italic('y')), cex.axis=1 )
points(E2, oscil(E2,th12,th23,th13,m21,m31))



print(ybest)
print(bestvalue)




