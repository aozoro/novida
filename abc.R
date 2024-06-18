C_A <- function(C,A,tipo,...){
  if(tipo=="exp"){
    aux<-C*exp(-A/C)
  }
  if(tipo=="LN"){
    aux<-C*(1-pnorm((log(A)-(mu+sigma^2))/sigma))-A*(1-pnorm((log(A)-mu)/sigma))
  }
  if(tipo=="Gamma"){
    aux<-C*(1-pgamma(A,a+1,b))-A*(1-pgamma(A,a,b))
  }
  return(aux)
}

a = 1.058925; b = 0.003134
mu = 5.357209; sigma = 0.8325546
C = 300
A = 500
R = 300
M = 900
alfa = 0.8

CX <- C_A(C,R,"LN",mu,sigma)+R*(1-pnorm((log(R)-mu)/sigma));CX
1-CX/C

f<-coverage(dexp,pexp,deductible=500,coinsurance=0.9,limit=2000,per.loss=TRUE)
curve(dlnorm(x,mu,sigma),xlim=c(0,3000),ylim=c(0,0.001))
curve(f(x,mu,sigma),xlim=c(0.01,3000),col="red",add=TRUE)

Pb<-100*medn/sum(Pi.av*b);Pb
Pc<-b/100*Pb;Pc
RSAL<-(mean(b)-b[1])/(b[nf]-b[1]);RSAL

CV<-(sum(Pi.av*(b-mean(b))ˆ2)ˆ.5)/mean(b);CV



###3

#Supongamos que la variable aleatoria N, número de siniestros en un periodo, puede 
#tomar valores 0, 1 y 2, con probabilidades 0,7; 0,2 y 0,1, respectivamente.  Supongamos 
#que la variable aleatoria X, cuantía de un siniestro, puede tomar valores 1 y 2, con 
#probabilidades 0,8 y 0,2, respectivamente. 
#Se pide, para la variable aleatoria S, coste total por póliza y periodo: 

##Calcular el valor de la función de distribución de S para los valores 0,1,2,3,4. 
  
#cargar el package actuar 
library(actuar) 
#crear vector con probabilidades N 
probn<-c(0.7,0.2,0.1) 
#crear vector con probabilidades X 
probx<-c(0,0.8,0.2) 
#cálculo de función distribución S 
Fscon<-aggregateDist("convolution",model.freq=probn,model.sev=probx) 
#cálculo valores de función distribución S 
Fscon(0:4) 

###b) Calcular la probabilidad que S tome valores 0,1,2,3,4.


 
#cálculo probabilidades S 
c(Fscon(0),diff(Fscon(0:4))) 



#2) Supongamos que la variable aleatoria X, cuantía de un siniestro, puede tomar valores 1 
#y 2, con probabilidades 0,3 y 0,7, respectivamente. 
#Supongamos que la variable aleatoria N, número de siniestros en un periodo: 
#• sigue una distribución de Poisson, con parámetro  0,8 λ= 
#• sigue una distribución Binomial Negativa con parámetros 1,0588655 a= y 
#0,1362035 b= 
#Se pide, para la variable aleatoria S, coste total por póliza y periodo: 
 
####a) Aplicando el método de recurrencia de Panjer, calcular la probabilidad que S 
#tome valores 0,1,2,3,4. 

#crear vector con probabilidades X 
probx<-c(0,0.3,0.7) 
#cálculo de función distribución S con Poisson compuesta 
FsPois<-aggregateDist("recursive",model.freq="poisson",lambda=0.8,model.sev=probx) 
#cálculo valores de función distribución S 
FsPois(0:4) 
#cálculo probabilidades S con Poisson compuesta 
R> FsPois(0) 
dpois(0,lambda=0.8) 
diff(FsPois(0:4)) 

#cálculo de función distribución S con Binomial negativa compuesta 
FsNbinom<-aggregateDist("recursive",model.freq="negative binomial", size=1.0588655, prob=1/(1+0.1362035),model.sev=probx) 
FsNbinom(0:4) 

#cálculo probabilidades S con Binomial negativa compuesta 
FsNbinom(0) 
dnbinom(0,size=1.0588655,prob=1/(1+0.1362035)) 
diff(FsNbinom(0:4)) 








#####b) Aplicando la transformada rápida de Fourier, calcular la probabilidad que S 
#tome valores 0,1,2,3,4. 

n<-32;p<-rep(0,n);p[2:3]<-c(0.3,0.7);lambda<-0.8;a<-1.0588655;b<-0.1362035 
fsFourierPoisson<-Re(fft(exp(lambda*(fft(p)-1)),inverse=TRUE))/n 
fsFourierPoisson[1:5] 

fsFourierNbinom<-Re(fft((1/(1-b*(fft(p)-1))^a),inverse=TRUE))/n 
fsFourierNbinom[1:5] 


#Supongamos que el número de siniestros sigue una distribución de Poisson de media 
#0,95 y que la cuantía de un siniestro tiene los siguientes momentos ordinarios: 
 
#• 1
# 5,375 α= 
#• 2
# 45,89 α= 
#• 3
# 488,05 α= 
#• 4
# 5899,64285 α= 
 
#Sabiendo que la variable cuantía de un siniestro está en unidades de 100 euros, 
#calcular: 
 
####a) El valor estimado para la esperanza, varianza y  coeficiente de asimetría de la 
#variable aleatoria S, coste total por póliza y periodo. 

 #crear variables con parámetro lambda y momentos de X 
lambda<-0.95;medx<-5.375;medx2<-45.89;medx3<-488.05 
#cálculo media, varianza y ceoficiente asimetría de S 
meds<-lambda*medx;meds 
vars<-lambda*medx2;vars 
assims<-(lambda*medx3)/(lambda*medx2)^(3/2);assims 



####b) A través de las aproximaciones Normal-Power y Wilson-Hilferty, calcular las 
#siguientes probabilidades: 
#• Probabilidad  que el coste total supere los 800 euros. 
#• Probabilidad que el coste total sea inferior a 100 euros. 
#• Probabilidad que el coste  total sea inferior a [ ] 1,05ES ⋅ . 
#• Probabilidad que el coste total supere los 2000  euros.


#cálculo de función distribución S con Normal-Power 
FsNP<-aggregateDist("npower",moments=c(meds,vars,assims)) 
#cálculo valores de función distribución S  
1-FsNP(8) 
FsNP(1) 
FsNP(1.05*meds) 
1-FsNP(20) 


