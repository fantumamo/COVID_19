
########################################################



### R code for  COVID 19 modeling in Ethiopia ###

##########################################################
library('reshape2')
library('ggplot2')
library(deSolve)
library(dplyr)
covid_model=function(current_timepoint, state_values, parameters)
{	
  
  ### state Variables ###
  S=state_values[1] # susceptible
  V=state_values[2] # Vaccinated 
  E=state_values[3] # Exposed
  Ia=state_values[4] # Asymptomatic Infectious
  Is=state_values[5] # Symptomatic Infectious
  Ih =state_values[6] # Isolated Infectious
  Iu =state_values[7] # Intensive care unit
  R=state_values[8]   # Recovered
  D=state_values[9]    # Died
  with(
    as.list(parameters),
    {
      
      ### computing derivative ###
      
      dS= -((beta*(Is+psi*Ia))/N)*S
      dV=tau*S
      dE= ((beta*(Is+psi*Ia))/N)*S -((sigma)*E)
      dIa= (sigma)*E -(gamma*Ia + eta*Ia)
      dIs=((sigma)*E +eta*Ia)  - (delta*Is)
      dIh= (delta*Is) - (epsilon*Ih+kappa*Ih)
      dIu= (epsilon*Ih)-(alpha*Iu + nu*Iu)
      dR= (gamma*Ia)  + (nu*Iu)
      dD= (kappa*Ih) + (nu*Iu)
      
      ### Combine results ###
      results=c(dS, dV, dE,dIa,dIs,dIh,dIu,dR,dD)
      list(results)
    }
    
  )
}
### setting parameter values###

beta= 0.725 #Transmission rate
tau= 0.0129 #The rate at which Susceptible population get vaccinated
gamma= 0.20 #Recovery rate among asymptomatic infectious individuals
sigma=0.1923 ##Rate of progression from exposed to symptomatic and asymptomatic infectious
eta=0.40 # rate of release from the asymptomatic to symptomatic infectious 
delta=0.8 ##Rate of isolation among symptomatic infectious individuals
kappa=0.015 #Death rate among isolated covid 19 patients
epsilon =0.01 #Rate at which Isolated patients' needs admission to ICU	
alpha=0.80  #Rate of Recovery after getting critical
nu= 0.092 #Death rate after getting critical
psi =0.5 # Transmission coefficient for asymptomatic individuals

parameter_list=c(beta=beta, tau=tau,  psi=psi, gamma=gamma, sigma=sigma,eta=eta, delta=delta, kappa=kappa, epsilon=epsilon, alpha=alpha,nu=nu)

### Initial values for state variables###

S=0.91790431309
v=0.0124966959
E=0.0230867011
Ia=0.0207780310
Is=0.0023086701
Ih=0.0003833999
Iu=0.00000393011
R=0.023003500
D=0.0000347588
N=S+v+E+Ia+Is+Ih+Iu+R+D
print(N)
initial_values=c(S=S, v=v, E=E, Ia=Ia, Is=Is, Ih=Ih, Iu=Iu, R=R, D=D)
# N=S+v+E+Ia+Is+Ih+Iu+R+D

### Output Time points ###
times=seq(0, 365, by = 7)

### Simulate the SVEIaIsIhIuRD transmission and solve the system of differential equations numerically with lsoda in the  deSolve package ###

output=as.data.frame(lsoda(initial_values, times, covid_model, parameter_list))

print(output)
tail(output)


### Plot the above Result####

melted = melt(output, id.vars="time")
ggplot(data=melted, aes(x=time, y=value, color=variable)) +
  geom_line(size=0.9)


### fill results with vectors###
SV=output$S
VV=output$v
EV=output$E
IaV=output$Ia
IsV=output$Is
IhV=output$Ih
IuV=output$Iu
RV=output$R
DV=output$D
Vtime=output$time
NV=SV+VV+EV+IaV+IsV+IhV+IuV+RV+DV

output_N <- output
output_N$N_t<- NV

### let us compute the incidence ###
output_incidence<-output %>% mutate(incidence = ((IsV+IaV- lag(IsV+IaV))/NV)*100000)
print(output_incidence$incidence)
tail(output_incidence)
plot(output_incidence$time,output_incidence$incidence,type = 'l',)

#### model fitting and parameter estimation using least square estimation technique###
# now find the values of beta and kappa that give the smallest RSS, which represents the best fit to the data.
# save the data in the git folder for us to access
ethiopian.data.I <- read.csv("C:/Users/User/Desktop/ethiopian data I.csv")
print(ethiopian.data.I)
Infected<-ethiopian.data.I$Infected
times_of_Infected<-ethiopian.data.I$date

outputt <- function(time, state, parameters) {
  par <- as.list(c(time, state, parameters))
  with(par, {
    dS= -((beta*(Is+psi*Ia))/N)*S
    dV=tau*S
    dE= ((beta*(Is+psi*Ia))/N)*S -((sigma)*E)
    dIa= (sigma)*E -(gamma*Ia + eta*Ia)
    dIs=((sigma)*E +eta*Ia)  - (delta*Is)
    dIh= (delta*Is) - (epsilon*Ih+kappa*Ih)
    dIu= (epsilon*Ih)-(alpha*Iu + nu*Iu)
    dR= (gamma*Ia)  + (nu*Iu)
    dD= (kappa*Ih) + (nu*Iu)
    results=c(dS, dV, dE,dIa,dIs,dIh,dIu,dR,dD)
    list(results)
  })
}
library(deSolve)

S=0.91790431309
v=0.0124966959
E=0.0230867011
Ia=0.0207780310
Is=0.0023086701
Ih=0.0003833999
Iu=0.00000393011
R=0.023003500
D=0.0000347588

### Initial state values for the above differential equations###
initial_values=c(S=S, v=v, E=E, Ia=Ia, Is=Is, Ih=Ih, Iu=Iu, R=R, D=D)
parameters<-c(0.725,  0.015); names(parameters) <- c("beta"   ,   "kappa")
days <- seq_along(ethiopian.data.I$Infected)

solution <- ode(initial_values, times = days, func = outputt, parms = parameters)
N=1
S=0.91790431309
v=0.0124966959
E=0.0230867011
Ia=0.0207780310
Is=0.0023086701
Ih=0.0003833999
Iu=0.00000393011
R=0.023003500
D=0.0000347588


initial_values=c(S=S, v=v, E=E, Ia=Ia, Is=Is, Ih=Ih, Iu=Iu, R=R, D=D)
## RSS as an R function
RSS <- function(parameters){
  names(parameters) <- c("beta",  "kappa") # parameters must be named
  solution <- ode(initial_values, times = days, func = outputt, parms = parameters)
  I <- solution[, 5] # fifth column of ODE solution
  
  return(sum(Infected - I)^2)
}
optimal_sol <- optim(c(0.5, 0.5), RSS, method = "L-BFGS-B",lower = c(0, 0), upper = c(1, 1))

fitted_pars <- setNames(optimal_sol$par, c("beta", "kappa"))
print(round(fitted_pars,3))









### Effective reproduction number(Re) calculation using Next generation matrix approach after we find F and V matrix.
### Finally  the spectral radius of the matrix F*V inverse gives the following result### 

Re=beta*psi/gamma
print(Re)

