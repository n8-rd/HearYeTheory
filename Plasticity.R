##A model of selection on a phenotypically plastic trait under
##spatially heterogeneous selection (an environmental cline)

#growth rate reaction norm parameters
a <- 1.0 #intercept 
c <- 0.33 #slope

#temperature gradient
temps <- seq(0.1,20,0.5) # range of temperatures

#fitness function parameters
B = 3.33
alpha = 0.02

#fitness function
get.fitness <- function(t){
    #argument is temperture
    GAA <- a + c*t
    GAa <- c*t
    Gaa <- -a + c*t
    wAA <- max(0, (1 - alpha*(B-GAA)^2))
    wAa <- max(0, (1 - alpha*(B-GAa)^2))
    waa <- max(0, (1 - alpha*(B-Gaa)^2))
    return(c(wAA,wAa,waa))
}

#selection function
selection <- function(p, t){
    #arguments are frequency of allele a, and temperature
    q <- 1 - p
    ws <- get.fitness(t)
    wAA <- ws[1]
    wAa <- ws[2]
    waa <- ws[3]
    p.prime <- ((wAA*p^2) + (wAa*p*q)) / ((wAA*p^2) + (2*wAa*p*q) + (waa*q^2))
    return(p.prime)
}

ps <- numeric()
gens <- seq(1,100,1)

for (t in temps) {
    p <- 0.5 #start p off at 50%
    for (g in gens){
        p <- selection(p, t)
    }
    ps <- c(ps, p)
}

#plot how p after 100 generations of selection varies across temperatures
plot(ps ~ temps, type='l', col='magenta', lwd=3, xlab='Temperature', ylab='Equilibrium freq. of allele A', main='Adaptiveness of A varies with temperature', ylim=c(0,1))
