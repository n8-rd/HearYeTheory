#Gomulkiewicz and Holt 95 model(s) of evolutionary rescue

#change in environment at time t0
#absolute fitness W0 < 1 so population starts to crash
#density at first is N0
#demographic change is Nt+1 = Wt*Nt under a simple exponetial model of growth
#critical density below which demographic stochasticity leads to extinction Nc
#how long to reach Nc: tE = (log(Nc) - log(N0))/log(W0)
#mean fitness increaes by delta each 
#how long to reach Wt > 1: tR = (1-W0)/delta
#if tR < tE rescure happens

#presistence is decided by a race between two processes: 
#demographic collapse and evolutionary resuce

#An example
N = 500 #starting population size
Nc = 75 #threshold level below which extinction is likely
W = 0.7 #starting fitness, right after onset of stressful environment
delta = 0.006 #fixed amount that fitness increases each generation
gens = 100 #how many generations to look at
K = 510 #the carrying capacity for a model of logistic population growth

#stressand rescue function
rescue <- function(N, W, delta, K, gens){
    Ns = rep(0, gens)
    for (i in 1:gens){
        W.prime = W+delta
        #simple exponential growth
        #N.prime = N*W.prime
        #or logistic growth
        N.prime = N + (W.prime-1)*N*(1-(N/K))
        Ns[i] <- N.prime
        W = W.prime
        N = N.prime
    }
    return(Ns)
}
par(mfrow=c(1,3))
Ns1 <- rescue(500,0.7,0.004,510,gens)
plot(Ns1~seq(1,gens), type='l', lwd=3, col='red', main='Delta = 0.004', xlab='Generations', ylab='Population size', ylim=c(0,510))
abline(h=Nc, lwd=2, lty=2)
text(5,100,'Nc')

Ns2 <- rescue(500,0.7,0.006,510,gens)
plot(Ns2~seq(1,gens), type='l', lwd=3, col='red', main='Delta = 0.006', xlab='Generations', ylab='Population size', ylim=c(0,510))
abline(h=Nc, lwd=2, lty=2)
text(5,100,'Nc')

Ns3 <- rescue(500,0.7,0.008,510,gens)
plot(Ns3~seq(1,gens), type='l', lwd=3, col='red', main='Delta = 0.008', xlab='Generations', ylab='Population size', ylim=c(0,510))
abline(h=Nc, lwd=2, lty=2)
text(5,100,'Nc')
