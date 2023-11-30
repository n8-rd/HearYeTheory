##Some implementations of classical population growth models in R

##First a model of exponential growth
PureBirth <- function(R, current){
    nex <- R*current
    return(nex)
    }
    
time.steps = seq(1,20)
initial.size = 2 # we start out with just two parents
R = 1.5 #this is our first growth rate parameter, the number of progeny produced by each parent
pop.size <- rep(0,20)
pop.size[1] <- initial.size

for (i in 2:length(time.steps)){
    curr <- pop.size[i-1]
    nex <- PureBirth(R, curr)
    pop.size[i] <- nex
    }

par(mfrow=c(1,2))
plot(pop.size ~ time.steps, type='l', lwd=2, col='green', main='Exponential Growth, R=1.5')

##Next a model of logistic growth
LogisticGrowth <- function(r, K, current){
    nex <- current + r*current*(1-(current/K))
    return(nex)
    }
    
r = 0.5 #this is a different growth rate parameter. r = (R-1). It's basically the per capita change in the number of individuals from one generation to the next.
K = 30 #this is the carrying capacity for the population

for (i in 2:length(time.steps)){
    curr <- pop.size[i-1]
    nex <- LogisticGrowth(r, K, curr)
    pop.size[i] <- nex
    }
    
plot(pop.size ~ time.steps, type='l', lwd=2, col='green', main='Logistic Growth, r=0.5, K=30')

