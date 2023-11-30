#evolve sex ration in a structured metapopulation
#start by implementing the Wilson and Colwell model of sex ratio evolution
#Then we'll modify it to model the evolution of dispersiveness

##model parameters
Seeder <- 2 #size of founder group
Patches <- 150 #number of host individuals that need colonizing
E <- 7 #number of eggs laid by each female that survive to adulthood
k1 <- 0.33 #proportion of male progeny of fAA or fAa females (A is dominant)
k2 <- 0.5 #proportion of male progeny of faa females
G <- 5 #number of generations between dispersal
Dispersals <- 30 #number of rounds of dispersal

##Model variables
N <- 0 # global population size
K <- 0.5 # starting global proportion of males in global population
FA <- 0.1 # starting global frequency of allele A in the global population


#Function for within-group selection
evolin <- function(Nf,Nm,p1.hat){
  N <- Nf + Nm
  E <- 7
  k1 <- 0.33 #proportion of male progeny of fAA or fAa females (A is dominant)
  k2 <- 0.5 #proportion of male progeny of faa females
  
  #get initial genotype frequencies
  p1 <- p1.hat
  p2 <- 1-p1.hat
  p11 <- p1^2
  p12 <- 2*p1*p2
  p22 <- p2^2
  P = p11 + p12 + p22
  
  #get expected number of females
  Nf.p11.b <- Nf*E*(1-k1) * p1 * (p11 + .5*p12)
  
  Nf.p12.b <- Nf*E*((1-k1) * p11 * p2 
    + (1-k1) * .5*p12
    + (1-k2) * p22 * p1)
  
  Nf.p22.b <- Nf*E*((1-k2) * p2 * p22 
    + (1-k1) * .5*p12)

  Nf.b <- Nf.p11.b + Nf.p12.b + Nf.p22.b
  
  #get expected numbr of males
  Nm.p1.b <- Nf*E*k1*(p11 + .5*p12)
  Nm.p2.b <- Nf*E*(.5*p12*k1 + p22*k2)
  Nm.b <- Nm.p1.b + Nm.p2.b
  
  N.b <- Nf.b + Nm.b
  K.b <- Nm.b / N.b
  
  #get expected genotype frequencies for each sex
  p11.b <- Nf.p11.b / Nf.b
  p12.b <- Nf.p12.b / Nf.b
  p22.b <- Nf.p22.b / Nf.b
  p1.b <- Nm.p1.b / Nm.b
  p2.b <- Nm.p2.b / Nm.b
  
  #get overall feq of allele A1
  p.hat.prime <- (Nf.b*(2*p11.b + p12.b) 
    + Nm.b*p1.b) / (2*Nf.b + Nm.b)
  
  return(c(N.b,Nf.b,Nm.b,p.hat.prime))
}

#Function to initialze host group with dispersal from global pool
seed.group <- function(Seeder,K,FA){
    Fa <- 1-FA
    propFem <- 1-K
    Nf <- ceiling(rbinom(1,Seeder,propFem))
    if (Nf < 1){
        Nf <- 1
    }
    Nm <- Seeder - Nf
    if (Nm < 1){
        Nm <- 1
        Nf <- Nf - 1
    }
    alleles <- (2*Nf) + Nm
    local.fA <- rbinom(1,alleles,FA) / alleles
    local.fa <- 1-local.fA
    return(data.frame(cbind(Nf,Nm,local.fA,local.fa)))
}

#track global fA, K, and N
freqs <- numeric(length=Dispersals)
cakes <- numeric(length=Dispersals)
snakes <- numeric(length=Dispersals)

for (i in 1:Dispersals){ #for each round of dispersal
  global.Ns <- numeric(length=Patches)
  global.fAs <- numeric(length=Patches)
  global.Ks <- numeric(length=Patches)
  print(FA)
  for (j in 1:Patches){ #for each patch
    ##initialize each patch with a random sample from the global pool
    founders <- seed.group(Seeder,K,FA)
    Nf <- founders$Nf[1]
    Nm <- founders$Nm[1]
    #track within-group variables
    patch.N <- Seeder
    patch.fA <- founders$local.fA[1]
    #patch.K <- Nm/(Nf+Nm)
    #do selection within each path for G generations
    for (x in 1:G){
        states <- evolin(Nf=Nf,Nm=Nm,p1.hat=patch.fA)
        N.prime <- states[1]
        fA.prime <- states[4]
        #K.prime <- states[2]
        patch.N <- N.prime
        patch.fA <- fA.prime
        #patch.K <- K.prime
        Nf <- states[2]
        Nm <- states[3]
        }
    global.Ns[j] <- patch.N
    global.fAs[j] <- patch.fA
    #global.Ks[j] <- patch.K
  }
  N <- sum(global.Ns)
  FA <- weighted.mean(global.fAs,w=global.Ns)
  #K <- weighted.mean(global.Ks,w=global.Ns)
  freqs[i] <- FA
  #cakes[i] <- K
  snakes[i] <- N
}

gen <- seq(1,Dispersals,1)
plot(gen, freqs, type='l', col='magenta', lwd=2, xlab='generation', ylab='frequency of allele A', main='Evolution of biased Sex ratio')


#fAAf <- (fA^2)
#fAaf <- 2*fA*fa
#faaf <- 1 - (fAAf + fAaf)
#fAm <- fA
#fam <- fa
#NmA <- Nm*fAm
#Nma <- Nm - NmA
#NfAA <- Nf*fAAf
#NfAa <- Nf*fAaf
#Nfaa <- Nf*faaf


#fA = 1/N #frequency of allele that biases sex ratio (or dispersal)
#fa = 1 - fA # frequency of null allele
#fAA = fA^2
#fAa = 2*fA*fa
#faa = fa^2
