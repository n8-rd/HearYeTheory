#evolve dispersiness in a structured metapopulation
#Uses the Wilson and Colwell (1981) model of the evolution of sex ratios as
# a springboard. But relaxes many of their assumptions.

R = 2 #per capitate growth rate
kA = 0.5 #proportion of A-allele bearing progeny that migrate
ka = 0.3 #proportion of a-allele bearing progeny that migrate
D = 0.1 #host death rate per generation
M = 0.1 #host imigration rate per generation
NM = 0 #size of migrant pool
P = 0.2 #frequency of A allele in migrant pool
H = 150 #number of host patches
#Hex = 1000 #host life expecancy
C = 1e+5 #host carrying capacity
Seed = 10 #number of individuals per host at the start
G = 200 #numer of generations

#function for haploid within-group selection
within <- function(N,p){
    #do logistic population growth
    r = 1 + R*(1-N/C)
    N.prime = N*r #number of progeny
    Nm.kA = N*r*p*kA #number of kA migrants
    Nm.ka = N*r*(1-p)*ka #number of ka migrants
    Nl.kA = N*r*p*(1-kA) #number of kA residents
    Nl.ka = N*r*(1-p)*(1-ka) #number of Ka residents
    Nm = Nm.kA + Nm.ka
    Nl = Nl.kA + Nl.ka
    if (Nl != 0){
        phat.l = Nl.kA / Nl
    } else {
        phat.l = 0
    }
    if (Nm !=0){
        phat.m = Nm.kA / Nm
    } else {
        phat.m = 0
    }
    return(data.frame(cbind(Nm,Nl,phat.l,phat.m)))
    }

#go for G generations
Hn <- rep(10,H) # within-host population sizes
Hp <- rep(0.2,H) # within-host A allele frequencies
Mn <- rep(0,H) # emmigrants from each host
Mp <- rep(0,H) # frequency of A allele in each host's emmigrants
Ha <- rpois(H,2) #age of hosts

Ps <- numeric(length=G)

for (i in 1:G){ #for each generation
    #print(c(i,P))
    for (j in 1:H){ #for each host patch
        #stochastic host death
        if (rbinom(1,1,D)){
            Hn[j] <- 0
            Hp[j] <- 0
        } 
#        #deterministic host death
#        if (Ha[j] >= Hex){
#            Hn[j] <- 0
#            Hp[j] <- 0
#        } 
        else { #the host survives
            #immigration
            if (rbinom(1,1,M)){ #host gets a migrant
                Nhd <- Hp[j]*Hn[j] #previous number of high dispersers
                Hn[j] <- Hn[j] + 1
                if (rbinom(1,1,(1-P))){ #migrant is high disperser
                    Nhd <- Nhd + 1
                    newp <- Nhd / Hn[j]
                    Hp[j] <- newp
                } else { #migrant is low disperser
                    newp <- Nhd / Hn[j]
                    Hp[j] <- newp
                }
            }
            #do within host offspring production
            x <- within(Hn[j],Hp[j])
            Hn[j] <- x$Nl #offspring replace their parents
            Mn[j] <- x$Nm
            Hp[j] <- x$phat.l
            Mp[j] <- x$phat.m
        }
    }
    #hosts get older
    Ha <- Ha + 1
    #New Migrant pool is sum of all host group emmigrants
    NM <- sum(Mn)
    #New frequency of allele A in Migrant pool is weighted mean
    #of host group frequencies
    P <- weighted.mean(Mp,w=Mn)
    Ps[i] <- P
}

gen <- seq(1,G,1)
plot(gen, Ps, type='l', col='magenta', lwd=2, xlab='generation', ylab='frequency of allele A', main='Evolution of dispersal')



    
