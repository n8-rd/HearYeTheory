#Here's the Leven 1953 model of course-grained selection over space
#We've got a di-allelic locus that affects viability differently in different
#resource patches. p is global proportion of allele S. q is the global
#proportion of allele F. In each patch, dispersers survive to adulthood
#according to patch-specific genotype viabilities

#we'll put our patch-specific genotype viabilities in a matrix
viabilities <- c(0.696, 1, 1.424, 0.898, 1, 1.012, 1.519, 1, 0.88, 0.913, 1, 0.976)
vm <- data.frame(matrix(viabilities, nrow=4, byrow=T))
colnames(vm) <- c('SS', 'SF', 'FF')
rownames(vm) <- c('p1', 'p2', 'p3', 'p4')
print(vm)

#Each patch contributes a fixed proportion, ci, of offspring to dispersal pool
cis <- c(0.4,0.05,0.5,0.05)

#here's a function that caclulate the global change in p, which is 
#basically the weighted average of local, within patch changes, weighted by
#the ci terms
get.global.delta <- function(p, vm, cis){
    q <- 1-p
    xs <- numeric(length=4) #a vector for the local changes in p across 4 patches
    for (i in 1:4){ #for each resource patch
        Vi <- vm[i,1] # get the local viability of SS genotypes
        Wi <- vm[i,3] # get the local viability of FF genotypes
        ## In this model heterozygote viabilities are set to 1
        ci <- cis[i] #get the proportion that the local patch contributed to the global pool of dispersers
        Wi.bar <- (p^2)*Vi + 2*p*q + (q^2)*Wi # this is average local viability
        p.prime <- ((p^2)*Vi + p*q) / Wi.bar # p in the next generation. This is just like in our simple diploid model of viability selection
        delta <- p.prime - p #the change in p from when oddspring settle to when adults reproduce
        x <- ci * delta #weight the change in p by the proportion of offspring the ocal patch contributes to the dispersal pool
        xs[i] <- x
    }
    #del <- p*q*sum(xs) 
    del <- sum(xs) #The global change in p
    return(del)
}


ps <- seq(0,1,0.01)
gds <- numeric(length=length(ps))

for (i in 1:length(ps)){
    p1 <- ps[i]
    gd <- get.global.delta(p1, vm, cis)
    gds[i] <- gd
}

plot(gds~ ps, type='l', lwd=3, col='purple', ylab='Change in S', xlab='Frequency of allele S', main='Soft selection in a mixed environment')
abline(h=0, lty=2, lwd=2, col='gray')

