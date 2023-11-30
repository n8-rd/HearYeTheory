d <- read.csv('burgerLynch.csv',header=F)
par(mfrow=c(2,1))
plot(d$V2 ~ d$V1, type='l', lwd=3, xlab='generation', ylab='Lag', main='Lag of mean phenotype behind O', col='purple')
plot(d$V3 ~ d$V1, type='l', lwd=3, xlab='generation', ylab='N', main='Evolution of Population Size', col='green')
