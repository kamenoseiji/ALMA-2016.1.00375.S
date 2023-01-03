CoreSPIX <- -0.5
JetSPIX <- -3
FFAspec <- function(Sc, alpha, tau0, nu){
	# nu : frequency [GHz]
	# alpha : spectral index d log(S) / d log(nu)
	# Sc : Unabsorbed flux density at 1 GHz
	# tau0 : FFA optical depth at 1 GHz (= 0.46 ne^2 Te^-3/2 dL)
	return(Sc * nu^alpha * exp(-tau0* nu^{-2.1}))
}

K53 <- function(x){ besselK(x, 5/3) }
F <- function(x){
	y <- numeric(length(x))
	for(index in 1:length(x)){
		y[index] <- x[index]* integrate(K53, lower = x[index], upper = 1000)$value
	}
	return( y )
}

SSAspec <- function(Sc, alpha, peakFreq, nu){
	return( Sc* (nu/peakFreq)^2.5 * (1.0 - exp( -(nu/peakFreq)^(alpha - 2.5))) )
}

SynchrotronBreakSpec <- function(Sc, alpha, breakFreq, nu){
	return( Sc* ( nu^(alpha) * exp(- nu/breakFreq) ))
	# return( Sc* ( nu^(alpha) * exp(- nu/breakFreq) + nu^(alpha - 2.5) * (1.0 - exp(- nu/breakFreq))))
}

FreqPoints <- c( 88.8,  86.9,  86.9,  85.0,  97.0,  98.7,  98.6, 100.4, 130.1, 128.2, 128.3, 126.5, 138.3, 140.1, 139.9, 141.7, 354.1, 356.0, 357.9, 351.0, 343.2, 344.2, 341.8, 338.0, 231.7, 213.1)
FluxPoints <- c(1.190, 1.205, 1.205, 1.221, 1.125, 1.112, 1.114, 1.100, 0.875, 0.883, 0.883, 0.891, 0.843, 0.836, 0.837, 0.831, 0.443, 0.445, 0.446, 0.441, 0.443, 0.443, 0.443, 0.444, 0.491, 0.492)
FluxError <-  c( 0.004256,  0.004256,  0.004256,  0.004256,  0.002932,  0.002932,  0.002932, 0.002932, 0.002390, 0.002390, 0.000157, 0.000157, 0.001359, 0.001359, 0.001359, 0.001359, 0.001137, 0.001137, 0.001137, 0.001137, 0.000901, 0.000901, 0.000901, 0.000901, 0.001359, 0.001359)
ContDF <- data.frame(freq=FreqPoints, flux=FluxPoints, err=FluxError)

fit <- nls(formula = flux ~ a* freq^(CoreSPIX)* exp(-b* freq^(-2.1)) + c* freq^(JetSPIX)* exp(-d* freq^(-2.1)), data=ContDF, weight=1/err, start=list(a=20, b=100000, c=2e+05, d=10000))

freq <- 10^seq(1, 3, by=0.01)
#CoreParam <- c(0.725, -0.1, 2.5e+04)
#JetParam <- c(1.9e+05, -2.5, 10000) 

CoreParam <- c(coef(fit)['a'], CoreSPIX, coef(fit)['b'])
JetParam <- c(coef(fit)['c'], JetSPIX, coef(fit)['d']) 

Score <- FFAspec(CoreParam[1], CoreParam[2], CoreParam[3], freq)
Sjet  <- FFAspec(JetParam[1], JetParam[2], JetParam[3], freq)

Colors <- c('red', 'blue')
Labels <- c('Jet', 'Core')
pdf('ContSpec.pdf', width=7, height=7)
plot(freq, Score + Sjet, type='l', log='xy', xlim=c(80,400), ylim=c(0.25, 1.25), xlab='Frequency [GHz]', ylab='Flux Density [Jy]', main='NGC 1052 Continuum Spectrum')
arrows(ContDF$freq, ContDF$flux + 5*ContDF$err, ContDF$freq, ContDF$flux - 5*ContDF$err, length=0)
points(ContDF, pch=20)
lines(freq, Sjet, type='l', col=Colors[1])
lines(freq, Score, type='l',  col=Colors[2]) 
legend("topright", legend = Labels, col = Colors, lty=1)
dev.off()