boltzmannPlot <- function(DF, RGB=c(0.2, 0.0, 0.8)){
	fit <- nls( formula=N ~ N0 * exp(-Ek/T), data=data.frame(N=DF$PDF, Ek=DF$E_k), weight=1.0/DF$ePDF^2, start=list(N0=0.6, T=400))
	N0 <- coef(summary(fit))[1,1]	# Best N0
	T0 <- coef(summary(fit))[2,1]	# Best T0
	eN0 <- summary(fit)$coefficients[3]	# Error in N0
	fit <- nls( formula=N ~ N0 * exp(-Ek/T), data=data.frame(N=DF$PDF, Ek=DF$E_k), weight=1.0/DF$ePDF^2, start=list(T=T0))
	T0 <- summary(fit)$coefficients[1]
	eT0 <- summary(fit)$coefficients[2]
	minT <- T0 - eT0
	maxT <- T0 + eT0
	Eplot <- seq(2, 90)
	UpperPDF <- (N0 + eN0) * exp(- Eplot / maxT)
	LowerPDF <- (N0 - eN0) * exp(- rev(Eplot) / minT)
	modelPDF <- N0 * exp(- Eplot / T0)
	polygon( c(Eplot, rev(Eplot)), c(UpperPDF, LowerPDF), col=rgb(RGB[1], RGB[2], RGB[3], 0.2), border=NA)
	lines( Eplot, modelPDF, col=rgb(RGB[1], RGB[2], RGB[3], 0.75), lwd=2) ; text(80, modelPDF[length(modelPDF)], sprintf('%.0f ± %.0f K', T0, eT0), cex=0.75, col=rgb(RGB[1], RGB[2], RGB[3]))
	points( DF$E_k, DF$PDF, pch=20, col=rgb(RGB[1], RGB[2], RGB[3]))
	arrows( DF$E_k, DF$PDF - DF$ePDF,  DF$E_k,  DF$PDF + DF$ePDF, length=0, col=rgb(RGB[1], RGB[2], RGB[3]))
	cat(sprintf('T = %.0f ± %.0f K / N = %.2f ± %.2f\n', T0, eT0, N0, eN0))
	return
}


#-------- SO line parameters
# J_low : lower J number
# N_low : lower N number
# E_k   : lower energy state (K)
# EW    : Equivalent width [km/s]
# ER    : error of EW
# restfreq : rest frequency (GHz)
SO2211  <- list(J_low=1, J_up=2, N_low=1, N_up=2, E_k=15.1819, EW=1.282, ER=0.067, restfreq=86.09395, name='SO J(N) = 2(2) - 1(1) v=0')
SO3221  <- list(J_low=2, J_up=3, N_low=1, N_up=2, E_k= 4.4600, EW=2.328, ER=0.109, restfreq=99.29987, name='SO J(N) = 3(2) - 2(1) v=0')
SO3322  <- list(J_low=2, J_up=3, N_low=2, N_up=3, E_k=19.3137, EW=2.396, ER=0.113, restfreq=129.13892, name='SO J(N) = 3(3) - 2(2) v=0') # subtract SO2 12(2,10)-11(2,10)
SO5544  <- list(J_low=4, J_up=5, N_low=4, N_up=5, E_k=33.7748, EW=4.833, ER=0.075, restfreq=215.22065, name='SO J(N) = 5(5) - 4(4) v=0')
SO8776  <- list(J_low=7, J_up=8, N_low=6, N_up=7, E_k=64.8931, EW=4.677, ER=0.503, restfreq=340.71416, name='SO J(N) = 8(7) - 7(6) v=0')
SO8877  <- list(J_low=7, J_up=8, N_low=7, N_up=8, E_k=70.9573, EW=8.128* 7.9/(7.9+1.2), ER=0.211, restfreq=344.31061, name='SO J(N) = 8(8) - 7(7) v=0') # correction to exclude H15CN contamination
#SO8978  <- list(J_low=7, J_up=8, N_low=8, N_up=9, E_k=62.1444, EW=8.781, ER=0.290, restfreq=346.52848, name='SO J(N) = 8(9) - 7(8) v=0')
SO8978  <- list(J_low=7, J_up=8, N_low=8, N_up=9, E_k=62.1444, EW=8.556, ER=0.106, restfreq=346.52848, name='SO J(N) = 8(9) - 7(8) v=0') # using EquiWidth(DF, c(1250, 1310, 1700, 1740) )


multiTransition <- function(lineList){
	lineNum <- length(lineList)
	JL <- JU <- NL <- NU <- EK <- EW <- ER <- numeric(lineNum)
	for(line_index in 1:lineNum){
		JL[line_index] <- lineList[[line_index]]$J_low
		JU[line_index] <- lineList[[line_index]]$J_up
		NL[line_index] <- lineList[[line_index]]$N_low
		NU[line_index] <- lineList[[line_index]]$N_up
		EK[line_index] <- lineList[[line_index]]$E_k
		EW[line_index] <- lineList[[line_index]]$EW
		ER[line_index] <- lineList[[line_index]]$ER
	}
	lineDF <- data.frame(J_up=JU, J_low=JL, N_up=NU, N_low=NL, E_k=EK, EW=EW, ER=ER)
	lineDF$PDF <- lineDF$EW / (2* lineDF$N_low + 1)
	lineDF$ePDF <- lineDF$ER / (2* lineDF$N_low + 1)
	return( lineDF )
}

lineDF <- multiTransition( list(SO2211,SO3221,SO3322,SO5544,SO8776,SO8877,SO8978) )

pdf('SOboltzmann.pdf', width=10, height=7)
plot( lineDF$E_k, lineDF$PDF, pch=20, ylim=c(0, 0.9), xlim=c(0.0,80), xlab=expression('E'[k]), ylab='EW / (2 N + 1)', main='SO Boltzmann Diagram', type='n')

lineDF <- multiTransition( list(SO5544,SO8776,SO8877,SO8978) )
boltzmannPlot(lineDF, c(0.2, 0.0, 0.8))
lineDF <- multiTransition( list(SO2211,SO3221,SO3322) )
boltzmannPlot(lineDF, c(0.8, 0.0, 0.2))


dev.off()