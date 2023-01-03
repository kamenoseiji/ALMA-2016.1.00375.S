setwd('./')
source('LineList.R')
plotLine <- function(lineSpec, velPlot, Vofs, Vreg, bottom, height){
	# velPlot <- seq(-100, 500, by=100)
	# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	cvel <- 299792.458
	top <- bottom + height
	freq <- (velPlot - Vofs - cvel) * lineSpec$restfreq / (Vreg + cvel)
	# cat(freq)
	segments( min(freq), bottom, max(freq), bottom )
	segments( freq, rep(bottom, length(freq)), freq, rep(top, length(freq)))
	for( tick_index in 1:length(velPlot)){
		if(velPlot[tick_index] == 0){
			text_sd <- 'Vsys'
		} else {
			text_sd <- sprintf(' %4.0f', velPlot[tick_index])
		}
		text( freq[tick_index], top, text_sd, adj=0, cex=0.3, srt=90)
	}
	text( mean(freq), bottom, lineSpec$name, pos=1, offset=0.2, cex=0.3)
}

plotTau <- function(DF, lineSpec, velRange, Vofs, Vreg){
	# velPlot <- seq(-100, 500, by=100)
	# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	cvel <- 299792.458
	DF$VLSR = Vofs + cvel + DF$freq/lineSpec$restfreq * (Vreg + cvel)
	plotDF <- DF[((DF$VLSR >= velRange[1]) & (DF$VLSR <= velRange[2])),]
	plotDF <- plotDF[order(plotDF$VLSR),]
	chGap <- c(0.0, diff(plotDF$VLSR))
	plot(plotDF$VLSR, plotDF$Tau, pch=20, cex=0.5, xlab='LSR Velocity [km/s]', ylab='Optical Depth', main=lineSpec$name, xlim=velRange)
	lines(plotDF$VLSR - 0.5*chGap, plotDF$Tau, type='s', lwd=0.5)
	abline(h=0, lty=2, lwd=0.5, col='grey')
	peakIndex <- which.max(plotDF$Tau)
	text_sd <- sprintf('%s TauMax = %f at %f km/s\n', lineSpec$name, plotDF$Tau[peakIndex], plotDF$VLSR[peakIndex])
	cat(text_sd)
	return(plotDF)
}

freqRead <- function(prefixList){
	fileNum <- length(prefixList)
	for(file_index in 1:fileNum){
		prefix <- prefixList[file_index]
		freqFileName <- sprintf('%s.freq.txt', prefix)
		tempFreqDF  <- read.table(freqFileName);  tempFreqDF$file <- prefix
		if( file_index == 1 ){
			freqDF  <- tempFreqDF
		} else {
			freqDF <- rbind(freqDF, tempFreqDF)
		}
	}
	names(freqDF) <- c('freq', 'flux', 'file')
	return(freqDF)
}

velocRead <- function(prefixList){
	fileNum <- length(prefixList)
	spix <- numeric(fileNum)
	for(file_index in 1:fileNum){
		prefix <- prefixList[file_index]
		velocFileName <- sprintf('%s.veloc.txt', prefix)
		tempVelocDF <- read.table(velocFileName); tempVelocDF$file <- prefix
		if( file_index == 1 ){
			velocDF  <- tempVelocDF
		} else {
			velocDF <- rbind(velocDF, tempVelocDF)
		}
	}
	names(velocDF) <- c('veloc', 'flux', 'file')
	return(velocDF)
}

specAlign <- function(DF){
	fileList <- unique(DF$file)
	fileNum <- length(fileList)
	contFlux <- numeric(fileNum)
	centerFreq  <- numeric(fileNum)
	DF$weight <- 0.0
	DF$fileIndex <- 0
	#-------- Global spectral index
	for(file_index in 1:fileNum){
		blFreq <- get(fileList[file_index])
		DF[DF$file == fileList[file_index],]$fileIndex <- file_index
		blCHindex <- which( ((DF$freq >= blFreq[1] - 0.001) & (DF$freq <= blFreq[2] + 0.001)) |  ((DF$freq >= blFreq[3] - 0.001) & (DF$freq <= blFreq[4] + 0.001)))
		DF$weight[blCHindex] <- 1.0
		DF[DF$fileIndex != file_index,]$weight <- 0.0
		centerFreq[file_index] <- mean(DF[DF$fileIndex == file_index,]$freq)
		contFlux[file_index] <- sum(DF$flux* DF$weight) / sum(DF$weight)
	}
	fit <- lm( formula=logS~logF, data=data.frame( logF=log(centerFreq), logS=log(contFlux) ))
	SPIX <- coef(fit)[[2]]
	FluxRef <- exp(coef(fit)[[1]])
	DF$refSpec <- FluxRef *  exp(SPIX* log(DF$freq))
	DF$correction <- 0.0
	#-------- Local (SPW) alignment
	for(file_index in 1:fileNum){
		blFreq <- get(fileList[file_index])
		blCHindex <- which( ((DF$freq >= blFreq[1] - 0.001) & (DF$freq <= blFreq[2] + 0.001)) |  ((DF$freq >= blFreq[3] - 0.001) & (DF$freq <= blFreq[4] + 0.001)))
		DF$weight <- 0.0
		DF$weight[blCHindex] <- 1.0
		refFreq <- DF$freq - centerFreq[file_index]
		fit <- lm(formula=S~F, data=data.frame(F=refFreq, S=DF$flux - DF$refSpec), weight=DF$weight )
		DF[DF$fileIndex == file_index,]$correction <- coef(fit)[[1]] + coef(fit)[[2]] * (DF[DF$fileIndex == file_index,]$freq - centerFreq[file_index])
	}
	#-------- Correction
	DF$flux <- DF$flux - DF$correction
	DF$Tau <- -log(DF$flux / DF$refSpec)
	return(DF)	
}

Dopp <- function(FD, VD, fRestList){
	fileList <- unique(FD$file)
	fileNum <- length(fileList)
	D1 <- numeric(fileNum)
	D2 <- numeric(fileNum)
	for(file_index in 1:fileNum){
		index <- which(FD$file == fileList[file_index])
		fit <- lm(formula=y~x, data=data.frame(x=FD$freq[index]/fRestList[file_index], y=VD$veloc[index]))
		D1[file_index] <- coef(fit)[1]
		D2[file_index] <- coef(fit)[2]
	}
	return( matrix(c(D1, D2), ncol=2) )
}

readPlot <- function(fileList, fRestList, listofLines, labelY, plotRange, Title){
	cvel <- 299792.458		# Velocity of light, km/s
	freqDF  <- freqRead(fileList)
	velocDF <- velocRead(fileList)
	Dfacts <- Dopp(freqDF, velocDF, fRestList) - cvel
	Dofs <- median(Dfacts[,1])
	Dreg <- median(Dfacts[,2])		# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	DF <- specAlign(freqDF)
	fileNum <- length(fileList)
	plot( DF$freq, DF$flux, type='n', xlim=plotRange[1:2], ylim=plotRange[3:4], xlab='Topocentric Frequency [GHz]', ylab='Flux Density [Jy]', main = Title)
	for( file_index in 1:fileNum){
		index <- which( DF$file == fileList[file_index])
		chSep <- median(diff(DF$freq[index]))
		lines( DF$freq[index]-0.5*chSep, DF$flux[index], pch=20, col=file_index, type='s')
		points( DF$freq[index], DF$flux[index], pch=20, col='black', cex=0.1)
	}
	lines(DF$freq, DF$refSpec, col='grey')
	velPlot <- c(1200, 1500, 1800)
	index <- 1
	for(lineID in listofLines){
		plotLine(lineID, velPlot, Dofs, Dreg, labelY[index], 0.002)
		index <- index + 1
	}
	for(lineID in listofLines){
		lineDF <- plotTau(DF, lineID, c(1100,2000), Dofs, Dreg)
		lineFile <- sprintf('%s.data', gsub("[ |/]", "_", lineID$name))
		write.table(lineDF[order(lineDF$VLSR),], lineFile, row.names=F, col.names=T, quote=F)
	}
	return(DF)
}

source('FileBaseline.R')
#-------- Band3 LSB
pdf('NGC1052B6LSB.pdf')
DF <- readPlot(
	c('NGC_1052_B3SPW0.2016.1.00375.S', 'NGC_1052_B3SPW1.2016.1.00375.S'),		# fileList
	c(88.6318470, 86.0939500),													# fRestList
	list(HCN10v0, HCO10, SO2211, H13CN10),				# listofLines
	labelY <- c(1.20, 1.20, 1,20, 1.20),
	plotRange <- c(85, 89, 1.15, 1.23),
	Title <- 'NGC 1052 Band-3 LSB'
)
dev.off()
