setwd('./')
source('LineList.R')
plotLine <- function(lineSpec, lineFreq, contLevel, height=0.005, Vsys=1492.0, Vofs=0.0, Vreg=0.0, bottom=0.0){
	cvel <- 299792.458		# Velocity of light, km/s
	velPlot <- seq(-300, 300, by=100)
	# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	#cvel <- 299792.458
	

	segments( lineFreq, contLevel + height*0.2, lineFreq, contLevel + height*0.4, lwd=2, col='red' )
	segments( lineFreq, contLevel + height*0.4, lineFreq + lineSpec$offset, contLevel + height*0.7, lwd=2, col='red' )
	segments( lineFreq + lineSpec$offset, contLevel + height*0.7, lineFreq + lineSpec$offset, contLevel + height*0.9, lwd=2, col='red' )
	text( lineFreq + lineSpec$offset, contLevel + height, lineSpec$label, pos=4, offset=0.1, cex=0.8, srt=90)
	# text( freq, bottom + height, lineSpec$name, pos=1, offset=0.2, cex=1, srt=90)
	# cat(freq)
	
	#if(bottom > 0.0){
	#	freq <- (Vsys - Vofs - cvel) * lineSpec$restfreq / (Vreg + cvel)
	#	top <- bottom + height
	#	segments( min(freq), bottom, max(freq), bottom )
	#	segments( freq, rep(bottom, length(freq)), freq, rep(top, length(freq)))
	#	for( tick_index in 1:length(velPlot)){
	#		if(velPlot[tick_index] == 0){
	#			text_sd <- 'Vsys'
	#		} else {
	#			text_sd <- sprintf(' %4.0f', velPlot[tick_index])
	#		}
	#		text( freq[tick_index], top, text_sd, adj=0, cex=0.3, srt=90)
	#	}
	#	#text( mean(freq), bottom, lineSpec$name, pos=1, offset=0.2, cex=0.3)
	#}
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
	text_sd <- sprintf('%s ATT=%f TauMax = %f at %f km/s Flux=%f Cont=%f\n', lineSpec$name, plotDF$ATT[peakIndex], plotDF$Tau[peakIndex], plotDF$VLSR[peakIndex], plotDF$flux[peakIndex], plotDF$refSpec[peakIndex])
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

SpliceOverlap <- function(DF, flagDF, fileIndexA, fileIndexB, flagIndex){
	fileList <- unique(DF$file)
	chSpc <- abs(median(diff(flagDF$freq))) 	# channel spacing
	overlap_min <- min(flagDF[flagDF$fileBit == flagIndex,]$freq)
	overlap_max <- max(flagDF[flagDF$fileBit == flagIndex,]$freq)
	cat(sprintf('Overlap %.1f - %.1f GHz\n', overlap_min, overlap_max))
	fileDFa <- DF[DF$file == fileList[fileIndexA],]
	fileDFb <- DF[DF$file == fileList[fileIndexB],]
	OffsetA <- mean(fileDFa[((fileDFa$freq > overlap_min) & (fileDFa$freq < overlap_max)),]$flux)
	OffsetB <- mean(fileDFb[((fileDFb$freq > overlap_min) & (fileDFb$freq < overlap_max)),]$flux)
	DF[DF$file == fileList[fileIndexA],]$flux <- DF[DF$file == fileList[fileIndexA],]$flux - 0.5*(OffsetA - OffsetB)
	DF[DF$file == fileList[fileIndexB],]$flux <- DF[DF$file == fileList[fileIndexB],]$flux - 0.5*(OffsetB - OffsetA)
	return(DF)
}
#-------- Alignment for multiple spectral files
specAlign <- function(DF){
	fileList <- unique(DF$file)
	fileNum <- length(fileList)
	contFlux <- numeric(fileNum)
	centerFreq  <- numeric(fileNum)
	# localSlope <- numeric(fileNum)
	DF$fileIndex <- 0
	#-------- Global spectral index
	for(file_index in 1:fileNum){
		blFreq <- get(fileList[file_index])
		DF[DF$file == fileList[file_index],]$fileIndex <- file_index
		blCHindex <- which(
				((DF$freq >= blFreq[1] - 0.001) & (DF$freq <= blFreq[2] + 0.001) & (DF$fileIndex == file_index)) |
				((DF$freq >= blFreq[3] - 0.001) & (DF$freq <= blFreq[4] + 0.001) & (DF$fileIndex == file_index)) )
		DF$weight <- 0.0; DF$weight[blCHindex] <- 1.0
		centerFreq[file_index] <- mean(DF[DF$fileIndex == file_index,]$freq)
		fit <- lm(formula=S~F, data=data.frame( S=DF$flux[blCHindex], F=DF$freq[blCHindex] - centerFreq[file_index]), weight=DF$weight[blCHindex] )
		contFlux[file_index] <- coef(fit)[[1]]
	}
	fit <- lm( formula=logS~logF, data=data.frame( logF=log(centerFreq), logS=log(contFlux) ))
	SPIX <- coef(fit)[[2]]
	FluxRef <- exp(coef(fit)[[1]])
	DF$refSpec <- FluxRef *  exp(SPIX* log(DF$freq))
	DF$correction <- 0.0
	#-------- Local (SPW) alignment
	for(file_index in 1:fileNum){
		blFreq <- get(fileList[file_index])
		blCHindex <- which(
				((DF$freq >= blFreq[1] - 0.001) & (DF$freq <= blFreq[2] + 0.001) & (DF$fileIndex == file_index)) |
				((DF$freq >= blFreq[3] - 0.001) & (DF$freq <= blFreq[4] + 0.001) & (DF$fileIndex == file_index)) )
		DF$weight <- 0.0; DF$weight[blCHindex] <- 1.0
		refFreq <- DF$freq - centerFreq[file_index]
		fit <- lm(formula=S~F, data=data.frame(F=refFreq, S=DF$flux - DF$refSpec), weight=DF$weight )
		DF[DF$fileIndex == file_index,]$correction <- coef(fit)[[1]] + coef(fit)[[2]] * (DF[DF$fileIndex == file_index,]$freq - centerFreq[file_index])
	}
	cat(sprintf('Std.Err. in FluxRef : %f\n', sd(DF$correction)))
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

readPlot <- function(fileList, refFileList, fRestList, listofLines, plotRange, Title, Vsys=1492.0){
	cvel <- 299792.458		# Velocity of light, km/s
	freqDF  <- freqRead(fileList)
	velocDF <- velocRead(fileList)
	Dfacts <- Dopp(freqDF, velocDF, fRestList) - cvel
	Dofs <- median(Dfacts[,1])
	Dreg <- median(Dfacts[,2])		# VLSR = Dofs + cvel + fobs/frest * (Dreg + cvel)
	DF <- specAlign(freqDF)
	freqRatio <- (Dreg + cvel)/(Vsys - Dofs - cvel)
	DF$restfreq <- DF$freq* freqRatio
    if(length(refFileList) > 0){
        refDF   <- freqRead(refFileList)
        refDF <- specAlign(refDF)
        refOffset <- min(DF$flux) - median(refDF$flux) - 0.01
        # cat(sprintf("Ref Offset = %f\n", refOffset))
    }
	fileNum <- length(fileList)
	#-------- Plot frame
	plot( DF$freq, DF$flux, type='n', axes = TRUE, xlim=plotRange[1:2], ylim=plotRange[3:4], xlab='Topocentric Frequency [GHz]', ylab='Flux Density [Jy]')
	#axis(1)
	#axis(2)
	#-------- Plot spectrum
	for( file_index in 1:fileNum){
		index <- which( DF$file == fileList[file_index])
		chSep <- median(diff(DF$freq[index]))
		lines( DF$freq[index]-0.5*chSep, DF$flux[index], pch=20, col='black', type='s')	# 'steps-mid'-like lines
        if(length(refFileList) > 0){
		    lines( refDF$freq[index]-0.5*chSep, refDF$flux[index] + refOffset, pch=20, col='cyan', type='s')	# Reference spectrum
		    points( DF$freq[index], DF$flux[index] + refOffset, pch=20, col='black', cex=0.2)			# data points
        }
		cat(sprintf('%.1f GHz : %.3f Jy\n', DF$freq[index[1]], DF$refSpec[index[1]]))
		cat(sprintf('%.1f GHz : %.3f Jy\n', DF$freq[index[length(index)]], DF$refSpec[index[length(index)]]))
	}
	#-------- Plot continuum level
	lines(DF$freq, DF$refSpec, col='grey')

	#-------- Plot continuum level
	#Vsys <- 1492.0
	for(lineID in listofLines){
		lineFreq <- (Vsys - Dofs - cvel) * lineID$restfreq / (Dreg + cvel)
		contLevel <- DF$refSpec[which.min(abs(DF$freq - lineFreq))]
		plotLine(lineID, lineFreq, contLevel+0.002, 0.005, Vsys, Dofs, Dreg, 0.45)
		
	}
	#-------- velocity scale
	lineID <- listofLines[[1]]
	velPlot <- c(-400, -300, -200, -100, 100, 200, 300, 400)
	minFreq <- 1e9; maxFreq <- -1e9
	lineFreq <- (Vsys - Dofs - cvel) * lineID$restfreq / (Dreg + cvel)
	contLevel <- DF$refSpec[which.min(abs(DF$freq - lineFreq))]
	for( velLine in velPlot){
		lineFreq <- (Vsys + velLine - Dofs - cvel) * lineID$restfreq / (Dreg + cvel)
		minFreq <- min(lineFreq, minFreq); maxFreq <- max(lineFreq, maxFreq)
		segments( lineFreq, contLevel + 0.004, lineFreq, contLevel + 0.005, lwd=1 )
		if( abs(velLine) > 300){
			text_vel <- sprintf('%+d', velLine)
			text( lineFreq, contLevel + 0.005, text_vel, pos=4, offset=0.1, cex=0.5, srt=90)
		}
	}
	segments( minFreq, contLevel + 0.004, maxFreq, contLevel + 0.004, lwd=1 )
	
	par(new = TRUE)
	plot( DF$restfreq, DF$flux, type='n', axes = FALSE, xlim=freqRatio* plotRange[1:2], ylim=plotRange[3:4], xlab='', ylab='', main = '')
	mtext("Rest Frequency [GHz]", side=3)
	axis(3)
	#title(main = Title, line=1)
	par(new = FALSE)
	
	for(lineID in listofLines){
		lineDF <- plotTau(DF, lineID, c(1100,2000), Dofs, Dreg)
		lineFile <- sprintf('%s.data', gsub("[ |/]", "_", lineID$name))
		write.table(lineDF[order(lineDF$VLSR),], lineFile, row.names=F, col.names=T, quote=F)
	}
	return(DF)
}

source('FileBaseline.R')
#-------- Band3 LSB
pdf('NGC1052B3LSB.pdf', width=10, height=7)
DF <- readPlot(
	c('NGC_1052_B3SPW0.2016.1.00375.S', 'NGC_1052_B3SPW1.2016.1.00375.S'),		# fileList
	c('J0243_0550_B3SPW0.2016.1.00375.S', 'J0243_0550_B3SPW1.2016.1.00375.S'),		# fileList
	c(88.6318470, 86.0939500),													# fRestList
	list(HCN10v0, HCO10, SO2211, H13CN10, SiO21v0, SiO21v1),				# listofLines
	plotRange <- c(85, 89, 1.15, 1.24),
	Title <- 'NGC 1052 Band-3 LSB',
	Vsys <- 1492.0
)
dev.off()

#-------- Band3 USB
pdf('NGC1052B3USB.pdf', width=10, height=7)
DF <- readPlot(
	c('NGC_1052_B3SPW2.2016.1.00375.S', 'NGC_1052_B3SPW3.2016.1.00375.S'),		# fileList
	c('J0243_0550_B3SPW2.2016.1.00375.S', 'J0243_0550_B3SPW3.2016.1.00375.S'),		# fileList
	c(97.98095, 99.29987),													# fRestList
	list(SO3221, SO4544, CS21),				# listofLines
	plotRange <- c(96.9, 100.45, 1.08, 1.133),
	Title <- 'NGC 1052 Band-3 USB'
)
dev.off()

#-------- Band4 LSB
pdf('NGC1052B4LSB.pdf', width=10, height=7)
DF <- readPlot(
	c('NGC_1052_B4SPW0.2016.1.00375.S', 'NGC_1052_B4SPW1.2016.1.00375.S'),		# fileList
	c('J0243_0550_B4SPW0.2016.1.00375.S', 'J0243_0550_B4SPW1.2016.1.00375.S'),		# fileList
	c(129.13892, 128.60513),													# fRestList
	list(SO2_12_12, SO3322, SO2_12_11, SO2_10_10, SiO32v0, SiO32v1),				# listofLines
	plotRange <- c(126.4, 130.1, 0.847, 0.905),
	Title <- 'NGC 1052 Band-4 LSB'
)
dev.off()

#-------- Band4 USB
pdf('NGC1052B4USB.pdf', width=10, height=7)
DF <- readPlot(
	c('NGC_1052_B4SPW2.2016.1.00375.S', 'NGC_1052_B4SPW3.2016.1.00375.S'),		# fileList
	c('J0243_0550_B4SPW2.2016.1.00375.S', 'J0243_0550_B4SPW3.2016.1.00375.S'),		# fileList
	c(139.76660, 140.83950),													# fRestList
	list(SO1312, SO1413, SO1514, SO1615, SO2_6_6),				# listofLines
	plotRange <- c(138.25, 141.65, 0.81, 0.855),
	Title <- 'NGC 1052 Band-4 USB'
)
dev.off()
#-------- Band7 HCN
pdf('NGC1052B7HCN.pdf', width=10, height=7)
DF <- readPlot(
	c('B7COSPW2', 'B7COSPW3', 'B7HCNSPW0', 'B7HCNSPW1'),		# fileList
	c(),		# fileList
	c(356.734223, 357.921987, 354.505473, 353.622753),													# fRestList
	list(HCN43, HCN_354460, H26alpha, HCO43, HCO_356549, HCN_356256, HCO_358242, SO2_359151, SO2_357963, SO2_357926, SO2_357892, SO2_357672, SO2_357581, SO2_357388, SO2_357241, SO2_357165),				# listofLines
	plotRange <- c(351, 358, 0.38, 0.465),
	Title <- 'NGC 1052 Band-7 HCN'
)
DF <- readPlot(
	c('B7COSPW0', 'B7COSPW1', 'B7HCNSPW2', 'B7HCNSPW3'),		# fileList
	c(),		# fileList
	c(345.79599, 344.916247, 342.88285, 340.71416),													# fRestList
	list(CS76, CO32, H13CN43, SO8877, H15CN43, SO_346528, SO7867, CN32a, CN32b, SO2_341276, SO2_341403, SO2_341674, SO_23_21),				# listofLines
	plotRange <- c(338, 345.2, 0.38, 0.465),
	Title <- 'NGC 1052 Band-7 HCN'
)
dev.off()

#-------- Band6 SO5_5-4_4
pdf('NGC1052B6SO.pdf', width=10, height=7)
DF <- readPlot(
	c('B6SPW0', 'B6SPW1', 'B6SPW2', 'B6SPW3'),		# fileList
	c(),		# fileList
	c(230.5380, 231.90093, 217.10498, 215.220653),													# fRestList
	list(SO54),				# listofLines
	plotRange <- c(213.2, 214.9, 0.47, 0.51),
	Title <- 'NGC 1052 Band-6 SO'
)
dev.off()
