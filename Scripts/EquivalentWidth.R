EquiWidth <- function(DF, baseLine=c(1180, 1270, 1750, 1970)){
	#-------- baseline range
	blCHindex <- which(
		((DF$VLSR >= baseLine[1]) & (DF$VLSR <= baseLine[2])) |
		((DF$VLSR >= baseLine[3]) & (DF$VLSR <= baseLine[4])) )
	#-------- line integration range
	lineCHindex <- which( (DF$VLSR >= baseLine[2]) & (DF$VLSR <= baseLine[3]) )
	#-------- baseline subtraction
	if( diff(range(DF$VLSR[blCHindex])) > abs(baseLine[3] - baseLine[2])){	#-------- double-side baseline
		BLfit <- lm( formula=Tau~V, data=data.frame( V=DF$VLSR[blCHindex] - 0.5*(baseLine[2] + baseLine[3]),  Tau=DF$Tau[blCHindex]) )	
		DF$TauC <- DF$Tau - (coef(BLfit)[1] + coef(BLfit)[2]* (DF$VLSR - 0.5*(baseLine[2] + baseLine[3])))
		BLerr <- coef(summary(BLfit))[1,2]
	} else {
		DF$TauC <- DF$Tau - mean(DF$Tau[blCHindex])
		BLerr <- sd(DF$Tau[blCHindex]) / sqrt(length(blCHindex) - 1)
	}
	#-------- peak optical depth
	peakIndex <- which.max(DF$Tau)
	text_sd <- sprintf('TauMax = %f (%f) at %f km/s ', DF$Tau[peakIndex], sd(DF$Tau[blCHindex]), DF$VLSR[peakIndex])
	cat(text_sd)
	#-------- Integrate optical depth
	lineWidth <- diff(DF$VLSR)
	range <- 1:length(lineCHindex)
	EW <- DF$TauC[lineCHindex[range]] %*% lineWidth[range]	# Equivalent Width
	ER <- BLerr* sum(lineWidth[range])	# Error caused by baseline uncertainty
	return(c(EW, ER))
}


setwd('/Volumes/SSD/CASAimaging/2016.1.00375.S/Spec')

fileList <- c(
	'HCN_J=1-0.data',
	'HCO+_J=1-0.data',
	'CS_J=2-1.data',
	'H13CN_J=1-0.data',
	'SO_J(N)_=_2(2)_-_1(1)_v=0.data',
	'SO_J(N)_=_3(2)_-_2(1)_v=0.data',
	'SO_J(N)_=_3(3)_-_2(2)_v=0.data',
	'SO_J(N)_=_4(5)_-_4(4)_v=0.data',
	'SO_J(N)_=_8(8)_-_7(8)_v=0.data',
	'SO_J(N)_=_13(13)_-_13(12)_v=0.data',
	'SO_J(N)_=_14(14)_-_14(13)_v=0.data',
	'SO_J(N)_=_15(15)_-_15(14)_v=0.data',
	'SO_J(N)_=_16(16)_-_16(15)_v=0.data',
	'SO_5_5-4_4.data',
	'SO_J_N=7_8-6_7.data',
	'SO_J_N=8_7-7_6.data',
	'SO_J_N=8_8-7_7.data',
	'SO_J_N=8_9-7_8.data',
	'SO2_6(2,4)-6(1,5)_v=0.data',
	'SO2_10(2,8)-10(1,9)_v=0.data',
	'SO2_12(1,11)-11(2,10)_v=0.data',
	'SO2_12(2,10)-12(1,11)_v=0.data')

for(fileName in fileList){
	DF <- read.table(fileName, header=T)
	EW <- EquiWidth(DF, c(1180, 1270, 1750, 1970) )
	text_sd <- sprintf("%s : EW = %.3f km/s  err = %.3f km/s\n", fileName, EW[1], EW[2])
	cat(text_sd)
}

#-------- Results
#> source('/Volumes/SSD/CASAimaging/2016.1.00375.S/Spec/EquivalentWidth.R')
#HCN_J=1-0.data : EW = 4.945 km/s  err = 0.066 km/s
#HCO+_J=1-0.data : EW = 1.832 km/s  err = 0.100 km/s
#CS_J=2-1.data : EW = 0.468 km/s  err = 0.131 km/s
#H13CN_J=1-0.data : EW = 0.640 km/s  err = 0.095 km/s
#SO_J(N)_=_2(2)_-_1(1)_v=0.data : EW = 1.282 km/s  err = 0.067 km/s
#SO_J(N)_=_3(2)_-_2(1)_v=0.data : EW = 2.328 km/s  err = 0.109 km/s
#SO_J(N)_=_3(3)_-_2(2)_v=0.data : EW = 3.018 km/s  err = 0.113 km/s
#SO_J(N)_=_4(5)_-_4(4)_v=0.data : EW = 0.718 km/s  err = 0.109 km/s
#SO_J(N)_=_8(8)_-_7(8)_v=0.data : EW = 0.231 km/s  err = 0.071 km/s
#SO_J(N)_=_13(13)_-_13(12)_v=0.data : EW = -0.140 km/s  err = 0.075 km/s
#SO_J(N)_=_14(14)_-_14(13)_v=0.data : EW = 0.150 km/s  err = 0.050 km/s
#SO_J(N)_=_15(15)_-_15(14)_v=0.data : EW = 0.138 km/s  err = 0.059 km/s
#SO_J(N)_=_16(16)_-_16(15)_v=0.data : EW = -0.151 km/s  err = 0.093 km/s
#SO2_6(2,4)-6(1,5)_v=0.data : EW = 0.025 km/s  err = 0.074 km/s
#SO2_10(2,8)-10(1,9)_v=0.data : EW = 0.597 km/s  err = 0.071 km/s
#SO2_12(1,11)-11(2,10)_v=0.data : EW = 2.845 km/s  err = 0.112 km/s
#SO2_12(2,10)-12(1,11)_v=0.data : EW = 0.648 km/s  err = 0.073 km/s