source('GaussConv.R')
FWHM_Sigma <- 2*sqrt(2*log(2))
plotSpec <- function(TauFileList, labels, colors, plotRange, Vsys=1492){
    velRange <- c(-500, 500)
    GaussRange <- c(-242, 258)
    GaussVel <- seq(GaussRange[1], GaussRange[2], by=10)
    GaussDF <- data.frame(V=GaussVel, Tau=as.numeric(0))
    plot(velRange, plotRange, type='n', xlab='Velocity- Vsys [km/s]', ylab='Optical Depth')
    abline(h=0)
    abline(v=0, lty=2, col='gray')
    #-------- Each line transition
    for(file_index in 1:length(TauFileList)){
        fileName <- TauFileList[file_index]
        DF <- read.table(fileName, header=T)
        DF$V <- DF$VLSR - Vsys
        flagIndex <- which(abs(diff(DF$V)) < 0.5*abs(median(diff(DF$V)))) + 1
        if(length(flagIndex) > 0){ DF <- DF[-flagIndex,]}
        chSep <- median(diff(DF$V))
        DF$Tau <- DF$Tau - quantile(DF$Tau, 0.05)
        lines(c(DF$V-0.5*chSep, DF$V[nrow(DF)]+0.5*chSep), c(DF$Tau, DF$Tau[nrow(DF)]), col=colors[file_index], type='s')
        points(DF$V, DF$Tau, pch=20, col=colors[file_index], cex=0.5)
        GaussDF$Tau <- GaussDF$Tau + gaussConv(DF$V, DF$Tau, GaussDF$V, median(abs(diff(DF$V))))
    }
    GaussDF$Tau <- GaussDF$Tau / length(TauFileList)
    legend('topleft', legend=labels, col=colors, pch=20, lty=1)
    return(GaussDF)
}


B67FileList <- c(
    'HCN_J=4-3_v=0.data',       # HCN J=4-3
    'SO_J_N=8_9-7_8.data',      # SO 8_9 - 7_8
    'SO_J_N=8_8-7_7.data',      # SO 8_8 - 7_7
    'SO_J_N=8_7-7_6.data',      # SO 8_7 - 7_6
    'SO_5_5-4_4.data')          # SO 5_5-4_4
B34FileList <- c(
    'HCN_J=1-0.data',                 # HCN J=1-0
    'SO_J(N)_=_3(3)_-_2(2)_v=0.data', # SO 3_3-2_2
    'SO_J(N)_=_3(2)_-_2(1)_v=0.data', # SO 3_3-2_2
    'SO_J(N)_=_2(2)_-_1(1)_v=0.data') # SO 2_2-1_1

pdf('OpticalDepthPlot.pdf', width=7, height=5)
velRange <- c(-500, 500)
plotRange <- c(-0.01, 0.16)
#---- submm plot and fit
submmDF <- plotSpec(B67FileList, c('HCN 4-3', expression('SO 8'[9]*'-7'[8]), expression('SO 8'[8]*'-7'[7]), expression('SO 8'[7]*'-7'[6]), expression('SO 5'[5]*'-4'[4])), c('black', 'red', 'pink', 'magenta', 'orange'), plotRange, 1492)
submmFit <- nls(formula=Tau~amp* exp(-0.5*((V - centerV)/sigmaV)^2), data=submmDF, start=list(amp=max(submmDF$Tau), centerV=0, sigmaV=100))
text_sd <- sprintf('submm : %.4f +- %.4f |  %.1f +- %.1f |  %.1f +- %.1f\n', summary(submmFit)$coefficients[1], summary(submmFit)$coefficients[4], summary(submmFit)$coefficients[2], summary(submmFit)$coefficients[5], FWHM_Sigma*summary(submmFit)$coefficients[3], FWHM_Sigma*summary(submmFit)$coefficients[6])
cat(text_sd)

#---- mm plot and fit
plotRange <- c(-0.002, 0.03)
mmDF <- plotSpec(B34FileList, c('HCN 1-0', expression('SO 3'[3]*'-2'[2]), expression('SO 3'[2]*'-2'[1]), expression('SO 2'[2]*'-1'[1])), c('black', 'purple', 'blue', 'orange'), plotRange,1492)
mmFit <- nls(formula=Tau ~ amp1* exp(-0.5*((V - centerV1)/sigmaV1)^2)
                       + amp2* exp(-0.5*((V - centerV2)/sigmaV2)^2)
                       + amp3* exp(-0.5*((V - centerV3)/sigmaV3)^2)
                       + amp4* exp(-0.5*((V - centerV4)/sigmaV4)^2)
                       + amp5* exp(-0.5*((V - centerV5)/sigmaV5)^2), 
            data=mmDF, start=list(
                amp1=0.014, centerV1=168, sigmaV1=27,
                amp2=0.007, centerV2=108, sigmaV2=23,
                amp3=0.005, centerV3=0, sigmaV3=90,
                amp4=0.002, centerV4=-52, sigmaV4=70,
                amp5=0.003, centerV5=-112, sigmaV5=20))
for(index in 1:5){
    text_sd <- sprintf('mm %d : %.4f +- %.4f |  %.1f +- %.1f |  %.1f +- %.1f\n', index,
        summary(mmFit)$coefficients[index*3 - 2], summary(mmFit)$coefficients[(index + 5)*3 - 2],
        summary(mmFit)$coefficients[index*3 - 1], summary(mmFit)$coefficients[(index + 5)*3 - 1],
        FWHM_Sigma*summary(mmFit)$coefficients[index*3],     FWHM_Sigma*summary(mmFit)$coefficients[(index+5)*3])
    cat(text_sd)
}
dev.off()
#---- mm plot and fit
pdf('OpticalDepthComponents.pdf', width=7, height=7)
GaussRange <- c(-260, 260)
plotRange <- c(-0.002, 0.07)
lineColors <- c('magenta', 'orange', 'darkgreen', 'skyblue', 'blue')
plot(GaussRange, plotRange, type='n', xlab='Velocity- Vsys [km/s]', ylab='Mean Optical Depth')
lines( submmDF$V, predict(submmFit), col='navy')
lines( mmDF$V, predict(mmFit), col='red')
for(index in 1:5){
    compDF <- data.frame(V = mmDF$V, Tau=mmDF$Tau)
    compDF$Tau <- summary(mmFit)$coefficients[index*3 - 2]* exp(-0.5*((compDF$V - summary(mmFit)$coefficients[index*3 - 1])/summary(mmFit)$coefficients[index*3])^2)
    lines(compDF, col=lineColors[index], lty=2, lwd=2)
}
points(submmDF$V, submmDF$Tau, pch=20, col='navy', cex=0.5)
points(mmDF$V, mmDF$Tau, pch=20, col='red', cex=0.5)
legend('topleft', legend=c('sub-mm', 'mm'), col=c('navy','red'), pch=20, lty=1)
dev.off()
