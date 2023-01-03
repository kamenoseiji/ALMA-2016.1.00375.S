prefix = 'uid___A002_Xc27dd7_X2e85'
'''
os.system('rm -rf ' + prefix + '_flagonline.txt')
importasdm(prefix)
plotants(vis=prefix + '.ms', figfile=prefix + '_plotants.png') #REFANT = DA41
listobs(vis=prefix + '.ms', listfile=prefix + '.listobs')
#Fields: 3
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0238+1636          02:38:38.930100 +16.36.59.27500 ICRS    0        7136072
#  1    none J0243-0550          02:43:12.469466 -05.50.55.29615 ICRS    1        3323086
#  2    none NGC_1052            02:41:04.798510 -08.15.20.75170 ICRS    2       28480670
'''
BPCAL = 'J0238+1636'
PHCAL = 'J0243-0550'
TARGET= 'NGC_1052'
REFANT = 'DA41'
msfile = prefix + '.ms'
calMS = prefix + '.cal.ms'
'''
flagdata(vis=prefix+'.ms', mode='manual', spw='5~32', autocorr=True, flagbackup=False)
flagmanager(vis=prefix+'.ms', mode='save', versionname='Apriori')
#-------- Tsys and WVR calibration table generation
os.system('rm -rf cal.tsys')
gencal(vis=prefix+'.ms', caltype='tsys', caltable='cal.tsys')
plotms('cal.tsys', xaxis = 'freq', yaxis = 'tsys', spw='17:4~124,19:4~124,21:4~124,23:4~124', iteraxis = 'antenna')
# WVR
os.system('rm -rf cal.WVR')
wvrgcal(vis = msfile, caltable = 'cal.WVR', wvrspw=[4], segsource=True, toffset = 0, tie=[PHCAL + ',' + BPCAL + ',' + TARGET], statsource = TARGET)
plotms('cal.WVR', xaxis='time',yaxis='phase',iteraxis='antenna', coloraxis='spw')
from almahelpers_localcopy import tsysspwmap
tsysmap = tsysspwmap(vis=prefix+'.ms', tsystable='cal.tsys', tsysChanTol=1)
for FIELD in [PHCAL, TARGET, BPCAL]:
    applycal(vis=msfile, field=FIELD, flagbackup=False, spw='25,27,29,31', interp=['nearest', 'nearest'], gaintable=['cal.tsys', 'cal.WVR'], gainfield=FIELD, spwmap=[tsysmap, []], calwt=True)
#
os.system('rm -rf ' + calMS + '*')
split(vis=msfile, outputvis=calMS, datacolumn='corrected', spw='25,27,29,31')
#-------- Phase Cal for bandpass
os.system('rm -rf P0*')
gaincal(vis=calMS, caltable='P0', spw='*:4~478', field=BPCAL, selectdata=True, solint='int', refant=REFANT, calmode='p',gaintype='G',minsnr=4)
plotms('P0', xaxis = 'time', yaxis = 'phase', coloraxis='spw', iteraxis = 'antenna')
#-------- Bandpass Cal
os.system('rm -rf B0*')
bandpass(vis=calMS, caltable='B0', gaintable='P0', field=BPCAL, scan='3', spw='*', minblperant=5, minsnr=5, solint=['inf','4ch'], bandtype='B', fillgaps=1, refant=REFANT, solnorm = True)
plotms('B0', xaxis='freq', yaxis='phase', plotrange = [0,0,-70,70], iteraxis='antenna', coloraxis='corr')
plotms('B0', xaxis='freq', yaxis='amp', plotrange = [0,0,0,1.2], iteraxis='antenna', coloraxis='corr')
#-------- Phase Cal for all
os.system('rm -rf P1*')
gaincal(vis=calMS, caltable='P1', spw='*:1~478', field=PHCAL + ',' + BPCAL + ',' + TARGET, gaintable=['P0','B0'], gainfield=[BPCAL, BPCAL], interp=['nearest', 'nearest'],  selectdata=True, solint='int', refant=REFANT, gaintype='T', calmode='p', minsnr=4)
plotms('P1', xaxis='time', yaxis='phase', iteraxis='antenna')
#-------- Flux cal
#Provincia.local[kameno]4: Rscript ~/ALMA_SV/polQuery.R -D2017/07/23 -F138 J0238+1636
#Provincia.local[kameno]5: more CalQU.data 
#J0238+1636 0.856487649877853 -0.0100622687606874 0.00193140504874338
os.system('rm -rf G0*')
setjy( vis = calMS, field=BPCAL, spw='0,1,2,3', standard='manual', fluxdensity=[0.856487649877853, -0.0100622687606874, 0.00193140504874338, 0.0], spix = [-0.536715,0], reffreq = '138.0GHz', usescratch=False)
gaincal(vis = calMS, caltable = 'G0', spw ='*', field = PHCAL + ',' + BPCAL + ',' + TARGET, minsnr=5.0, solint='10s', selectdata=True, solnorm=False, refant = REFANT, gaintable = ['P0', 'B0','P1'], gainfield=[BPCAL, BPCAL, ''], calmode = 'a', gaintype='T', spwmap=[[0,1,2,3],[0,1,2,3],[0,1,2,3]], interp=['nearest','nearest', 'nearest'])
plotms('G0', xaxis = 'time', yaxis = 'amp', plotrange = [], iteraxis = 'antenna')
fluxscale(vis=calMS, caltable='G0', fluxtable='G0.flux', reference=BPCAL, transfer=PHCAL + ',' + TARGET, refspwmap=[0,1,2,3])
#2022-08-28 02:36:20 INFO fluxscale	##########################################
#2022-08-28 02:36:20 INFO fluxscale	##### Begin Task: fluxscale          #####
#2022-08-28 02:36:20 INFO fluxscale	fluxscale( vis='uid___A002_Xc27dd7_X2e85.cal.ms', caltable='G0', fluxtable='G0.flux', reference=['J0238+1636'], transfer=['J0243-0550,NGC_1052'], listfile='', append=False, refspwmap=[0, 1, 2, 3], gainthreshold=-1.0, antenna='', timerange='', scan='', incremental=False, fitorder=1, display=False )
#2022-08-28 02:36:20 INFO fluxscale	****Using NEW VI2-driven calibrater tool****
#2022-08-28 02:36:20 INFO fluxscale	Opening MS: uid___A002_Xc27dd7_X2e85.cal.ms for calibration.
#2022-08-28 02:36:20 INFO fluxscale	Initializing nominal selection to the whole MS.
#2022-08-28 02:36:20 INFO fluxscale	Beginning fluxscale--(MSSelection version)-------
#2022-08-28 02:36:20 INFO fluxscale	 Found reference field(s): J0238+1636
#2022-08-28 02:36:20 INFO fluxscale	 Found transfer field(s):  J0243-0550 NGC_1052
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for J0243-0550 in SpW=0 (freq=1.29149e+11 Hz) is: 0.158234 +/- 0.000962423 (SNR = 164.413, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for J0243-0550 in SpW=1 (freq=1.27384e+11 Hz) is: 0.15891 +/- 0.0012842 (SNR = 123.743, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for J0243-0550 in SpW=2 (freq=1.39191e+11 Hz) is: 0.1553 +/- 0.00111506 (SNR = 139.275, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for J0243-0550 in SpW=3 (freq=1.40819e+11 Hz) is: 0.154668 +/- 0.00137988 (SNR = 112.088, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for NGC_1052 in SpW=0 (freq=1.29149e+11 Hz) is: 1.00751 +/- 0.000865041 (SNR = 1164.7, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for NGC_1052 in SpW=1 (freq=1.27384e+11 Hz) is: 1.00947 +/- 0.000993298 (SNR = 1016.28, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for NGC_1052 in SpW=2 (freq=1.39191e+11 Hz) is: 1.00231 +/- 0.000969366 (SNR = 1033.99, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Flux density for NGC_1052 in SpW=3 (freq=1.40819e+11 Hz) is: 1.00165 +/- 0.00106519 (SNR = 940.352, N = 46)
#2022-08-28 02:36:26 INFO fluxscale	 Fitted spectrum for J0243-0550 with fitorder=1: Flux density = 0.156768 +/- 4.70999e-05 (freq=134.004 GHz) spidx: a_1 (spectral index) =-0.260377 +/- 0.00704833 covariance matrix for the fit:  covar(0,0)=1.36916e-05 covar(0,1)=9.66057e-05 covar(1,0)=9.66057e-05 covar(1,1)=0.0399514
#2022-08-28 02:36:26 INFO fluxscale	 Fitted spectrum for NGC_1052 with fitorder=1: Flux density = 1.00519 +/- 0.00024826 (freq=134.004 GHz) spidx: a_1 (spectral index) =-0.0741322 +/- 0.00567801 covariance matrix for the fit:  covar(0,0)=2.32364e-07 covar(0,1)=1.16805e-06 covar(1,0)=1.16805e-06 covar(1,1)=0.000651145
#2022-08-28 02:36:27 INFO fluxscale	Storing result in G0.flux
#2022-08-28 02:36:27 INFO fluxscale	Writing solutions to table: G0.flux
#2022-08-28 02:36:28 INFO fluxscale	Task fluxscale complete. Start time: 2022-08-27 22:36:19.663617 End time: 2022-08-27 22:36:27.522811
#2022-08-28 02:36:28 INFO fluxscale	##### End Task: fluxscale            #####
applycal(vis= calMS, flagbackup=False, spw='*', field='', interp=['nearest','nearest', 'nearest'], gainfield = [BPCAL,BPCAL, '',''], gaintable=['P0', 'B0', 'P1','G0.flux'], spwmap=[[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3]])
plotms(calMS, xaxis='time', yaxis='amp', ydatacolumn='corrected', coloraxis='corr', correlation='XX,YY', antenna='*&', field=TARGET, spw='0', avgchannel='480')
plotms(vis=calMS, xaxis='freq', yaxis='amp', ydatacolumn='corrected', selectdata=True, field=BPCAL, averagedata=True, avgchannel='', avgtime='10000', avgscan=True, avgbaseline=True, coloraxis='corr')
os.system('rm -rf ' + TARGET + '.ms*') 
os.system('rm -rf ' + PHCAL + '.ms*') 
os.system('rm -rf ' + BPCAL + '.ms*') 
split(vis=calMS, outputvis=TARGET + '.ms', field=TARGET, datacolumn='corrected', spw='')
split(vis=calMS, outputvis=PHCAL + '.ms', field=PHCAL, datacolumn='corrected', spw='')
split(vis=calMS, outputvis=BPCAL + '.ms', field=BPCAL, datacolumn='corrected', spw='')
plotms(vis=BPCAL + '.ms', xaxis='freq', yaxis='phase', ydatacolumn='corrected', selectdata=True, field=BPCAL, averagedata=True, avgtime='10000', avgscan=True, avgbaseline=False, coloraxis='corr')
#-------- Imaging calibrator
os.system('rm -rf ' + BPCAL + '_SPW0.*')
tclean(vis=BPCAL + '.ms', imagename=BPCAL + '_SPW0', spw='0', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='129.13892GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=256, cell='0.05arcsec', weighting='natural', pbcor=False)
imview(BPCAL + '_SPW0.image')
#-------- Imaging target continuum
os.system('rm -rf ' + TARGET + '_cont.*')
os.system('rm -rf G3.'+ TARGET + '_SPW*') 
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW0', field=TARGET, selectdata=True, spw='0', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[0.88, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW0', xaxis='time',yaxis='amp', spw='0', iteraxis='antenna', coloraxis='corr')
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW1', field=TARGET, selectdata=True, spw='1', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[0.89, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW1', xaxis='time',yaxis='amp', spw='1', iteraxis='antenna', coloraxis='corr')
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW2', field=TARGET, selectdata=True, spw='2', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[0.84, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW2', xaxis='time',yaxis='amp', spw='2', iteraxis='antenna', coloraxis='corr')
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW3', field=TARGET, selectdata=True, spw='3', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[0.83, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW3', xaxis='time',yaxis='amp', spw='3', iteraxis='antenna', coloraxis='corr')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='0',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW0')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='1',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW1')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='2',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW2')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='3',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW3')
os.system('rm -rf ' + TARGET + '_B4_cont*') 
tclean(vis=TARGET + '.ms', imagename=TARGET + '_B4_cont', spw='*:1~478', specmode='mfs', nterms=1, niter=1000, threshold='0.0mJy', imsize=[768,384], gain=0.1, cell='0.05arcsec', weighting='natural', selectdata=True, interactive=True, pbcor=True)
#-------- Imaging PHCAL
os.system('rm -rf ' + PHCAL + '_SPW0.*')
tclean(vis=PHCAL + '.ms', imagename=PHCAL + '_SPW0', spw='0', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='129.13892GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural', pbcor=False) # (SO 3_3 - 2_2)
os.system('rm -rf ' + PHCAL + '_SPW1.*')
tclean(vis=PHCAL + '.ms', imagename=PHCAL + '_SPW1', spw='1', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='128.60513GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural', pbcor=False) # SO2 12(2,10)-12(1,11) v=0
os.system('rm -rf ' + PHCAL + '_SPW2.*')
tclean(vis=PHCAL + '.ms', imagename=PHCAL + '_SPW2', spw='2', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='139.76660GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural', pbcor=False) # NH2CN 7(6,1) - 6(6.0)
os.system('rm -rf ' + PHCAL + '_SPW3.*')
tclean(vis=PHCAL + '.ms', imagename=PHCAL + '_SPW3', spw='3', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='140.83950GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural', pbcor=False) # H2CO 2(1,2) - 1(1,1)
imview(PHCAL + '_SPW0.image')
imview(PHCAL + '_SPW1.image')
imview(PHCAL + '_SPW2.image')
imview(PHCAL + '_SPW3.image')
#-------- Imaging Target
os.system('rm -rf ' + TARGET + '_SPW0.*')
tclean(vis=TARGET + '.ms', imagename=TARGET + '_SPW0', spw='0', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='129.13892GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256, 256], cell='0.05arcsec', weighting='natural', pbcor=False) # (SO 3_3 - 2_2)
os.system('rm -rf ' + TARGET + '_SPW1.*')
tclean(vis=TARGET + '.ms', imagename=TARGET + '_SPW1', spw='1', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='128.60513GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256, 256], cell='0.05arcsec', weighting='natural', pbcor=False) # SO2 12(2,10)-12(1,11) v=0
os.system('rm -rf ' + TARGET + '_SPW2.*')
tclean(vis=TARGET + '.ms', imagename=TARGET + '_SPW2', spw='2', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='139.76660GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256, 256], cell='0.05arcsec', weighting='natural', pbcor=False) # NH2CN 7(6,1) - 6(6.0)
os.system('rm -rf ' + TARGET + '_SPW3.*')
tclean(vis=TARGET + '.ms', imagename=TARGET + '_SPW3', spw='3', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='140.83950GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256, 256], cell='0.05arcsec', weighting='natural', pbcor=False) # H2CO 2(1,2) - 1(1,1)

imview(TARGET + '_SPW0.image')
imview(TARGET + '_SPW1.image')
imview(TARGET + '_SPW2.image')
imview(TARGET + '_SPW3.image')
'''
