prefix = 'uid___A002_Xc2bb44_X1baf'
os.system('rm -rf ' + prefix + '_flagonline.txt')
importasdm(prefix)
plotants( vis=prefix + '.ms', figfile=prefix + '_plotants.png') #REFANT = DA49
listobs(vis=prefix + '.ms', listfile=prefix + '.listobs')
#Fields: 4
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0238+1636          02:38:38.930100 +16.36.59.27500 ICRS    0        7346205
#  1    none J0243-0550          02:43:12.469466 -05.50.55.29615 ICRS    1        5059350
#  2    none J0245-1107          02:45:24.952113 -11.07.16.81108 ICRS    2        1583010
#  3    none NGC_1052            02:41:04.798510 -08.15.20.75170 ICRS    3       12757545
BPCAL = 'J0238+1636'
PHCAL = 'J0243-0550'
CHECK = 'J0245-1107'
TARGET= 'NGC_1052'
REFANT = 'DA53'
msfile = prefix + '.ms'
calMS = prefix + '.cal.ms'
gcalMS = prefix + '.gcal.ms'
flagdata(vis=prefix+'.ms', mode='manual', spw='5~32', autocorr=True, flagbackup=False)
flagmanager(vis=prefix+'.ms', mode='save', versionname='Apriori')
#-------- Tsys and WVR calibration table generation
os.system('rm -rf cal.tsys')
gencal(vis=prefix+'.ms', caltype='tsys', caltable='cal.tsys')
plotms('cal.tsys', xaxis = 'freq', yaxis = 'tsys', spw='17:4~124,19:4~124,21:4~124,23:4~124', iteraxis = 'antenna', coloraxis='corr')
# WVR
os.system('rm -rf cal.WVR')
wvrgcal(vis = msfile, caltable = 'cal.WVR', wvrspw=[4], segsource=True, toffset = 0, tie=[PHCAL + ',' + CHECK + ',' + BPCAL + ',' + TARGET], statsource = TARGET)
plotms('cal.WVR', xaxis='time',yaxis='phase',iteraxis='antenna')
from almahelpers_localcopy import tsysspwmap
tsysmap = tsysspwmap(vis=prefix+'.ms', tsystable='cal.tsys', tsysChanTol=1)
for FIELD in [PHCAL, TARGET, BPCAL]:
    applycal(vis=msfile, field=FIELD, flagbackup=False, spw='25,27,29,31', interp=['nearest', 'nearest'], gaintable=['cal.tsys', 'cal.WVR'], gainfield=FIELD, spwmap=[tsysmap, []], calwt=True)
#
applycal(vis=msfile, field=CHECK, flagbackup=False, spw='25,27,29,31', interp=['nearest', 'nearest'], gaintable=['cal.tsys', 'cal.WVR'], gainfield=[TARGET, CHECK], spwmap=[tsysmap, []], calwt=True)
flagdata(vis=prefix+'.ms', mode='manual', antenna='DA44', flagbackup=False)
os.system('rm -rf ' + calMS + '*')
split(vis=msfile, outputvis=calMS, datacolumn='corrected', spw='25,27,29,31')
#-------- Phase Cal for bandpass
os.system('rm -rf P0*')
gaincal(vis=calMS, caltable='P0', spw='*:4~478', field=BPCAL, selectdata=True, solint='int', refant=REFANT, calmode='p',gaintype='G',minsnr=4)
plotms('P0', xaxis = 'time', yaxis = 'phase', plotrange = [0,0,-180,180], iteraxis = 'antenna', coloraxis='spw')
#-------- Bandpass Cal
os.system('rm -rf B0*')
bandpass(vis=calMS, caltable='B0', gaintable='P0', field=BPCAL, scan='3', spw='*', minblperant=5, minsnr=5, solint=['inf','4ch'], bandtype='B', fillgaps=1, refant=REFANT, solnorm = True)
plotms('B0', xaxis='freq', yaxis='phase', plotrange = [0,0,-70,70], coloraxis='corr', iteraxis='antenna')
plotms('B0', xaxis='freq', yaxis='amp', plotrange = [0,0,0,1.2], coloraxis='corr', iteraxis='antenna')
#-------- Phase Cal for all
os.system('rm -rf P1*')
gaincal(vis=calMS, caltable='P1', spw='*:1~478', field=PHCAL + ',' + BPCAL + ',' + CHECK + ',' + TARGET, gaintable=['P0','B0'], gainfield=[BPCAL, BPCAL], interp=['nearest', 'nearest'],  selectdata=True, solint='int', refant=REFANT, gaintype='T', calmode='p', minsnr=4)
plotms(caltable='P1', xaxis='time', yaxis='phase', iteraxis='antenna')
#-------- Flux cal
#Provincia.local[kameno]4: Rscript ~/ALMA_SV/polQuery.R -D2017/07/23 -F97.5 J0238+1636
#Provincia.local[kameno]5: more CalQU.data 
#J0238+1636 1.04847824542899 -0.00643683165423746 0.00765566091768755
os.system('rm -rf G0*')
setjy( vis = calMS, field=BPCAL, spw='0,1,2,3', standard='manual', fluxdensity=[1.04847824542899, -0.00643683165423746, 0.00765566091768755, 0])
gaincal(vis = calMS, caltable = 'G0', spw ='*', field = PHCAL + ',' + CHECK + ',' + BPCAL + ',' + TARGET, minsnr=5.0, solint='10s', selectdata=True, solnorm=False, refant = REFANT, gaintable = ['P0', 'B0','P1'], gainfield=[BPCAL, BPCAL, ''], calmode = 'a', gaintype='T', spwmap=[[0,1,2,3],[0,1,2,3],[0,1,2,3]], interp=['nearest','nearest', 'nearest'])
plotms('G0', xaxis = 'time', yaxis = 'amp', plotrange = [], iteraxis = 'antenna')
fluxscale(vis=calMS, caltable='G0', fluxtable='G0.flux', reference=BPCAL, transfer=PHCAL + ',' + CHECK + ',' + TARGET, refspwmap=[0,1,2,3])
#2022-08-28 05:37:18 INFO fluxscale	##### Begin Task: fluxscale          #####
#2022-08-28 05:37:18 INFO fluxscale	fluxscale( vis='uid___A002_Xc2bb44_X1baf.cal.ms', caltable='G0', fluxtable='G0.flux', reference=['J0238+1636'], transfer=['J0243-0550,J0245-1107,NGC_1052'], listfile='', append=False, refspwmap=[0, 1, 2, 3], gainthreshold=-1.0, antenna='', timerange='', scan='', incremental=False, fitorder=1, display=False )
#2022-08-28 05:37:18 INFO fluxscale	****Using NEW VI2-driven calibrater tool****
#2022-08-28 05:37:18 INFO fluxscale	Opening MS: uid___A002_Xc2bb44_X1baf.cal.ms for calibration.
#2022-08-28 05:37:18 INFO fluxscale	Initializing nominal selection to the whole MS.
#2022-08-28 05:37:18 INFO fluxscale	Beginning fluxscale--(MSSelection version)-------
#2022-08-28 05:37:18 INFO fluxscale	 Found reference field(s): J0238+1636
#2022-08-28 05:37:18 INFO fluxscale	 Found transfer field(s):  J0243-0550 J0245-1107 NGC_1052
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0243-0550 in SpW=0 (freq=8.78549e+10 Hz) is: 0.190088 +/- 0.00592362 (SNR = 32.0898, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0243-0550 in SpW=1 (freq=8.59799e+10 Hz) is: 0.190654 +/- 0.00676397 (SNR = 28.1867, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0243-0550 in SpW=2 (freq=9.78969e+10 Hz) is: 0.186024 +/- 0.00575898 (SNR = 32.3015, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0243-0550 in SpW=3 (freq=9.95219e+10 Hz) is: 0.185622 +/- 0.00640307 (SNR = 28.9896, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0245-1107 in SpW=0 (freq=8.78549e+10 Hz) is: 0.0585376 +/- 0.000872263 (SNR = 67.11, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0245-1107 in SpW=1 (freq=8.59799e+10 Hz) is: 0.0588483 +/- 0.00112524 (SNR = 52.2983, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0245-1107 in SpW=2 (freq=9.78969e+10 Hz) is: 0.056845 +/- 0.00116526 (SNR = 48.7831, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for J0245-1107 in SpW=3 (freq=9.95219e+10 Hz) is: 0.0569467 +/- 0.00113513 (SNR = 50.1677, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for NGC_1052 in SpW=0 (freq=8.78549e+10 Hz) is: 1.1373 +/- 0.00260619 (SNR = 436.382, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for NGC_1052 in SpW=1 (freq=8.59799e+10 Hz) is: 1.13996 +/- 0.00293274 (SNR = 388.702, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for NGC_1052 in SpW=2 (freq=9.78969e+10 Hz) is: 1.12641 +/- 0.00280326 (SNR = 401.823, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Flux density for NGC_1052 in SpW=3 (freq=9.95219e+10 Hz) is: 1.12419 +/- 0.00299111 (SNR = 375.843, N = 44)
#2022-08-28 05:37:23 INFO fluxscale	 Fitted spectrum for J0243-0550 with fitorder=1: Flux density = 0.188085 +/- 6.96552e-05 (freq=92.6217 GHz) spidx: a_1 (spectral index) =-0.189754 +/- 0.00585914 covariance matrix for the fit:  covar(0,0)=0.000269791 covar(0,1)=-0.000195194 covar(1,0)=-0.000195194 covar(1,1)=0.358036
#2022-08-28 05:37:23 INFO fluxscale	 Fitted spectrum for J0245-1107 with fitorder=1: Flux density = 0.057791 +/- 7.60696e-05 (freq=92.6217 GHz) spidx: a_1 (spectral index) =-0.241216 +/- 0.0209266 covariance matrix for the fit:  covar(0,0)=8.48189e-05 covar(0,1)=0.000523992 covar(1,0)=0.000523992 covar(1,1)=0.113663
#2022-08-28 05:37:23 INFO fluxscale	 Fitted spectrum for NGC_1052 with fitorder=1: Flux density = 1.13194 +/- 0.000160409 (freq=92.6217 GHz) spidx: a_1 (spectral index) =-0.0927265 +/- 0.00223375 covariance matrix for the fit:  covar(0,0)=1.55677e-06 covar(0,1)=3.05587e-06 covar(1,0)=3.05587e-06 covar(1,1)=0.00205077
#2022-08-28 05:37:23 INFO fluxscale	Storing result in G0.flux
#2022-08-28 05:37:23 INFO fluxscale	Writing solutions to table: G0.flux
#2022-08-28 05:37:23 INFO fluxscale	Task fluxscale complete. Start time: 2022-08-28 01:37:18.066740 End time: 2022-08-28 01:37:23.221779
#2022-08-28 05:37:23 INFO fluxscale	##### End Task: fluxscale            #####
applycal(vis= calMS, flagbackup=False, spw='*', field='', interp=['nearest','nearest', 'nearest'], gainfield = [BPCAL,BPCAL, '',''], gaintable=['P0', 'B0', 'P1','G0.flux'], spwmap=[[0,1,2,3],[0,1,2,3],[0,1,2,3],[0,1,2,3]])
plotms(calMS, xaxis='time', yaxis='amp', ydatacolumn='corrected', coloraxis='corr', correlation='XX,YY', antenna='*&', field=TARGET, spw='0', avgchannel='480')
plotms(vis=calMS, xaxis='freq', yaxis='amp', ydatacolumn='corrected', selectdata=True, field=BPCAL, averagedata=True, avgchannel='', avgtime='10000', avgscan=True, avgbaseline=True, coloraxis='corr')
os.system('rm -rf ' + TARGET + '.ms*') 
os.system('rm -rf ' + PHCAL + '.ms*') 
os.system('rm -rf ' + BPCAL + '.ms*') 
os.system('rm -rf ' + CHECK + '.ms*') 
split(vis=calMS, outputvis=CHECK + '.ms', field=CHECK, datacolumn='corrected', spw='')
split(vis=calMS, outputvis=TARGET + '.ms', field=TARGET, datacolumn='corrected', spw='')
split(vis=calMS, outputvis=PHCAL + '.ms', field=PHCAL, datacolumn='corrected', spw='')
split(vis=calMS, outputvis=BPCAL + '.ms', field=BPCAL, datacolumn='corrected', spw='')
plotms(vis=BPCAL + '.ms', xaxis='freq', yaxis='phase', ydatacolumn='corrected', selectdata=True, field=BPCAL, averagedata=True, avgtime='10000', avgscan=True, avgbaseline=False, coloraxis='corr')
#-------- Imaging calibrator
os.system('rm -rf G3.'+BPCAL) 
gaincal(vis=BPCAL+'.ms', caltable='G3.'+BPCAL, field=BPCAL, selectdata=True, solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[1.04340668244966, 0.0, 0.0,0.0])
plotcal('G3.'+BPCAL, xaxis='time',yaxis='amp', spw='0,1,2,3', iteration='antenna', subplot=231)
applycal(vis= BPCAL + '.ms', flagbackup=False, spw='0,1,2,3',field=BPCAL, interp='nearest', gainfield = BPCAL, gaintable='G3.' + BPCAL)
os.system('rm -rf ' + BPCAL + '_SPW*') 
tclean(vis=BPCAL+'.ms', imagename=BPCAL + '_SPW0', field=BPCAL, spw='0', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='88.6318470GHz', niter=1000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # HCN J=1-0 v=0 88.6318470 GHz
imview(BPCAL + '_SPW0.image')
tclean(vis=BPCAL+'.ms', imagename=BPCAL + '_SPW1', field=BPCAL, spw='1', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='86.09395GHz', niter=1000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # SO 2(2) - 1(1) 86.0939500 GHz
#-------- Imaging phase calibrator
os.system('rm -rf G3.'+PHCAL) 
tclean(vis=PHCAL+'.ms', imagename=PHCAL + '_SPW0', field=PHCAL, spw='0', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='88.6318470GHz', niter=1000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # HCN J=1-0 v=0 88.6318470 GHz
os.system('rm -rf ' + PHCAL + '_SPW1*') 
tclean(vis=PHCAL+'.ms', imagename=PHCAL + '_SPW1', field=PHCAL, spw='1', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='86.0939500GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # SO 2(2) - 1(1) 86.0939500 GHz
os.system('rm -rf ' + PHCAL + '_SPW2*') 
tclean(vis=PHCAL+'.ms', imagename=PHCAL + '_SPW2', field=PHCAL, spw='2', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='97.98095GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # CS 2-1 97.98095 GHz
os.system('rm -rf ' + PHCAL + '_SPW3*') 
tclean(vis=PHCAL+'.ms', imagename=PHCAL + '_SPW3', field=PHCAL, spw='3', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='99.29987GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # SO J(N) = 3(2) - 2(1) v=0 99.29987 GHz
#-------- Imaging CHECK Source
os.system('rm -rf G3.'+CHECK) 
tclean(vis=CHECK+'.ms', imagename=CHECK + '_SPW0', field=CHECK, spw='0', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='88.6318470GHz', niter=1000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # HCN J=1-0 v=0 88.6318470 GHz
os.system('rm -rf ' + CHECK + '_SPW1*') 
tclean(vis=CHECK+'.ms', imagename=CHECK + '_SPW1', field=CHECK, spw='1', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='86.0939500GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # SO 2(2) - 1(1) 86.0939500 GHz
os.system('rm -rf ' + CHECK + '_SPW2*') 
tclean(vis=CHECK+'.ms', imagename=CHECK + '_SPW2', field=CHECK, spw='2', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='97.98095GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # CS 2-1 97.98095 GHz
os.system('rm -rf ' + CHECK + '_SPW3*') 
tclean(vis=CHECK+'.ms', imagename=CHECK + '_SPW3', field=CHECK, spw='3', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='99.29987GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=128, cell='0.05arcsec', weighting='natural') # SO J(N) = 3(2) - 2(1) v=0 99.29987 GHz
#-------- Target Continuum imaging
os.system('rm -rf G3.'+ TARGET + '_SPW*') 
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW0', field=TARGET, selectdata=True, spw='0', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[1.13730, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW0', xaxis='time',yaxis='amp', spw='0', iteraxis='antenna', coloraxis='corr')
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW1', field=TARGET, selectdata=True, spw='1', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[1.13996, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW1', xaxis='time',yaxis='amp', spw='1', iteraxis='antenna', coloraxis='corr')
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW2', field=TARGET, selectdata=True, spw='2', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[1.11854, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW2', xaxis='time',yaxis='amp', spw='2', iteraxis='antenna', coloraxis='corr')
gaincal(vis=TARGET+'.ms', caltable='G3.' + TARGET + '_SPW3', field=TARGET, selectdata=True, spw='3', solint='int,inf',refant=REFANT, gaintype='G', calmode='ap', smodel=[1.10656, 0.0, 0.0,0.0] )
plotms('G3.'+ TARGET + '_SPW3', xaxis='time',yaxis='amp', spw='3', iteraxis='antenna', coloraxis='corr')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='0',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW0')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='1',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW1')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='2',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW2')
applycal(vis= TARGET + '.ms', flagbackup=False, spw='3',field=TARGET, interp='nearest', gainfield = TARGET, gaintable='G3.' + TARGET + '_SPW3')
os.system('rm -rf ' + TARGET + '_B3_cont*') 
tclean(vis=TARGET + '.ms', imagename=TARGET + '_B3_cont', spw='*:1~478', specmode='mfs', nterms=1, niter=1000, threshold='0.0mJy', imsize=[768,384], gain=0.1, cell='0.05arcsec', weighting='natural', selectdata=True, interactive=True, pbcor=True)
os.system('rm -rf ' + TARGET + '_SPW0*') 
tclean(vis=TARGET+'.ms', imagename=TARGET + '_SPW0', field=TARGET, spw='0', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='88.6318470GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256,256], cell='0.05arcsec', weighting='natural') # HCN J=1-0 v=0 88.6318470 GHz
os.system('rm -rf ' + TARGET + '_SPW1*') 
tclean(vis=TARGET+'.ms', imagename=TARGET + '_SPW1', field=TARGET, spw='1', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='86.0939500GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256,256], cell='0.05arcsec', weighting='natural') # SO 2(2) - 1(1) 86.0939500 GHz
os.system('rm -rf ' + TARGET + '_SPW2*') 
tclean(vis=TARGET+'.ms', imagename=TARGET + '_SPW2', field=TARGET, spw='2', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='97.98095GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256,256], cell='0.05arcsec', weighting='natural') # CS 2-1 97.98095 GHz
os.system('rm -rf ' + TARGET + '_SPW3*') 
tclean(vis=TARGET+'.ms', imagename=TARGET + '_SPW3', field=TARGET, spw='3', specmode='cube', start=1, nchan=478, width=1, outframe='LSRK', veltype='radio', restfreq='99.29987GHz', niter=10000, gain=0.2, threshold='0.0mJy', interactive=True, selectdata=True, imsize=[256,256], cell='0.05arcsec', weighting='natural') # SO J(N) = 3(2) - 2(1) v=0 99.29987 GHz
imview(TARGET + '_SPW0.image')
imview(TARGET + '_SPW1.image')
imview(TARGET + '_SPW2.image')
imview(TARGET + '_SPW3.image')
