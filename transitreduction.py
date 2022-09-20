import phot_pipe as pp 
from astropy.coordinates import SkyCoord, Angle, EarthLocation
from astropy import units as u
from quick_image import read_image, display_image
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

targname = 'TIC70899085.01'     # Add the transit identifier
obstype = 'Transit'     # Always gonna be transit

paths = pp.get_paths(targname=targname,obstype=obstype)     # Getting paths to files

dates = pp.get_dates(targname=targname)     # Gets all dates for transit

date = dates[0] # Choose the desired date

# XO-2b
#rastr  = '07 48 06.5'
#decstr = '50 13 32.92'

# TIC398733009
rastr = '04:16:45.6'
decstr =  '-12:05:02.4'
midtrans = 2458786.860

# WASP 33
#rastr  = '05:25:35.66'
#decstr = '08:26:27.2'


coords = SkyCoord(rastr,decstr,unit=(u.hour,u.deg))
ra  = coords.ra.deg
dec = coords.dec.deg

midtrans = 2458775.841

# Chose the proper code for the band below

# Get reference stars from PanSTARRS catalog
refs = pp.get_panstarrs_refs(ra=rastr,dec=decstr,maxMag=15)

#ifiles,ict = pp.get_files(prefix=targname,tag='i_',date=date,suffix='solved.fits',raw=False,clean=False)
#rfiles,rct = pp.get_files(prefix=ftargname,tag='r_',date=date,suffix='solved.fits',raw=False,clean=False)
gfiles,gct = pp.get_files(prefix=targname,tag='g_',date=date,suffix='solved.fits',raw=False,clean=False)

#iinfo = pp.check_headers(ifiles,rastr=rastr,decstr=decstr,dosex=False)
#rinfo = pp.check_headers(rfiles,rastr=rastr,decstr=decstr,dosex=False)
ginfo = pp.check_headers(gfiles,rastr=rastr,decstr=decstr,dosex=False)

#dists = np.sqrt((iinfo['xpos'][0]-iinfo['xpos'])**2+(iinfo['ypos'][0]-iinfo['ypos'])**2)
#inds = np.array([0,4,18])
#dists = np.sqrt((iinfo['xpos']-1024)**2 + (iinfo['ypos']-1024)**2)
#angles = np.degrees(np.arctan(iinfo['ypos']-1024,iinfo['xpos']-1024))
#inds = np.array([np.argmax(dists),np.argmin(dists),np.argsort(dists)[ict/2]])

band = 'g'
files = gfiles
ct = gct

# Get calibration frames
cals = pp.get_cal_frames(date,band=band,targname=targname,obstype=obstype,archive=True)

# Make straight stack to check fringeing

'''
fringestack = []
for i in range(len(files)):
    f = files[i]
    im,h = read_image(f)
    cal = pp.calibrate_image(im,h,cals,rotated_flat=False)
    if len(fringestack) == 0:
        fringestack = cal.reshape((1,2048,2048))
    else:
        fringestack = np.append(fringestack,cal.reshape((1,2048,2048)),axis=0)
fringe = np.amin(fringestack,axis=0)

'''
# Stack images for a deeper image
stack,hstack = pp.all_stack(files,ra=ra,dec=dec,outname=targname+'_'+band+'_stacked',fwhm_max=4.0,
                            outpath=paths['output'],calarchive=False,dtheta_max=30,dpix_max=150,
                            obstype=obstype,targname=targname)


im, h = read_image(files[-1])
cal = pp.calibrate_image(im,h,cals,rotated_flat=False)
from photutils import Background2D, MedianBackground

bkg_estimator = MedianBackground()
bkg = Background2D(cal, box_size=(7, 7), filter_size=(15, 15), bkg_estimator=bkg_estimator)

display_image(cal,fignum=5)
plt.figure(6)
plt.imshow(bkg.background, cmap='Greys',origin='lower')
plt.colorbar()
display_image(cal-bkg.background,fignum=7)


targ_info = pp.make_targ_info(refs,ra=rastr,dec=decstr)

pp.display_targets(stack,hstack,targ_info,targname=targname,obstype=obstype)

ras, decs = pp.choose_refs(files[len(files)/2],rastr,decstr,bias=cals['bias'],dark=cals['dark'],
                           flat=cals['flat'],outdir=paths['output'],outfile='choose_refs_coords4.txt')

targ_info['RAdeg'] = np.append(targ_info['RAdeg'],ras[1:])
targ_info['DECdeg'] = np.append(targ_info['DECdeg'],decs[1:])
targ_info['gmag'] = np.append(targ_info['gmag'],np.nan*np.zeros(len(ras)-1))
targ_info['rmag'] = np.append(targ_info['rmag'],np.nan*np.zeros(len(ras)-1))
targ_info['imag'] = np.append(targ_info['imag'],np.nan*np.zeros(len(ras)-1))
targ_info['zmag'] = np.append(targ_info['zmag'],np.nan*np.zeros(len(ras)-1))

pp.display_targets(stack,hstack,targ_info,targname=targname,obstype=obstype)

# Ref 2 and 4 too close to edge
for key in targ_info.keys():
    targ_info[key] = np.delete(targ_info[key],4)


pp.check_apertures(stack,hstack,targ_info,index=0,aperture=15,skyrad=[30,40])
pp.check_apertures(stack,hstack,targ_info,index=1,aperture=25,skyrad=[35,45])
pp.check_apertures(stack,hstack,targ_info,index=2,aperture=25,skyrad=[35,45])
pp.check_apertures(stack,hstack,targ_info,index=3,aperture=17,skyrad=[35,45])
pp.check_apertures(stack,hstack,targ_info,index=4,aperture=25,skyrad=[35,45])
pp.check_apertures(stack,hstack,targ_info,index=5,aperture=25,skyrad=[35,45])


pp.display_targets(stack,hstack,targ_info,targname=targname,obstype=obstype)

pp.write_targ_info(targ_info,obstype=obstype,targname=targname)


targ_info = pp.read_targ_info(obstype=obstype,targname=targname)
pp.display_targets(stack,hstack,targ_info,targname=targname,obstype=obstype,write=True,
                   outfile=targname+'_targets.png')

phot = pp.night_phot(date,obstype=obstype,targname=targname,targ_info=targ_info,
                     rotated_flat=False,tagpre='-',clobber=False,
                     targRA=ra,targDec=dec,pamax=15,distmax=150,saturation=30000.0)

#phot_rot = pp.night_phot(date,obstype=obstype,targname=targname,targ_info=targ_info,
#                         rotated_flat=True,tagpre='_',clobber=False,phot_tag='_rot',
#                         targRA=ra,targDec=dec)


plt.ion()
plt.figure(1)
plt.clf()
band = 'g'
targ = []
bjd  = []
airmass = []
saturated = []
for obs in phot[band][targname]:
    targ = np.append(targ,obs['flux'][2])
    bjd = np.append(bjd,obs['bjd'])
    airmass = np.append(airmass,obs['airmass'])
    if 'saturation' in obs['flags']:
        sat = True
    else:
        sat = False    
    saturated= np.append(saturated,sat)
plt.plot(bjd,targ/np.median(targ),'o')


refflux = {}
i = 1
for k in phot[band].keys():
    if 'ref' in k:        
        lc = []
        bjd = []
        for obs in phot[band][k]:
            lc = np.append(lc,obs['flux'][2])
            bjd = np.append(bjd,obs['bjd'])
        plt.plot(bjd,lc/np.median(lc)+i*0.3,'o')
        refflux[k] = lc
        i +=1
plt.xlabel('BJD')
plt.ylabel('Relative Flux + offset')
plt.title('Raw Photometry of '+targname+' in '+band)
plt.savefig(paths['output']+date+'/RawPhot_'+band+'.png',dpi=300)


relflux = []
for i in range(len(targ)):
    reftot = 0.0
    for k in refflux.keys(): #['ref1','ref9','ref10']:
        reftot += refflux[k][i]
    relflux = np.append(relflux,targ[i]/reftot)
                        

#satargs, = np.where(saturated == 1)

plt.figure(2)
plt.clf()
reltime = (bjd-midtrans)*24
args = np.argsort(bjd)
reltime = reltime[args]
relflux = relflux[args]
airmass = airmass[args]
relflux /= np.median(np.append(relflux[0:25],relflux[-0:]))
plt.plot(reltime,relflux,'o')#,label='Not Saturated')
#plt.plot(reltime[satargs],relflux[satargs],'o',label='Saturated')
#plt.legend(numpoints=1,loc='lower left')
plt.xlabel('Time from Mid-Transit (hrs)')
plt.ylabel('Relative Flux')
plt.title('Transit Photometry of '+targname+' in '+band+'$^\prime$: '+date)
plt.axhline(y=1.0,linestyle='--',color='gray')
plt.axvline(x=0.0,linestyle='--',color='gray')
plt.savefig(paths['output']+date+'/DiffPhot_'+band+'.png',dpi=300)

std = np.std(np.append(relflux[0:10],relflux[-1:]),ddof=1)
# std = 0.5 ppt


plt.ion()
plt.figure(3)
plt.clf()
targ_g = []
bjd_g  = []
airmass_g  = []
for obs in phot['g'][targname]:
    targ_g = np.append(targ_g,obs['flux'][2])
    bjd_g = np.append(bjd_g,obs['bjd'])
    airmass_g = np.append(airmass_g,obs['airmass'])
plt.plot(bjd_g,targ_g/np.median(targ_g),'o')

refflux_g = {}
i = 1
for k in phot['g'].keys():
    if 'ref' in k:        
        lc = []
        for obs in phot['g'][k]:
            lc = np.append(lc,obs['flux'][2])
        plt.plot(bjd_g,lc/np.median(lc)+i*0.15,'o')
        refflux_g[k] = lc
        i +=1
plt.xlabel('BJD')
plt.ylabel('Relative Flux + offset')
plt.title('Raw Photometry of '+targname+" in $g^\prime$ band")
plt.savefig(paths['output']+date+'/RawPhot_g.png',dpi=300)


relflux_g = []
for i in range(len(targ_g)):
    reftot = 0.0
    for k in refflux_g.keys():
        reftot += refflux_g[k][i]
    relflux_g = np.append(relflux_g,targ_g[i]/reftot)

# outlier

arg = np.argmax(relflux_g)
relflux_g = np.delete(relflux_g,arg)
bjd_g = np.delete(bjd_g,arg)
airmass_g = np.delete(airmass_g,arg)

plt.figure(4)
plt.clf()
midtrans = 2458786.860
 # from swarthmore site
reltime_g = (bjd_g-midtrans)*24
args = np.argsort(reltime_g)
reltime_g = reltime_g[args]
relflux_g = relflux_g[args]
airmass_g = airmass_g[args]
relflux_g /= np.median(np.append(relflux_g[0:10],relflux_g[-15:]))
plt.plot(reltime_g,relflux_g,'o')
plt.plot(np.append(reltime_g[0:19],reltime_g[-15:]),np.append(relflux_g[0:19],relflux_g[-15:]),'o')

plt.xlabel('Time from Mid-Transit (hrs)')
plt.ylabel('Relative Flux')
plt.title('Transit Photometry of '+targname+' in $g^\prime$: '+date)
plt.axhline(y=1.0,linestyle='--',color='gray')
plt.axvline(x=0.0,linestyle='--',color='gray')
plt.savefig(paths['output']+date+'/DiffPhot_g.png',dpi=300)

std = np.std(np.append(relflux_g[0:10],relflux_g[-15:]),ddof=1)






plt.figure(5)
plt.clf()
plt.plot(reltime,relflux,'o')
#plt.plot(reltime_g,relflux_g,'o')

plt.xlabel('Time from Mid-Transit (hrs)')
plt.ylabel('Relative Flux')
plt.title('Transit Photometry of '+targname+': '+date)
plt.axhline(y=1.0,linestyle='--',color='gray')
plt.axvline(x=0.0,linestyle='--',color='gray')
plt.savefig(paths['output']+date+'/DiffPhot_both.png',dpi=300)


plt.figure(6)
plt.clf()
plt.plot(np.append(airmass[0:10],airmass[-50:]),np.append(relflux[0:10],relflux[-50:]),'o')

x = np.append(airmass[0:10],airmass[-50:])
y = np.append(relflux[0:10],relflux[-50:])

opt = np.polyfit(x,y,1)
fit = np.polyval(opt,airmass)
plt.plot(airmass,fit)
plt.xlabel('Airmass')
plt.ylabel('Relative Flux')
plt.title('Airmass Fit for '+targname+': '+date)
plt.savefig(paths['output']+date+'/AirmassFit_'+band+'.png',dpi=300)


relflux_corr = relflux/np.polyval(opt,airmass)


plt.figure(7)
plt.clf()
plt.plot(reltime,relflux_corr,'o',label='AstroPy')
#plt.plot(reltime_g,relflux_g_corr,'o')
plt.xlabel('Time from Mid-Transit (hrs)')
plt.ylabel('Relative Flux')
plt.title('Transit Photometry of '+targname+' in '+band+'$^\prime$: '+date)
plt.axhline(y=1.0,linestyle='--',color='gray')
plt.axvline(x=0.0,linestyle='--',color='gray')
plt.ylim(0.85,1.1)
data = pd.read_csv(paths['output']+'Measurements.tbl.csv')
bjd_aij = data['BJD_TDB'].values
flux_aij = data['rel_flux_T1'].values
relflux_aij = flux_aij/np.median(flux_aij[-50:])
plt.plot((bjd_aij-midtrans)*24.0,relflux_aij,'o',label='AIJ')
plt.legend(loc='lower right',numpoints=1)
plt.savefig(paths['output']+date+'/DiffPhot_r_corr_compare.png',dpi=300)


#aij = pd.read_csv(paths['output']+date+'/XO-2bg.csv',header=1,names=['bjd','flux'])


plt.figure(8)
plt.clf()
aij_norm = np.median(aij['flux'][0:20])
plt.plot((aij['bjd']-midtrans)*24,aij['flux']/aij_norm,'o',label='AIJ')
plt.plot(reltime_g,relflux_g_corr,'o',label='phot_pipe')
plt.xlabel('Time from Mid-Transit (hrs)')
plt.ylabel('Relative Flux')
plt.title('Transit Photometry of '+targname+' in $g^\prime$: '+date)
plt.axhline(y=1.0,linestyle='--',color='gray')
plt.axvline(x=0.0,linestyle='--',color='gray')
plt.legend(numpoints=1,loc='lower left')
plt.savefig(paths['output']+date+'/ComparePhot_g.png',dpi=300)