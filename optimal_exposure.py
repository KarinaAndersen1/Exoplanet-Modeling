'''
Get all KIC846 files in each band
Find all PanSTARRS refs
Loop through them
  Optimal aperture on each ref star
  Save: total counts, max counts, airmass, fwhm, integration time
  Graph: IM vs airmass, Peak counts as a function of mag and seeing at certain airmass, 2 std deviations below saturation
  Fit a 5 dimensional surface in space to data
    #add if it makes you happy, gonna soak up the sun to spotify

'''
# Imports and variables
import numpy as np 
import matplotlib.pyplot as plt
import phot_pipe as pp
from quick_image import *

targname= 'KIC8462852'
obstype= 'Tabby'
rastr = '20:06:15.45'
decstr = '+44:27:24.79'
bands = ['g', 'r', 'i', 'z', 'V']

# Get paths to Tabby data and all dates of observation
paths = pp.get_paths(targname=targname,obstype=obstype)
dates = pp.get_dates(targname=targname)

# Find Tabby's panstarrs refs
refs = pp.get_panstarrs_refs(ra=rastr,dec=decstr,maxMag=15)
x = []
y = []
for ra in len(refs['RAdeg']):
    x.append(refs['RAdeg'])
    y.append(refs['DECdeg'])

mag = []
maxcts = []
airmass = []
fwhm = []
inttime = []

for date in dates:
    for band in bands:
        files, ct = pp.get_files(prefix=targname, date=date, tag=band)
        maglist = band + 'band'
        for file in files:
            im,h = read_image(file,plot=False)
            for i in len(refs['RAdeg']):
                mag.append(refs[maglist])
                data = pp.optimal_aperture(im, h, x=x[i], y=y[i])
                maxcts.append(data['totcounts'])
                airmass.append(data['secz'])
                fwhm.append(data['fwhm'])
                inttime.append(data['exptime'])





               

            
            

def optimal_exposure(magnitude=none,band=none,seeing=none,airmass=none)