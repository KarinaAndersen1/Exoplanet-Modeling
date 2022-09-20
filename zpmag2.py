from recursive_glob import recursive_glob
import pickle
import os
import matplotlib.pyplot as plt
import numpy as np

ifiles = recursive_glob('/home/student/Cal/KIC8462852/', 'phot*_i.pck')
imag = []
flux = []
airmass = []
refs = ['KIC8462852', 'ref1', 'ref2', 'ref3', 'ref4', 'ref5', 'ref6', 'ref7', 'ref8', 'ref9']

for file in ifiles:
    dir_path = os.path.dirname(os.path.realpath(file))
    os.chdir(dir_path)
    try:
        data = pickle.load(open(file, 'rb'))
        for ref in refs:
            airmass.append(data[ref]['airmass'])
            flux.append(data[ref]['flux'])
            for i in range(len(data[ref]['flux'])):
                imag.append(data[ref]['imag'])
        print "file %s complete" % file
    except ValueError:
        print "Pickle being a jerk for %s" % file

fluxlist = []
airmasslist = []

for f in flux:
    for i in f:
        fluxlist.append(i)

for a in airmass:
    for i in a:
        airmasslist.append(i)

fluxlist = np.array(fluxlist)
airmasslist = np.array(airmasslist)
imag = np.array(imag)

IM = (imag + 2.5 * np.log(fluxlist))

plt.ion()
plt.figure()
plt.scatter(airmasslist,IM)