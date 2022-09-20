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
            rmag.append(data[ref]['rmag'])
        print "file %s complete" % file
    except ValueError:
        print "Pickle being dumb for %s" % file

lengthf = []
sumf = []
for f in range(len(flux)):
    lengthf.append(len(flux[f]))
    sumf.append(sum(flux[f]))

lengthf = np.array(lengthf)
sumf = np.array(sumf)
fluxavg = (sumf / lengthf)

lengtha = []
suma = []
for a in range(len(airmass)):
    lengtha.append(len(airmass[a]))
    suma.append(sum(airmass[a]))

lengtha = np.array(lengtha)
suma = np.array(suma)
airmassavg = (lengtha/suma)

imag = np.array(imag)
IM = (imag + 2.5 * np.log(fluxavg))

plt.ion()
plt.figure()
plt.scatter(airmassavg,IM)


