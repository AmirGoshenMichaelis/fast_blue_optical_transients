#! /usr/bin/python3
###
import enum
import os
import sys
import json
import numpy as np
from matplotlib import pyplot as pl

os.chdir('/home/amirm/code/fast_blue_optical_transients/Perley2018')
with open('AT2018bow_photometry_table.1808.00969.json', 'r') as fid:
    dd = json.load(fid)


mjd = np.array([])
mag = np.array([])
err = np.array([])

for i, filt in enumerate(dd['filt']):
    if 'V' == filt.upper():
        mjd = np.append(mjd, dd['mjd'][i])
        mag = np.append(mag, dd['mag'][i])
        err = np.append(mag, dd['err'][i])
mag = mag - 5.*np.log10(60e6/10)
pl.plot(mjd,mag, 'rx')
pl.gca().invert_yaxis()
pl.show()

np.savez('AT2018bow.npz', mjd=mjd, mag=mag, err=err)
