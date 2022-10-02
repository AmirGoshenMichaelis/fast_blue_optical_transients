#! /usr/bin/python3
###
import os
import sys
import json
# import numpy as np
# from matplotlib import pyplot as pl

os.chdir('/home/amirm/code/fast_blue_optical_transients/Perley2018')
dd = {'mjd':[], 'inst':[], 'filt':[], 'mag':[], 'ABmag':[], 'err':[], 'source':[]}
with open('AT2018bow_photometry_table.1808.00969.dat', 'r') as fid:
    for line in fid:
        if line[0]=='#':
            print('comment:',line)
            continue
        ll = line.split()
        if len(ll)<6:
            print('skipping data point:',line)
            continue
        try: #mjd       inst.   filt.  mag    ABmag  err  source
            mjd   = float(ll[0])
            inst  = ll[1]
            filt  = ll[2]
            mag   = float(ll[3])
            ABmag = float(ll[4])
            err   = float(ll[5])
        except ValueError:
            print('skipping data point:',line)
            continue
        dd['mjd'].append(mjd)
        dd['inst'].append(inst)
        dd['filt'].append(filt)
        dd['mag'].append(mag)
        dd['ABmag'].append(ABmag)
        dd['err'].append(err)
        if len(ll)>=7:
            dd['source'].append(ll[6])
        else:
            dd['source'].append(None)
        
with open('AT2018bow_photometry_table.1808.00969.json', 'w') as fid:
    json.dump(dd, fid)
###
# with open('AT2018bow_photometry_table.1808.00969.json', 'r') as fid:
#     dd = json.load(fid)

