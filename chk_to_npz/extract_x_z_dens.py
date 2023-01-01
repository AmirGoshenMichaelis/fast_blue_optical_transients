#! /usr/bin/python3
import os
import sys
import yt
import numpy as np

dir = ''
if len(dir)>0: os.chdir(dir)

for no in range(1110):
    fn = os.path.join(dir,'PN1_hdf5_chk_{:04}'.format(no))
    if not os.path.exists(fn):
        # print('skipping {}'.format(fn))
        continue
    ds = yt.load(fn)
    sl = ds.r[::1024j,0,::1024j]

    ctime = ds.current_time.value
    xx = sl['flash','x'].value
    zz = sl['flash','z'].value
    dens = sl['flash','dens'].value
    velx = sl['flash','velx'].value
    vely = sl['flash','vely'].value
    velz = sl['flash','velz'].value
    temp = sl['flash','temp'].value

    np.savez('xx_zz_dens_time_{:04}.npz'.format(no), X=xx, Z=zz, dens=dens, time=ctime, velx=velx, vely=vely, velz=velz, temp=temp)
