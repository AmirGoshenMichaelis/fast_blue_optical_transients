#! /usr/bin/python3

import os
import sys
import math
from xml.dom.pulldom import ErrorHandler
import numpy as np
from sklearn.model_selection import LeaveOneGroupOut
import yt
from yt import derived_field
from matplotlib import pyplot as pl
from astropy import constants as C
from astropy import units as U
import h5py

AU   = C.au.cgs.value
Rsun = C.R_sun.cgs.value
Msun = C.M_sun.cgs.value
day  = 24*3600
###
# cmap_array = ['Spectral', 'nipy_spectral', 'inferno', 'coolwarm', 'spring jet', 'gist_earth', 'summer', 'hot', 'Greys', 'rainbow']
chk_dir = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/npz_9/dd/'
odir    = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/npz_9/plots/'
ffn = [None,]*3
ffn[0] = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (7,11,15,19,22,26) ]
ffn[1] = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (44,60,63,66,68,71) ]
# ffn[2] = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (175,180,185,190) ]
ffn[2] = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (95,103,111,120) ]
fname_set = 3
fname = ffn[fname_set-1]
cmap = pl.cm.get_cmap('coolwarm', 512)
###
def plot_dens_ax(ax, no, fontsize=16):

    if not os.path.exists(fname[no]):
        raise Exception(f'File not found {fname[no]}')
    
    ldata = np.load(fname[no])
    xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    X       = ldata['X']/AU
    Y       = ldata['Z']/AU
    Z       = np.log10(ldata['dens'])
    Vx      = ldata['velx']
    Vy      = ldata['velz']
    ctime   = ldata['time']/day

    decimate = 50
    cf = ax.contourf(X, Y, Z, np.linspace(-17., -9, 512), cmap=cmap, extend='both')
    if no % 2 :
        ax.set_xlim(xlimit[0],xlimit[1]-1)
    else:
        ax.set_xlim(xlimit[0],xlimit[1])
    ax.set_ylim(ylimit[0],ylimit[1])
    ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.067), ylimit[1]*(1-0.2), f't={int(ctime)} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    return cf

def dens_(fontsize=16):
    xy = 'xz'
    nx=2
    ny=2
    figsize=(nx*4,ny*4)
    fig, axs = pl.subplots(ny,nx,figsize=figsize, dpi=150 , sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.01, hspace=0.05)
    for i in range(ny):
        for j in range(nx):
            print(f'Processing subplot {i} {j} {i*nx+j}')
            cf = plot_dens_ax(axs[i,j], i*nx+j, fontsize=fontsize)
            axs[i,j].tick_params(axis='both', labelsize=fontsize)
            axs[i,j].minorticks_on()
            if i==ny-1:
                axs[i,j].set_xlabel(r'${}/AU$'.format(xy[0].upper()), size=fontsize)
        axs[i,0].set_ylabel(r'${}/AU$'.format(xy[1].upper()), size=fontsize)

    fig.text(0.4, 0.886, r'$\log~\rho~[g~cm^{-3}]$', fontsize=fontsize)
    cax = pl.axes([0.91, 0.12, 0.022, 0.74])
    cbobj = fig.colorbar(cf, cax=cax, format='%2.f', ticks=np.arange(-17., -8., 1))
    cbobj.ax.tick_params(axis='both', labelsize=fontsize)

    ofn = os.path.join(odir,f'dens_{xy}_{ny}_{nx}_{fname_set}.png')
    fig.savefig(ofn, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)
    pl.close()
###
def main():
    dens_()
###
if __name__ == "__main__":
    main()
###
