#! /usr/bin/python3

import os
import sys
import math
import numpy as np
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
chk_dir = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/npz_2/dd/'
odir    = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/npz_2/'
ffn = [None, None]
ffn[0] = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (7,11,15,19,22,26) ]
ffn[1] = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (170,175,180,185,190,195) ]
scl = [None,None]
scl[0] = np.array([0.05]*6)
scl[1] = np.array([0.05]*6)
vvelscl = [None, None]
vvelscl[0] = [5000.e5, 5000.e5, 10000.e5, 10000.e5, 5000.e5, 5000.e5]
vvelscl[1] = [5000.e5, 5000.e5, 10000.e5, 10000.e5, 10000.e5, 10000.e5]

fname_set = 2
fname = ffn[fname_set-1]
scale = scl[fname_set-1]
vel_scale = vvelscl[fname_set-1]

cmap = pl.cm.get_cmap('Spectral', 512)
###
def plot_dens_ax(ax, no, fontsize=16):

    if not os.path.exists(fname[no]):
        raise Exception(f'File not found {fname[no]}')
    
    ldata = np.load(fname[no])
    xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    X       = ldata['X']/AU
    Y       = ldata['Z']/AU
    Vx      = ldata['velx']
    Vy      = ldata['velz']
    Vabs    = np.sqrt(Vx**2+Vy**2)
    # Z       = np.log10(Vabs)
    Z       = np.log10( np.sqrt( ldata['velx']**2+ldata['vely']**2+ldata['velz']**2 )   )
    ctime   = ldata['time']/day

    decimate = 50
    # print(Z.min(),Z.max())
    cf = ax.contourf(X, Y, Z, np.linspace(2.9, 9.7, 512), cmap=cmap, extend='both')
    Q  = ax.quiver(X[::decimate,::decimate],Y[::decimate,::decimate],Vx[::decimate,::decimate]/Vabs[::decimate,::decimate],Vy[::decimate,::decimate]/Vabs[::decimate,::decimate], scale=scale[no], units='xy',color=[0.1,0.1,0.1], edgecolors=[0.0,0.0,0.0])
    if no % 2 :
        ax.set_xlim(xlimit[0],xlimit[1]-1)
    else:
        ax.set_xlim(xlimit[0],xlimit[1])
    ax.set_ylim(ylimit[0],ylimit[1])
    ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.067), ylimit[1]*(1-0.2), f't={int(ctime)} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    return cf,Q

def dens_(fontsize=16):
    xy = 'xz'
    nx=2
    ny=3
    figsize=(nx*4,ny*4)
    fig, axs = pl.subplots(ny,nx,figsize=figsize, dpi=150 , sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.01, hspace=0.05)
    for i in range(ny):
        for j in range(nx):
            no = i*nx+j
            print(f'Processing subplot {i} {j} {no}')
            cf, Q = plot_dens_ax(axs[i,j], no, fontsize=fontsize)
            if i==ny-1:
                axs[i,j].set_xlabel(r'${}/AU$'.format(xy[0].upper()), size=fontsize)

            # axs[i,j].quiverkey(Q, 0.74, 0.95, vel_scale[no], f'{(vel_scale[no]/1e5):.0f} [km/s]', labelpos='S', fontproperties={'size': fontsize})
        axs[i,0].set_ylabel(r'${}/AU$'.format(xy[1].upper()), size=fontsize)

    fig.text(0.4, 0.886, r'$\left| \vec V \right| ~ \rm{[cm/s]}$', fontsize=fontsize)
    cax = pl.axes([0.91, 0.12, 0.022, 0.74])
    cbobj = fig.colorbar(cf, cax=cax, format='%2.1f', ticks=np.arange(2.9, 9.8, 1))

    ofn = os.path.join(odir,f'vtot_vec_{xy}_{ny}_{nx}_{fname_set}.png')
    fig.savefig(ofn, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)
    pl.close()
###
def main():
    dens_()
###
if __name__ == "__main__":
    main()
###
