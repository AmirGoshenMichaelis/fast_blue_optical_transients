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
chk_dir = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/npz_9/dd/'
odir    = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/npz_9/plots/'
fname = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (1, 45, 89) ]
cmap = pl.cm.get_cmap('coolwarm', 512)
###
def plot_dens_ax(ax, no, fontsize=16):

    if not os.path.exists(fname[no]):
        raise Exception(f'File not found {fname[no]}')
    
    ldata = np.load(fname[no])
    # xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    # ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    xlimit  = [-80.,80.]
    ylimit  = [-80.,80.]
    X       = ldata['X']/AU
    Y       = ldata['Z']/AU
    Z       = np.log10(ldata['dens'])
    Vx      = ldata['velx']
    Vy      = ldata['velz']
    ctime   = ldata['time']/day

    decimate = 50
    cf = ax.contourf(X, Y, Z, np.linspace(-16., -8, 512), cmap=cmap, extend='both')
    if no % 2 :
        ax.set_xlim(xlimit[0],xlimit[1]-1)
    else:
        ax.set_xlim(xlimit[0],xlimit[1])
    ax.set_ylim(ylimit[0],ylimit[1])
    ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.067), ylimit[1]*(1-0.2), f't={int(ctime)} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    return cf

def plot_vtot_ax(ax, no, fontsize=16):

    if not os.path.exists(fname[no]):
        raise Exception(f'File not found {fname[no]}')
    
    ldata = np.load(fname[no])
    # xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    # ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    xlimit  = [-80.,80.]
    ylimit  = [-80.,80.]
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
    cf = ax.contourf(X, Y, Z, np.linspace(6.5, 9.8, 512), cmap=cmap, extend='both')
    Q  = ax.quiver(X[::decimate,::decimate],Y[::decimate,::decimate],Vx[::decimate,::decimate]/Vabs[::decimate,::decimate],Vy[::decimate,::decimate]/Vabs[::decimate,::decimate], units='xy',color=[0.1,0.1,0.1], edgecolors=[0.0,0.0,0.0]) #scale=scale[no], 
    if no % 2 :
        ax.set_xlim(xlimit[0],xlimit[1]-1)
    else:
        ax.set_xlim(xlimit[0],xlimit[1])
    ax.set_ylim(ylimit[0],ylimit[1])
    ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.067), ylimit[1]*(1-0.2), f't={int(ctime)} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    return cf

def dens_(fontsize=18):
    xy = 'xz'
    nx=1
    ny=3
    figsize=(nx*6,ny*6)
    fig, axs = pl.subplots(ny,nx,figsize=figsize, dpi=150 , sharex=True, sharey=True)
    if type(axs) is not np.ndarray: # not isinstance(bb,np.ndarray)
        axs = np.array([[axs,],])
    elif len(axs.shape)==1:
        axs.shape=(axs.shape[0],1)

    fig.subplots_adjust(wspace=0.01, hspace=0.13)
    for i in range(ny):
        for j in range(nx):
            no_ij = i*nx+j
            print(f'Processing subplot {i} {j} {no_ij}')

            cf = plot_dens_ax(axs[i,j], no_ij, fontsize=fontsize)
            # if no_ij==0:
            #     axs[i,j].set_title(r'$\log~\rho~\rm{[g~cm^{-3}]}$', fontsize=fontsize)
            #     cax = pl.axes([0.89, 0.66, 0.022, 0.2])
            # if no_ij==1:
            #     cax = pl.axes([0.89, 0.39, 0.022, 0.2])
            # if no_ij==2:
            #     cax = pl.axes([0.89, 0.13, 0.022, 0.2])
            if no_ij==0:
                cax = pl.axes([0.89, 0.13, 0.022, 0.73])
                cbobj = fig.colorbar(cf, cax=cax, format='%2.f', ticks=np.arange(-16., -7., 1))
                cbobj.ax.tick_params(axis='both', labelsize=fontsize)

            if j==0:
                axs[i,j].set_ylabel(r'${}/AU$'.format(xy[1].upper()), size=fontsize)
            if i==ny-1:
                xlimit = axs[i,j].get_xlim()
                axs[i,j].set_xticks(np.arange(xlimit[0]+20,xlimit[1]+1,20))
                axs[i,j].set_xlabel(r'${}/AU$'.format(xy[0].upper()), size=fontsize)
            axs[i,j].tick_params(axis='both', labelsize=fontsize)
            axs[i,j].minorticks_on()

        axs[i,0].set_ylabel(r'${}/AU$'.format(xy[1].upper()), size=fontsize)

    ofn = os.path.join(odir,f'dens_{xy}_{ny}_{nx}.png')
    fig.savefig(ofn, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)
    pl.close()
###
def main():
    dens_()
###
if __name__ == "__main__":
    main()
###
