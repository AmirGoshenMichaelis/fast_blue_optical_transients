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

AU   = C.au.cgs.value
Rsun = C.R_sun.cgs.value
Msun = C.M_sun.cgs.value
day = 24*3600
###
# @derived_field(name=('gas', 'vtot'), units="cm/s", sampling_type="cell", force_override=True)
# def _vtot(field, data):
#     # import pdb; pdb.set_trace()
#     # data._debug()
#     vx = data[('flash','velx')].to('cm/s')
#     vy = data[('flash','vely')].to('cm/s')
#     vz = data[('flash','velz')].to('cm/s')
#     return np.sqrt(vx*vx+vy*vy+vz*vz)
###
def already_exist(fn, odir, no):
    for src in ['xz']:
        ofn = os.path.join(odir,f'temp_{src}_{no:04}.png')
        if not os.path.exists(ofn):
            return False
    return True
###
def plot_temp_ax(param, ax, no=1, fontsize=14):
    ldata = np.load(param['fn'])
    xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    X       = ldata['X']/AU
    Y       = ldata['Z']/AU
    Z       = np.log10(ldata['temp'])
    Vx      = ldata['velx']
    Vy      = ldata['velz']
    Vz      = ldata['velz']
    Vtot    = np.sqrt( Vx**2 + Vy**2 + Vz**2 )
    Vx      = Vx/Vtot
    Vy      = Vy/Vtot
    ctime   = ldata['time']/day


    # cmap_array = ['Spectral', 'nipy_spectral', 'inferno', 'coolwarm', 'spring jet', 'gist_earth', 'summer', 'hot', 'Greys', 'rainbow']
    cmap = pl.cm.get_cmap('coolwarm', 512)
    decimate = 20
    cf = ax.contourf(X, Y, Z, np.linspace(3.5, 4.3, 512), cmap=cmap, extend='both')
    Q  = ax.quiver(X[::decimate,::decimate],Y[::decimate,::decimate],Vx[::decimate,::decimate],Vy[::decimate,::decimate], scale=0.1, units='xy',color=[0.1,0.1,0.1], edgecolors=[0.0,0.0,0.0])
    ax.set_xlim(xlimit[0],xlimit[1])
    ax.set_ylim(ylimit[0],ylimit[1])
    ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.03), ylimit[1]*(1-0.1), f't={ctime:.0f} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    return cf,Q
###
def temp_(param, plane='xy', fontsize=16):
    nx=1
    ny=1
    xy = plane
    figsize=(8,8)
    fig, axs = pl.subplots(ny,nx,figsize=figsize, dpi=150 , sharex=True, sharey=True)
    if type(axs) is not np.ndarray: # not isinstance(bb,np.ndarray)
        axs = np.array([[axs,],])
    elif len(axs.shape)==1:
        axs.shape=(axs.shape[0],1)
    fig.subplots_adjust(wspace=0.01, hspace=0.05)
    for i in range(ny):
        for j in range(nx):
            print(f'Processing subplot {i} {j} {i*nx+j}')
            cf, Q = plot_temp_ax(param, ax=axs[i,j], no=i*nx+j, fontsize=fontsize)
            if i==ny-1:
                axs[i,j].set_xlabel(r'${}/R_\odot$'.format(xy[0].upper()), size=fontsize)
            # if i==0 and j==nx-1:
            #     yloc = 1.09
            #     vel_scale = 1000.e5
            #     axs[i,j].quiverkey(Q, 0.74, yloc, vel_scale, f'{(vel_scale/1e5):.0f} [km/s]', labelpos='S', fontproperties={'size': fontsize})
        axs[i,0].set_ylabel(r'${}/R_\odot$'.format(xy[1].upper()), size=fontsize)

    fig.text(0.45, 0.9, r'$T~[K]$', fontsize=fontsize)
    cax = pl.axes([0.91, 0.12, 0.022, 0.74])
    cbobj = fig.colorbar(cf, cax=cax, format='%2.1f', ticks=np.arange(3.5,4.31,0.2))

    ofn = os.path.join(param['chk']['odir'],f'temp_{plane}_{param["no."]:04}.png')
    fig.savefig(ofn, dpi=300, facecolor='w', edgecolor='w', orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)
    pl.close()
###
def calc(param):
    temp_(param, plane='xz')
###
def main():
    param = dict()
    param['chk'] = dict()
    param['chk']['dir'] = '/home/amirm/code/fast_blue_optical_transients/dd/exp1/'
    param['chk']['odir'] = '/home/amirm/code/fast_blue_optical_transients/pp_temp/exp1/'
    param['chk']['range'] = [1, 655]
    param['chk']['fn'] = 'xx_zz_dens_time_{:04}.npz'
    param['no.'] = 185
    param['fn'] = os.path.join(param['chk']['dir'],param['chk']['fn'].format(param['no.']))

    os.chdir(param['chk']['dir'])
    calc(param)
###
if __name__ == "__main__":
    main()
    print(f'done {sys.argv[0]}')
###
