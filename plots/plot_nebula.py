#! /usr/bin/python3

import os
import sys
import math
import numpy as np
from matplotlib import pyplot as pl
from astropy import constants as C
from astropy import units as U

AU   = C.au.cgs.value
Rsun = C.R_sun.cgs.value
Msun = C.M_sun.cgs.value
day  = 24*3600
###
# cmap_array = ['Spectral', 'nipy_spectral', 'inferno', 'coolwarm', 'spring jet', 'gist_earth', 'summer', 'hot', 'Greys', 'rainbow']
cmap = pl.cm.get_cmap('coolwarm', 512)
###

def plot_dens(fn, odir, fontsize=16):
    zlimit = [-15.8, -10.8]
    if not os.path.exists(fn):
        raise Exception(f'File not found (plot_dens) {fn}')

    fig, ax = pl.subplots(1, 1,figsize=(8,8), dpi=150)
    # fig.subplots_adjust(wspace=0.01, hspace=0.05)

    ldata = np.load(fn)
    xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    # xlimit  = [-33.35, 33.35]
    # ylimit  = [-33.35, 33.35]
    X       = ldata['X']/AU
    Y       = ldata['Z']/AU
    Z       = np.log10(ldata['dens'])
    ctime   = int(ldata['time']/day)

    decimate = 50
    cf = ax.contourf(X, Y, Z, np.linspace(zlimit[0], zlimit[1], 512), cmap=cmap, extend='both')
    # ax.set_xlim(xlimit[0],xlimit[1])
    # ax.set_ylim(ylimit[0],ylimit[1])
    # ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.067), ylimit[1]*(1-0.2), f't={ctime} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    ax.tick_params(axis='both', labelsize=fontsize)
    ax.minorticks_on()
    ax.set_xlabel(r'$X/AU$', size=fontsize)
    ax.set_ylabel(r'$Z/AU$', size=fontsize)

    fig.text(0.4, 0.886, r'$\log~\rho~[\rm{g~cm^{-3}}]$', fontsize=fontsize)
    cax = pl.axes([0.91, 0.12, 0.022, 0.74])
    cbobj = fig.colorbar(cf, cax=cax, format='%2.1f', ticks=np.round(np.arange(zlimit[0],zlimit[1]+0.1,0.5),1))
    cbobj.ax.tick_params(axis='both', labelsize=fontsize)

    ofn = os.path.join(odir,f'nova_dens_{ctime}.png')
    fig.savefig(ofn, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)
    pl.close()
###

def plot_vtot(fn, odir, fontsize=16):
    zlimit = [4.6, 7.4]

    fig, ax = pl.subplots(1, 1,figsize=(8,8), dpi=150)
    # fig.subplots_adjust(wspace=0.01, hspace=0.05)

    if not os.path.exists(fn):
        raise Exception(f'File not found {fn}')
    
    ldata = np.load(fn)
    xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    # xlimit  = [-33.35, 33.35]
    # ylimit  = [-33.35, 33.35]
    X       = ldata['X']/AU
    Y       = ldata['Z']/AU
    Vx      = ldata['velx']
    Vy      = ldata['velz']
    Vabs    = np.sqrt(Vx**2+Vy**2)
    Z       = np.log10( np.sqrt( ldata['velx']**2+ldata['vely']**2+ldata['velz']**2 ) )
    ctime   = int(ldata['time']/day)

    decimate = 50
    cf = ax.contourf(X, Y, Z, np.linspace(zlimit[0], zlimit[1], 512), cmap=cmap, extend='both')
    Q  = ax.quiver(X[::decimate,::decimate],Y[::decimate,::decimate],Vx[::decimate,::decimate]/Vabs[::decimate,::decimate],Vy[::decimate,::decimate]/Vabs[::decimate,::decimate], units='xy',color=[0.1,0.1,0.1], edgecolors=[0.0,0.0,0.0])

    # ax.set_xlim(xlimit[0],xlimit[1])
    # ax.set_ylim(ylimit[0],ylimit[1])
    ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.067), ylimit[1]*(1-0.2), f't={ctime} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    ax.tick_params(axis='both', labelsize=fontsize)
    ax.minorticks_on()
    ax.set_xlabel(r'$X/AU$', size=fontsize)
    ax.set_ylabel(r'$Z/AU$', size=fontsize)

    fig.text(0.4, 0.886, r'$\log ~ \left| \vec V \right| ~ \rm{[cm/s]}$', fontsize=fontsize)
    cax = pl.axes([0.91, 0.12, 0.022, 0.74])
    cbobj = fig.colorbar(cf, cax=cax, format='%2.1f', ticks=np.round(np.arange(zlimit[0],zlimit[1]-0.1,0.1),1))
    cbobj.ax.tick_params(axis='both', labelsize=fontsize)

    ofn = os.path.join(odir,f'nova_vtot_{ctime}.png')
    fig.savefig(ofn, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)
    pl.close()
###

def plot_temp(fn, odir, fontsize=16):
    zlimit = [3., 3.6]

    if not os.path.exists(fn):
        raise Exception(f'File not found {fn}')
    
    fig, ax = pl.subplots(1, 1,figsize=(8,8), dpi=150)
    # fig.subplots_adjust(wspace=0.01, hspace=0.05)

    ldata = np.load(fn)
    xlimit  = np.array([ldata['X'].min(), ldata['X'].max()])/AU
    ylimit  = np.array([ldata['Z'].min(), ldata['Z'].max()])/AU
    # xlimit  = [-33.35, 33.35]
    # ylimit  = [-33.35, 33.35]
    X       = ldata['X']/AU
    Y       = ldata['Z']/AU
    Z       = np.log10(ldata['temp'])
    ctime   = int(ldata['time']/day)

    decimate = 50
    cf = ax.contourf(X, Y, Z, np.linspace(zlimit[0], zlimit[1], 512), cmap=cmap, extend='both')
    # ax.set_xlim(xlimit[0],xlimit[1])
    # ax.set_ylim(ylimit[0],ylimit[1])
    ax.set_aspect('equal')
    ax.text(xlimit[0]*(1-0.067), ylimit[1]*(1-0.2), f't={ctime} d', bbox=dict(boxstyle="round", fc="w"), fontsize=fontsize)

    ax.tick_params(axis='both', labelsize=fontsize)
    ax.minorticks_on()
    ax.set_xlabel(r'$X/AU$', size=fontsize)
    ax.set_ylabel(r'$Z/AU$', size=fontsize)

    fig.text(0.4, 0.886, r'$\log ~ T ~ \rm{[K]}$', fontsize=fontsize)
    cax = pl.axes([0.91, 0.12, 0.022, 0.74])
    cbobj = fig.colorbar(cf, cax=cax, format='%2.1f', ticks=np.round(np.arange(zlimit[0],zlimit[1]-0.1,0.1),1))
    cbobj.ax.tick_params(axis='both', labelsize=fontsize)

    ofn = os.path.join(odir,f'nova_temp_{ctime}.png')
    fig.savefig(ofn, dpi=150, facecolor='w', edgecolor='w', orientation='portrait', transparent=False, bbox_inches=None, pad_inches=0.1)
    pl.close()
###

def main():
    chk_dir = '/home/amirm/code/fast_blue_optical_transients/dd/t5'
    odir    = '/home/amirm/code/fast_blue_optical_transients/plots/'
    fname = [os.path.join(chk_dir,f'xx_zz_dens_time_{i:04}.npz') for i in (160, 160, 160) ]
    plot_dens(fname[0],odir)
    plot_temp(fname[1],odir)
    plot_vtot(fname[2],odir)
###

if __name__ == "__main__":
    main()
    print(f'done {sys.argv[0]}')
###
