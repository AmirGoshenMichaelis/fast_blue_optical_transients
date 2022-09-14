#! /usr/bin/python3
###
import os
import numpy as np
from matplotlib import pyplot as pl
from multiprocessing import Pool
Rsun = 6.957e10
AU   = 1.49597871e+13
day  = 24*3600
###
def remove_nan(x):
    return x[~np.isnan(x)]
###
def Rot(theta, phi):
    theta = theta*np.pi/180.
    phi = phi*np.pi/180.
    R = np.array(
        [[np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
         [-np.sin(phi)             , np.cos(phi)              , 0],
         [np.cos(phi)*np.sin(theta), np.sin(theta)*np.sin(phi), np.cos(theta)]]
    )
    return R
###
dir = '/home/amirm/code/fast_blue_optical_transients/dd/exp3/'
odir = '/home/amirm/code/fast_blue_optical_transients/photosphere/exp3/'
os.chdir(dir)
cmap = pl.cm.get_cmap('Spectral', 512)
###
def calc(no):
    fn = os.path.join(dir,'xx_zz_dens_time_{:04}.npz'.format(no))
    if not os.path.exists(fn):
        print('skipping {}'.format(fn))
        return
    dd = np.load(fn)
    X = dd['X']/AU
    Y = dd['Z']/AU
    Z = np.log10(dd['dens'])
    ctime = dd['time']/day
    # for fix_kappa in [0.06, 0.3, 1.5, 0.1, 0.2]:
    for fix_kappa in [0.1,]:
        fn = os.path.join(dir, 'TRL_eff_{:04}_fix_kappa_{}.npz'.format(no, fix_kappa) )
        if not os.path.exists(fn):
            print('skipping {}'.format(fn))
            continue
        dd = np.load(fn)

        t = dd['time']
        T = dd['Teff']
        R = dd['Reff']
        L = dd['Leff']
        XYZ   = dd['XYZeff']
        p1    = dd['p1']
        theta = dd['theta']
        phi   = dd['phi']

        Temp =0.0
        pr_phi = 0.
        # for pr_theta in [10, 30, 50, 80,]:
        for pr_theta in [10, 50, 80,]:
            fig, ax = pl.subplots()
            cf    = ax.contourf(X, Y, Z, np.linspace(-18., -8, 512), cmap=cmap, extend='both')
            fig.text(0.45, 0.9, '$\\rho~[\\rm{g/cm^3}]$', fontsize=14)
            cax = pl.axes([0.71, 0.12, 0.022, 0.74])
            cbobj = fig.colorbar(cf, cax=cax, format='%2.f', ticks=np.arange(-18., -7., 1))
            # cax.tick_params(axis='both', which='major', labelsize=14)

            j = np.argmin(np.abs(theta-pr_theta))
            k = np.argmin(np.abs(phi-pr_phi))
            pp = [~np.isnan(XYZ[j,k,:,l]) for l in range(3)]
            pp = pp[0]&pp[1]&pp[2]

            teff_pp = (~np.isnan(T[j,k,:])) & (T[j,k,:]>Temp)
            pp = pp[0]&pp[1]&pp[2]&teff_pp
            vec = XYZ[j,k,pp,:].transpose()
            R   = Rot(theta[j],phi[k])
            los_sys = np.dot(R,vec)
            pos = los_sys[2,:]>=0

            ax.plot(XYZ[j,k,pp,0][pos]/AU,XYZ[j,k,pp,2][pos]/AU, 'D', markerfacecolor='Cyan', markeredgecolor='Red', markersize=5)

            ax.set_aspect('equal')
            ax.text(ax.get_xlim()[0]*(1-0.067), ax.get_ylim()[1]*(1-0.15), f't={ctime:.0f} d', bbox=dict(boxstyle="round", fc="w"), fontsize=14)
            ax.set_xlabel('X [AU]', fontsize=14)
            ax.set_ylabel('Z [AU]', fontsize=14)

            ofn = 'kappa_{}/{}/kappa_{}_dens_photo_{}_PN1_{:04}.png'.format(fix_kappa, pr_theta, fix_kappa, pr_theta, no)
            fig.savefig(os.path.join(odir, ofn))
            pl.close(fig)

arg_list = list()
for no in range(1,300):
    arg_list.append( (no,) )

with Pool(8) as p:
    data = p.starmap(calc, arg_list)

