#! /usr/bin/python3
###
from cProfile import label
import os
import sys
import numpy as np
from matplotlib import pyplot as pl
from multiprocessing import Pool
Rsun = 6.957e10
AU   = 1.49597871e+13
day  = 24*3600
sigma_sb = 5.6703744191844314e-05
Lsun = 3.828e+33
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
exp_dir = sys.argv[2]
base_dir = '/home/amirm/code/fast_blue_optical_transients/'
dir  = os.path.join(base_dir,'dd',exp_dir)
odir = os.path.join(base_dir,'rph',exp_dir)
at2018cow_dir = os.path.join(base_dir,'Perley2018')
os.chdir(dir)
cmap = pl.cm.get_cmap('Spectral', 512)
Temp = float(sys.argv[1])*1e3
###
def calc(no, fix_kappa):
    fn = os.path.join(dir,'xx_zz_dens_time_{:04}.npz'.format(no))
    if not os.path.exists(fn):
        print('skipping {}'.format(fn))
        return
    dd = np.load(fn)
    X = dd['X']/AU
    Y = dd['Z']/AU
    Z = np.log10(dd['dens'])
    # T = dd['temp']
    ctime = dd['time']/day

    fn = os.path.join(dir, 'TRL_eff_{:04}_fix_kappa_{}.npz'.format(no, fix_kappa) )
    if not os.path.exists(fn):
        print('skipping {}'.format(fn))
        return
    dd = np.load(fn)

    t = dd['time']
    T = dd['Teff']
    R = dd['Reff']
    L = dd['Leff']
    XYZ   = dd['XYZeff']
    p1    = dd['p1']
    theta = dd['theta']
    phi   = dd['phi']

    pr_phi = 0.
    r_ph = []
    for pr_theta in [10, 30, 50, 80,]:
        # fig, ax = pl.subplots()
        # cf    = ax.contourf(X, Y, Z, np.linspace(-18., -8, 512), cmap=cmap, extend='both')
        # fig.text(0.45, 0.9, '$\\rho~[\\rm{g/cm^3}]$', fontsize=14)
        # cax = pl.axes([0.71, 0.12, 0.022, 0.74])
        # cbobj = fig.colorbar(cf, cax=cax, format='%2.f', ticks=np.arange(-18., -7., 1))
        # # cax.tick_params(axis='both', which='major', labelsize=14)

        j = np.argmin(np.abs(theta-pr_theta))
        k = np.argmin(np.abs(phi-pr_phi))
        pp = [~np.isnan(XYZ[j,k,:,l]) for l in range(3)]
        teff_pp = (~np.isnan(T[j,k,:])) & (T[j,k,:]>Temp)
        pp = pp[0]&pp[1]&pp[2]&teff_pp
        vec = XYZ[j,k,pp,:].transpose()
        R   = Rot(theta[j],phi[k])
        los_sys = np.dot(R,vec)
        pos = los_sys[2,:]>=0
        if np.sum(pos)==0:
            photosphere_radios = 1e-16
        else:
            teff_pos = T[j,k,pp][pos]
            # dr = (los_sys[0,:].max()-los_sys[0,:].min())/np.sum(pos)
            dr = np.append(0., np.diff(los_sys[0,pos]))
            dr = np.mean(dr)
            if dr==0: dr=1e-16
            # photosphere_radios = np.sum( np.abs(2.*np.pi*los_sys[0,pos]*dr) )*sigma_sb*(8224**4) 
            photosphere_radios = 4. * np.sum( np.abs(2.*np.pi*los_sys[0,pos]*dr) )*sigma_sb*(3e3**4) 
            # photosphere_radios = np.sqrt( np.sum( np.abs(2.*los_sys[0,pos]*dr) ) )
            # photosphere_radios = np.sqrt( np.sum( np.abs(los_sys[0,pos]*dr)*teff_pos**4/np.sum(teff_pos**4) ) )
            # photosphere_radios = np.max(los_sys[0,pos])
        r_ph.append( photosphere_radios )
        # xx = vec[0,pos]/AU
        # yy = vec[2,pos]/AU
        # ax.plot(xx,yy, 'D', markerfacecolor='Cyan', markeredgecolor='Red', markersize=5)

        # ax.set_aspect('equal')
        # ax.text(ax.get_xlim()[0]*(1-0.067), ax.get_ylim()[1]*(1-0.15), f't={ctime:.0f} d', bbox=dict(boxstyle="round", fc="w"), fontsize=14)
        # ax.set_xlabel('X [AU]', fontsize=14)
        # ax.set_ylabel('Z [AU]', fontsize=14)
        # ofn = 'kappa_{}_dens_photo_{}_PN1_{:04}.png'.format(fix_kappa, pr_theta, no)
        # fig.savefig(os.path.join(odir, ofn))
        # pl.close(fig)
    return ctime, r_ph
###
def do_calc(rph_npz_fn, fix_kappa):
    arg_list = list()
    for no in range(1,500,1):
        arg_list.append( (no, fix_kappa) )

    # data = []
    # for arg in arg_list:
    #     res = calc(*arg)
    #     data.append(res)

    with Pool(10) as p:
        data = p.starmap(calc, arg_list)

    first_nonNone_loc = 0
    for ii, info in enumerate(data): 
        if info is not None:
            first_nonNone_loc = ii
            break
    n_t = len(data)
    n_theta = len(data[first_nonNone_loc][1])
    r_ph = np.zeros((n_t,n_theta))
    t = np.zeros((n_t,1))
    for i,dd in enumerate(data):
        if dd is None:
            t[i] = np.nan
            r_ph[i,:] = np.nan
            continue
        t[i] = dd[0]
        r_ph[i,:] = dd[1]
    np.savez(rph_npz_fn, t=t, r_ph=r_ph)
###
def main():
    # fix_kappa_list = [0.005, 0.01, 0.06, 0.3, 1.5, 0.1, 0.2]
    fix_kappa_list = [0.005, 0.01,]
    for fix_kappa in fix_kappa_list:
        rph_npz_fn = os.path.join(odir, f'mag_{fix_kappa}_T_{Temp/1e3:.0f}e3.npz')
        if not os.path.exists(rph_npz_fn):
            do_calc(rph_npz_fn, fix_kappa)

    # M=m-5*np.log10(60e6)+5
    at2018 = np.load(os.path.join(at2018cow_dir,'AT2018bow.npz'))

    for fix_kappa in fix_kappa_list:
        rph_npz_fn = os.path.join(odir, f'mag_{fix_kappa}_T_{Temp/1e3:.0f}e3.npz')
        data = np.load(rph_npz_fn)
        t    = data['t']
        r_ph = data['r_ph']
        n_theta = r_ph.shape[1]
        n_t = r_ph.shape[0]

        key = [f'$\\theta = {theta}$' for theta in [10, 30, 50, 80,] ]

        fig, ax = pl.subplots(figsize=(11,7), dpi=150)
        # ax.plot(at2018['mjd']-at2018['mjd'][0]+70,at2018['mag'], 'D', label='Perley V')

        i = 1
        y = 4.64-2.5*np.log10(r_ph[:,i]/Lsun)
        # y = r_ph[:,i]/Lsun
        y[y>30]=np.nan
        # y[t<=178]=np.nan
        # y = r_ph[:,i]
        linestyle = '-'
        # if i<2:
        #     linestyle = '--'
        ax.plot(t,y, linestyle, label=key[i], linewidth=2.0)

        ylimit = ax.get_ylim()
        # ax.set_ylim(-22.0, -13.0)
        # ax.set_ylim(-26.5, -18.0)
        # ax.set_ylim(0, 9e31)
        # ax.set_ylim(ylimit[0], -16.0)
        # ax.set_ylim([-20, -17])
        ax.invert_yaxis()

        xlimit = ax.get_xlim()
        ax.set_xlim([xlimit[0], 175])
        # xlimit = [50, 110]
        # ax.set_xlim(xlimit)
        # ax.set_xticks(np.arange(xlimit[0],xlimit[1],5))

        ax.set_xlabel('t (day)', fontsize=19)
        ax.set_ylabel('Absolute Magnitude', fontsize=19)
        # ax.set_ylabel('Area ', fontsize=19)
        ax.legend(fontsize=19)
        ax.tick_params(axis='both', labelsize=19)
        # fig.savefig(os.path.join(odir, f'area_k_{fix_kappa}_T_{Temp/1e3:.0f}e3.png'))
        fig.savefig(os.path.join(odir, f'mag_spec_k_{fix_kappa}_T_{Temp/1e3:.0f}e3.png'))
        # pl.show()
        pl.close(fig)

###
if __name__ == "__main__":
    main()
    print(f'done {sys.argv[0]}')
###
