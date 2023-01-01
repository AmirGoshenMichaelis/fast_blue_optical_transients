#! /usr/bin/python3
###
# import glob
# import inspect
# import importlib
# import re
# from multiprocessing import Pool
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# from matplotlib import rc
# from astropy import time
# from scipy import interpolate
# from scipy import signal
# from scipy import stats
# import pyxsim
# import soxs
# from scipy import optimize
###
import enum
import getopt
import os
import sys
import json
import math
import itertools
from mpi4py import MPI 
import yt
# yt.enable_parallelism()
from yt import derived_field
import numpy as np
# from matplotlib import pyplot as pl
from astropy import constants as C
from astropy import units as U
# import h5py
Rsun = C.R_sun.cgs.value
Msun = C.M_sun.cgs.value
day  = 24*3600
Lsun = C.L_sun.cgs.value
sigma_sb = C.sigma_sb.cgs.value
###
class LineOfSight:
    _sl = None
    _lp = None
    los = None
    p1  = None
    r   = None
    def __init__(self, ds, theta, phi, p1) -> None:
        # https://yt-project.org/doc/analyzing/objects.html
        # https://yt-project.org/doc/_modules/yt/data_objects/selection_objects/ray.html
        # θ = theta*np.pi/180.0
        # φ = phi*np.pi/180.
        # self.los = np.array([np.sin(θ)*np.cos(φ),np.sin(θ)*np.sin(φ),np.cos(θ)])
        self.los = np.array([np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)])
        self.p1 = p1
        # R = np.sqrt(np.sum(ds.domain_right_edge.value**2))
        R = np.max(ds.domain_right_edge.value)
        start = self.p1-R*self.los
        end   = self.p1+R*self.los

        lp = np.array([start+dl*self.los for dl in np.linspace(0,2*R,1024)])
        pp = (lp[:,0]<=ds.domain_right_edge.value[0]) & \
             (lp[:,1]<=ds.domain_right_edge.value[1]) & \
             (lp[:,2]<=ds.domain_right_edge.value[2]) & \
             (lp[:,0]>=ds.domain_left_edge.value[0]) &  \
             (lp[:,1]>=ds.domain_left_edge.value[1]) &  \
             (lp[:,2]>=ds.domain_left_edge.value[2])
        if not np.any(pp):
            self.r   = []
            self.xyz = []
        start = lp[pp,:][0]
        end   = lp[pp,:][-1]
        # ds.domain_right_edge.value, ds.domain_left_edge.value
        self._sl = ds.r[start:end]
        X  = self._sl['flash','x'].value
        Y  = self._sl['flash','y'].value
        Z  = self._sl['flash','z'].value

        xyz = np.array([X-start[0],Y-start[1],Z-start[2]]).transpose()
        rxyz = np.linalg.norm(xyz,axis=1)

        self._lp = np.argsort(rxyz)

        self.r = rxyz[self._lp]
        self.xyz = np.array([X[self._lp],Y[self._lp],Z[self._lp]]).transpose()

    def sl(self,src):
        return np.array(self._sl[src].value)[self._lp]
###
class Kappa:
    opal = None
    def __init__(self, opal_table):
        with open(opal_table, 'r') as fh:
            opal = json.load(fh)
        self.opal = opal
    def __call__(self, rho, temp):
        # k = 0.2*(1+0.7) # 0.2(1+X)cm^2/gr
        logRho = np.log10(rho)
        logT   = np.log10(temp)
        logR = logRho - 3.0*logT + 18.0
        indR = np.argmin(np.abs(self.opal['logR']-logR))
        indT = np.argmin(np.abs(self.opal['logT']-logT))
        # print(f" Rho={rho:.2e} T={temp:.2e} kappa={10.**self.opal['table'][indT][indR]:.2e} logR={logR:.2e} logRho={logRho:.2e} logT={logT:.2e} logKappa={self.opal['table'][indT][indR]:.2e}")
        return 10.**self.opal['table'][indT][indR]
###
def already_exist(fn, odir, no):
    return False
    # ofn = os.path.join(odir, f'TRL_eff_{no:04}.npz')
    # if os.path.exists(ofn):
    #     return True
    # return False
###
def p1(alpha,theta,phi):
    p1 = np.array([np.cos(alpha)*np.cos(theta)*np.cos(phi) - np.sin(alpha)*np.sin(phi),
               np.sin(alpha)*np.cos(phi) + np.cos(alpha)*np.cos(theta)*np.sin(phi),
               -np.cos(alpha)*np.sin(theta)]) #.reshape((3,1))
    return p1
###
def calc_fk(param):
    ds = yt.load(param['fn'])

    all_obs = list()

    param['theta'] = np.array(param['theta'])
    param['phi'] = np.array(param['phi'])
    
    param['R'] = np.linspace(3.*Rsun, np.min(ds.domain_right_edge.value), 100)
    param['alpha_rad'] = np.array([0.,np.pi])
    
    ni = len(param['theta'])
    nj = len(param['phi'])
    nk = len(param['R'])*len(param['alpha_rad'])

    Teff = np.zeros((ni,nj,nk))*np.nan
    Reff = np.zeros((ni,nj,nk))*np.nan
    Leff = np.zeros((ni,nj,nk))*np.nan
    XYZeff  = np.zeros((ni,nj,nk,3))*np.nan
    p1_list = np.zeros((ni,nj,nk,3))*np.nan
    vel_eff = np.zeros((ni,nj,nk,3))*np.nan

    time = np.array([ds.current_time.to('d').value,])

    theta_rad = param['theta']*np.pi/180.
    phi_rad = param['phi']*np.pi/180.

    dtheta = np.pi/len(param['alpha_rad'])
    dphi   = np.pi/len(param['alpha_rad']) 

    # dtheta_list = np.diff(theta_rad)
    # dtheta_list = np.append(dtheta_list[0],dtheta_list) if theta_rad.size>1 else np.array([np.pi/180.,])

    # dphi_list = np.diff(phi_rad)
    # dphi_list = np.append(dphi_list[0],dphi_list) if phi_rad.size>1 else np.array([np.pi/180.,])

    # opal_high = Kappa(param['Opal_Table']['high'])
    # opal_low  = Kappa(param['Opal_Table']['low'])

    # opacity_units = np.array(['COND_VAR [cm^2/s]', 'OPAC_VAR[cm^2/gr]', 'trans_opac [1/cm]', 'kappa=c/3/COND_VAR[1/cm]', 'opac=c/3/COND_VAR/DENS_VAR[cm^2/gr]', 'RadTrans.F90 line 305'], dtype='S40')
    # for i,theta,dtheta in zip(range(ni),theta_rad,dtheta_list):
    #     for j,phi,dphi in zip(range(nj),phi_rad,dphi_list):
    fix_kappa = param['kappa']
    for i,theta in zip(range(ni),theta_rad):
        for j,phi in zip(range(nj),phi_rad):
            k = -1
            for rr in param['R']:
                for alpha in param['alpha_rad']:
                    k = k + 1
                    p_shift = rr*p1(alpha,theta,phi)
                    ray = LineOfSight(ds,theta,phi,p_shift)

                    if ray.r.size<2:
                        print('ray size is to small {} at: theta {} phi {} R {} alpha {}'.format(ray.r.size,theta*180./np.pi,phi*180./np.pi,rr/Rsun,alpha*180./np.pi))
                        continue

                    temp  = ray.sl( ('flash', 'temp') )
                    dens  = ray.sl( ('flash', 'dens') )
                    # opac  = ray.sl( ('flash', 'opac') ) # in units of cm^2/gr
                    # size_of_opac = len(dens)
                    # opac = np.zeros((size_of_opac,))
                    # for rho,T,ind_op in zip(dens,temp,range(size_of_opac)):
                    #     if T < param['Opal_Table']['T']:
                    #         opac[ind_op] = opal_low(rho, T)
                    #     else:
                    #         opac[ind_op] = opal_high(rho, T)
                    # kappa = opac*dens # in units of 1/cm

                    # ind_op_de = calc_optical_depth(kappa, ray.r)
                    ind_op_de = calc_optical_depth(fix_kappa*dens, ray.r)
                    if np.isnan(ind_op_de):
                        print('Opacity not reach value of 2/3 at: theta {} phi {} R {} alpha {}'.format(theta*180./np.pi,phi*180./np.pi,rr/Rsun,alpha*180./np.pi))
                        continue

                    p1_list[i,j,k,:]=p_shift
                    Teff[i,j,k] = temp[ind_op_de]
                    Reff[i,j,k] = np.sqrt(np.sum(ray.xyz[ind_op_de,:]**2))
                    XYZeff[i,j,k,:] = ray.xyz[ind_op_de,:]
                    vel_eff[i,j,k,0] = ray.sl( ('flash', 'velx') )[ind_op_de]
                    vel_eff[i,j,k,1] = ray.sl( ('flash', 'vely') )[ind_op_de]
                    vel_eff[i,j,k,2] = ray.sl( ('flash', 'velz') )[ind_op_de]
                    theta_xyz = np.arccos(ray.xyz[ind_op_de,2]/Reff[i,j,k])
                    cos_theta_xyz = np.dot(ray.xyz[ind_op_de,:],ray.los)/(np.linalg.norm(ray.los)*np.linalg.norm(ray.xyz[ind_op_de,:]))
                    if cos_theta_xyz <= 0:
                        Leff[i,j,k] = 0.
                    else:
                        Leff[i,j,k] = sigma_sb * Teff[i,j,k]**4 * Reff[i,j,k]**2 * np.sin(theta_xyz)*dtheta*dphi * cos_theta_xyz

    ofn = os.path.join(param['chk']['odir'], 'TRL_eff_{:04}_fix_kappa_{}.npz'.format(param['no.'], fix_kappa))
    np.savez(ofn, time=time, Teff=Teff, Reff=Reff, Leff=Leff, XYZeff=XYZeff, p1=p1_list, vel_eff=vel_eff, theta=param['theta'], phi=param['phi']) 
    return 'Write {}'.format(ofn)
###
def calc(param):
    for param['kappa'] in param['fix_kappa_list']:
        calc_fk(param)
###
def calc_optical_depth(kappa, r):
    dr = np.diff(r)
    dr = np.append(dr[0],dr)
    opt_dep = kappa*dr
    opacity = np.cumsum(opt_dep[::-1])
    ind = np.sum(opacity<=(2./3.))
    if ind>=opacity.size :
        return np.nan
    return opacity.size-1-ind

    # if ind>=opacity.size :
    #     ind=opacity.size-1
    # return opacity[::-1], opacity.size-1-ind
###
def main():
    param = dict()
    param['chk'] = dict()
    param['chk']['dir'] = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/'
    param['chk']['odir'] = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/'
    param['chk']['fn'] = 'PN1_hdf5_chk_{:04}'
    param['no.'] = 64 # 61
    param['fn'] = os.path.join(param['chk']['dir'],param['chk']['fn'].format(param['no.']))

    param['Opal_Table'] = dict()
    param['Opal_Table']['low']  = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/lowT_fa05_gs98_z0.02_x0.7.data.json'
    param['Opal_Table']['high'] = '/home/amirm/code/flash_disk_v4.6.2/obj-res/akashi/gs98_z0.02_x0.7.data.json'
    param['Opal_Table']['T'] = 6.5e3

    param['theta'] = list(range(0,91,10)) # [30., 70.,]
    param['phi']   = [0.,]

    param['fix_kappa_list'] = [0.06, 0.3, 1.5]

    os.chdir(param['chk']['dir'])
    calc(param)
###
if __name__ == "__main__":
    main()
    print('done {}'.format(sys.argv[0]))
###
###
