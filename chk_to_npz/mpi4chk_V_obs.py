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
# import h5py
# yt.enable_parallelism()
# import pyxsim
# import soxs
# from scipy import optimize
###
import getopt
import os
import sys
import json
import math
import itertools
from mpi4py import MPI 
import yt
from yt import derived_field
import numpy as np
# from matplotlib import pyplot as pl
from astropy import constants as C
from astropy import units as U
Rsun = C.R_sun.cgs.value
Msun = C.M_sun.cgs.value
day = 24*3600
# import V_obs_calc
import V_obs_calc_resample as V_obs_calc
###
def calc(param):
    print('Processing {} ...'.format(param["fn"]))
    return V_obs_calc.calc(param)
###
def already_exist(fn, odir ,no):
    return V_obs_calc.already_exist(fn, odir ,no)
###
def get_default_param():
    param = dict()
    param['chk'] = dict()
    param['chk']['dir'] = '/data/home/chelouche/makashi/FLASH_lobes/object'
    param['chk']['odir'] = '/data/home/chelouche/makashi/FLASH_lobes/observ'
    param['chk']['range'] = [1, 260]
    param['chk']['fn'] = 'PN1_hdf5_chk_{:04}'

    param['Opal_Table'] = dict()
    param['Opal_Table']['low']  = '/data/home/chelouche/makashi/FLASH_lobes/observ/lowT_fa05_gs98_z0.02_x0.7.data.json'
    param['Opal_Table']['high'] = '/data/home/chelouche/makashi/FLASH_lobes/observ/gs98_z0.02_x0.7.data.json'
    param['Opal_Table']['T'] = 6.0e3

    param['theta'] = list(range(0,91,10))
    param['phi']   = [0.,]

#    param['fix_kappa_list'] = [0.06, 0.3, 1.5]
    param['fix_kappa_list'] = [0.06, 0.1, 0.2, 0.3, 1.5]

    return param
###
def get_arg_list(param_file):
    if os.path.exists(param_file):
        with open(param_file) as f:
            param = json.load(f)
    else:
        print('Param file - {} - is missing. Using default param.'.format(param_file))
        param = get_default_param() 
    arg_list = list()
    for i in range(param['chk']['range'][0],param['chk']['range'][1]+1):
        _param = dict(param)
        _param['no.'] = i
        _param['fn'] = os.path.join(_param['chk']['dir'],_param['chk']['fn'].format(_param['no.']))
        if not os.path.exists(_param['fn']):
            # print(f'File {fn_str} does not exists! continue ...')
            continue
        if already_exist(_param['fn'],_param['chk']['odir'],_param['no.']):
            print('file {} already exist. skiping file ...'.format(_param['fn']))
            continue
        arg_list.append( (_param,) )
    return arg_list
###
def print_usage():
    print("Usage:  {} --param or -p parameter file in json format".format(sys.argv[0]))
    print("")
    print("Usage: {} --help or -h print this help".format(sys.argv[0]))
###
def main():
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank == 0:
        param_file = ''
        try:
            optlist, args = getopt.getopt(sys.argv[1:], 'p:h', ['param=', 'help'])
        except getopt.GetoptError as err:
            print(err)
            print_usage()
            sys.exit(2)
        for opt, val in optlist:
            if opt in ("-p", "--param"):
                param_file = val
            elif opt in ("-h", "--help"):
                print_usage()
                return
            else:
                assert False, "unhandled option"
        arg_list = get_arg_list(param_file)

        m = int(math.ceil(float(len(arg_list))/size))
        for chunk in range(1,size):
            comm.send(arg_list[chunk*m:(chunk+1)*m], dest=chunk)
        local_arg_list = arg_list[0:m]
    else:
        local_arg_list = comm.recv(source=0)
    
    res = list(map(lambda x: calc(*x), local_arg_list))
    allres = comm.gather(res, root=0)
    if rank==0:
        allres = list(itertools.chain.from_iterable(allres))
        print('post processing all results:')
        print(allres)
###
if __name__ == "__main__":
    main()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    print('done {} #{}/{}'.format(sys.argv[0],rank,size))
###

