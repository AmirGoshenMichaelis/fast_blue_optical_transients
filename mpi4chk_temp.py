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
# import yt
# from yt import derived_field
# import numpy as np
# from matplotlib import pyplot as pl
# from astropy import constants as C
# from astropy import units as U
# Rsun = C.R_sun.cgs.value
# Msun = C.M_sun.cgs.value
# day = 24*3600
import plot_temp
###
def calc(param):
    print(f'Processing {param["fn"]} ...')
    return plot_temp.calc(param)
###
def already_exist(fn, odir ,no):
    return plot_temp.already_exist(fn, odir ,no)
###
def get_default_param():
    param = dict()
    param['chk'] = dict()
    param['chk']['dir'] = '/home/amirm/code/fast_blue_optical_transients/dd/exp3/'
    param['chk']['odir'] = '/home/amirm/code/fast_blue_optical_transients/pp_temp/exp3/'
    param['chk']['range'] = [1, 300]
    param['chk']['fn'] = 'xx_zz_dens_time_{:04}.npz'

    return param
###
def get_arg_list(param_file):
    if os.path.exists(param_file):
        with open(param_file) as f:
            param = json.load(f)
    else:
        print(f'Param file - {param_file} - is missing. Using default param.')
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
    print(f'done {sys.argv[0]} #{rank}/{size}')
###
