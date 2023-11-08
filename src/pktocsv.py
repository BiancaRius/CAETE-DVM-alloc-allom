import os
import joblib
from pathlib import Path
from caete_module import global_par as gp
import caete_module as mod
import _pickle as pkl
import bz2
from numpy import asfortranarray
from pandas import read_csv

import cfunits
import numpy as np
import pandas as pd
import cftime as cf

import plsgen as pls
from aux_plot import get_var

#get the number of PLSs
npls = input('how many PLSs?')

#get pls table
pls_table = pd.read_csv("../outputs/pls_attrs-{}.csv".format(npls))

#get the run name
run_name = input('run name: ')

grd = input("Which gridcell (lat-long)")

#get the run path
path = f"../outputs/{run_name}/gridcell{grd}"

grd_name = f"gridcell{grd}"

# Navigate to the folder
os.chdir(path)

#get the spin
spin = input('which spin?')

#get the pkz file
def read_pkz(spin):
    with open("spin{}.pkz".format(spin), 'rb') as fh:
        dt = joblib.load(fh)
    print(dt.keys()) # list the available keys for the ouputs
    return dt

file = read_pkz(spin)


def pkz2csv(file, path, grd_name):
    print(path)
    print(grd_name)

pkz2csv(file, path, grd_name)



def apply_spin(grid):
    """pre-spinup use some outputs of daily budget (water, litter C, N and P) to start soil organic pools"""
    w, ll, cwd, rl, lnc = grid.bdg_spinup(
        start_date="20000102", end_date="20050102")
    grid.sdc_spinup(w, ll, cwd, rl, lnc)
    return grid


# def get_spin(grd: mod.plot, spin) -> dict:
#     import joblib
#     if spin < 10:
#         name = f'spin0{spin}.pkz'
#     else:
#         name = f'spin{spin}.pkz'
#     with open(grd.outputs[name], 'rb') as fh:
#         spin_dt = joblib.load(fh)
#     return spin_dt

# spin_dt = get_spin(grd, spin)
# Get the number of PLSs


pkz = '/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/t1_0510/gridcell175-235/spin19.pkz'
output_folder = '/home/bianca/bianca/CAETE-alloc-allom/outputs/t1_0510/gridcell175-235/'

def pkz2csv(file_path, csv_file):
    # Criar o diretório se não existir
    os.makedirs(Path(output_folder), exist_ok=True)

    with open(file_path, 'rb') as fh:
        dt = joblib.load(fh)
    
    df = pd.DataFrame(dt)

    # Construir o caminho para o arquivo CSV no diretório especificado
    csv_path = Path(output_folder) / csv_file
    df.to_csv(csv_path, index=False)

# Usar a função
pkz2csv(pkz, 't.csv')

# import os
# import joblib
# import pickle as pkl
# import copy
# import bz2
# from pathlib import Path
# import multiprocessing as mp
# import pandas as pd
# import numpy as np 
# import joblib
# from caete_module import global_par as gp
# import caete_module as mod

# #get the number of PLSs
# npls = gp.npls

# pkz = '/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/t1_0510/gridcell175-235/spin19.pkz'
# output_folder = '/home/bianca/bianca/CAETE-alloc-allom/outputs/t1_0510/gridcell175-235/output_folder/'

# def pkz2csv(file_path, csv_file):

#     # Criar o diretório se não existir
#     os.makedirs(Path(file_path).parent, exist_ok=True)

#     with open(pkz, 'rb') as fh:
#        dt = joblib.load(fh)
    
#     df = pd.DataFrame(dt)

#     csv_path = output_folder + csv_file
#     df.to_csv(csv_path, index=False)


# t = pkz2csv(pkz, 't.csv')


    
#     #    print(dt.keys()) # list the available keys for the ouputs
#     #    # 'calendar', 'time_unit', 'sind', 'eind'
#     #    print(dt['calendar'])
#     #    print(dt['time_unit'])
#     #    print(dt['sind'])
#     #    print(dt['eind'])

# # pkz_folder = '/home/bianca/bianca/CAETE-DVM-alloc-allom/outputs/t1_0510/gridcell175-235/'
# # output_folder = '/home/bianca/bianca/CAETE-alloc-allom/outputs/t1_0510/gridcell175-235/output_folder'

# def pkz2csv(pkz, rpath, mod, scen):

#     from cftime import num2pydate
#     spin_dt = read_pkz(pkz)
#     s = num2pydate(spin_dt['sind'], spin_dt['time_unit'], spin_dt['calendar'])
#     e = num2pydate(spin_dt['eind'], spin_dt['time_unit'], spin_dt['calendar'])
#     start = s.strftime("%Y%m%d")
#     end   = e.strftime("%Y%m%d") 
#     idxT1 = pd.date_range(start, end, freq='D', closed=None)
#     idx = idxT1.to_series()


 
