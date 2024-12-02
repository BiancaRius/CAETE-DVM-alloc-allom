import os
import shutil
import multiprocessing as mp
from pathlib import Path
import joblib
import numpy as np
from post_processing import write_h5
import h52nc


pls_number = '3000' #input('how many PLSs? ')
run_name = 'base_run' #input('run name: ')


while True:
    server = input('Are you running in the server? y/n ')

    if server == 'y':
        main_path = '/home/amazonfaceme/biancarius/CAETE-DVM-alloc-allom/'
        break
    if server == 'n':
        main_path = '/home/bianca/bianca/CAETE-DVM-alloc-allom/'
        break

start_cond_filename = f'CAETE_STATE_START_{run_name}_.pkz'

attr_filename = f'{main_path}/outputs/{run_name}/gridcell175-235/pls_attrs-{pls_number}.csv'

run_path = Path(f"{main_path}/outputs/{run_name}/{start_cond_filename}")
pls_path = Path(f"{main_path}/outputs/{run_name}/{attr_filename}")

with open(run_path, 'rb') as fh:
    init_conditions = joblib.load(fh)

new_run_name = input('what will be the name of your new run? ')

# new outputs folder
dump_folder = Path(f"{new_run_name}")

for gridcell in init_conditions:
    gridcell.clean_run(dump_folder, "init_cond")
    gridcell.pr -= gridcell.pr * 0.5
    # prevent negative values
    gridcell.pr[np.where(gridcell.pr < 0.0)[0]] = 0.0
    assert np.all(gridcell.pr >= 0.0)

from caete import run_breaks_hist as rb

def zip_gridtime(grd_pool, interval):
    res = []
    for i, j in enumerate(grd_pool):
        res.append((j, interval[i % len(interval)]))
    return res


def apply_funX(grid, brk):
    grid.run_caete(brk[0], brk[1])
    return grid

n_proc = mp.cpu_count()

for i, brk in enumerate(rb):
    print(f"Applying model to the interval {brk[0]}-{brk[1]}")
    init_conditions = zip_gridtime(init_conditions, (brk,))




