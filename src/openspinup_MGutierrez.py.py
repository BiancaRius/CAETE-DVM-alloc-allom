import pktocsv_allspins as p 
import time_series as t
import os
import joblib
import pandas as pd
import numpy as np

#### Put here the path of your folder
main_path = f'/home/amazonfaceme/biancarius/CAETE-DVM-alloc-allom/outputs/'

run_name = input('what is the run name?')

grd = '186-239' ### identify Manaus gridcell
    
path = f"../outputs/{run_name}/gridcell{grd}"

start_year = 1979
end_year   = 1990 ## ATENÇAO: data do último ano +1


run_breaks_hist1 = []
run_breaks_hist2 = []


for year in range(start_year, end_year, 1):
    #Crie as datas de início e fim no formato 'YYYYMMDD'
    start_date = f"{year}0101"
    end_date = f"{year}1231"

    print(start_date)
    print(end_date)
    
    #     # Obtenha o número do spin 
    #     spin_id = str((year - start_year) // 1 + 1).zfill(2)

    #     # Adicione a tupla à lista run_breaks_hist
    #     run_breaks_hist1.append((start_date, end_date, spin_id))


    # for year in range(start_year, end_year, 1):
    #     #Crie as datas de início e fim no formato 'YYYYMMDD'
    #     start_date = f"{year}0101"
    #     end_date = f"{year}1231"
    

    #     # Adicione a tupla à lista run_breaks_hist
    #     run_breaks_hist2.append((start_date, end_date))


    # # # Process spins 1 to ..
    # for date_range in run_breaks_hist1:

    #     start_date, end_date, spin_id = date_range
    #     file = p.read_pkz(int(spin_id), run_name, grd_name, grd_acro)
    #     p.pkz2csv(file, path, grd_name, run_name, int(spin_id), date_range, grd_acro)

    # # # Navigate to the specified folder to access spins
    # os.chdir(f'{main_path}{run_name}/{grd_name}/')

    # for date_range in run_breaks_hist2:
    #     print('Joining together all time series, dates, and spins =====',date_range)

    # print('Plotting')
    # t.join_plot(start_date, end_date, run_breaks_hist2, main_path, run_name, grd_name)



