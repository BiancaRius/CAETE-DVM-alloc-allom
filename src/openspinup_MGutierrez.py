
import pandas as pd
import joblib

# Caminho para o arquivo .pkz
path = f"../outputs"
run_name = "test_spinup_marcela3"
grd = '186-239'
start_year = 1979
end_year = 1989

# Abrir o arquivo PKZ
with open(f"{path}/{run_name}/gridcell{grd}/spin35.pkz", 'rb') as fh:
    dt = joblib.load(fh)

# Variável de interesse
variables = ['photo', 'ar', 'npp', 'lai', 'f5', 'evapm',
                     'cleaf', 'cwood', 'croot', 'csap', 'cheart',
                     'csto', 'wue', 'ls']

# Gerar a sequência de datas
date_index = pd.date_range(start=f"{start_year}-01-01", end=f"{end_year}-12-31", freq='D')
num_days = len(date_index)

# Inicializa um dicionário com a coluna de datas
data_dict = {'Date': date_index}

# Adiciona cada variável ao dicionário
for var in variables:
    series = dt.get(var, [])
    if len(series) != num_days:
        raise ValueError(f"A variável '{var}' tem {len(series)} valores, mas o esperado é {num_days}.")
    data_dict[var] = series

# Cria o DataFrame
df = pd.DataFrame(data_dict)


# Salvar o DataFrame em um arquivo CSV
output_path = f"{path}/{run_name}/gridcell{grd}/npp_data_{run_name}.csv"
df.to_csv(output_path, index=False)

print(f"Arquivo CSV salvo em: {output_path}")


