import os
import shutil
from glob import glob

while True:
    grd_acro = input('Gridcell acronym [AFL, ALP, FEC, MAN, CAX, NVX]: ')

    if grd_acro == 'ALP':
        grd = '188-213'
        break
    elif grd_acro == 'FEC':
        grd = '200-225'
        break
    elif grd_acro == 'MAN':
        grd = '186-239'
        break
    elif grd_acro == 'CAX':
        grd = '183-257'
        break
    elif grd_acro == 'NVX':
        grd = '210-249'
        break
    elif grd_acro == 'AFL':
        grd = '199-248'
        break
    else:
        print('This acronym does not correspond')
        break

def mover_csv(origem, destino, grd_acro, grd):
    # Certifique-se de que a pasta de destino exista
    if not os.path.exists(destino):
        os.makedirs(destino)

    # Cria o padrão do nome do arquivo usando o prefixo
    padrao_nome_arquivo = f"concatenated_series_{grd_acro}_10prec_7y.csv"

    caminho_origem = os.path.join(origem, padrao_nome_arquivo)

    # Obtém a lista de arquivos correspondentes ao padrão
    arquivos = glob(caminho_origem)

    if not arquivos:
        print(f'Nenhum arquivo encontrado em {caminho_origem}')
        return

    # Move os arquivos para a pasta de destino
    for arquivo in arquivos:
        novo_caminho = os.path.join(destino, os.path.basename(arquivo))
        shutil.move(arquivo, novo_caminho)
        print(f'Movendo {arquivo} para {novo_caminho}')

    print('Operação concluída com sucesso.')

if __name__ == "__main__":
    origem = f'/home/amazonfaceme/biancarius/CAETE-DVM-alloc-allom/outputs/{grd_acro}/experiments/10perc_reduction/{grd_acro}_10prec_7y/gridcell{grd}/'
    destino = f'/home/amazonfaceme/biancarius/CAETE-DVM-alloc-allom/outputs/{grd_acro}/experiments/10perc_reduction'

    mover_csv(origem, destino, grd_acro, grd)
