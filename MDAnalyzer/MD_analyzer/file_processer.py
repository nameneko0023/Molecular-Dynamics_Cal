from tqdm import tqdm
import numpy as np
import os

import folder_processer as folder

#==========csvファイルを読み込む==========

def csv_reader(file_path):
    rfile = np.loadtxt(
        fname=file_path,
        dtype=float,  #float型として読み込む
        delimiter=",",  #コンマ区切りのデータとして読み込む
    )
    return rfile

#========================================

#==========textファイルを読み込む==========

def txt_reader(file_path):
    rfile = np.loadtxt(
        fname=file_path,
        dtype=float,
        delimiter="\n",  #タブ区切りのデータとして読み込む
        skiprows=0
    )
    return rfile

#========================================