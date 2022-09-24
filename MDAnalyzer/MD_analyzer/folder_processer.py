import tkinter as tk
from tkinter import filedialog as fd
import natsort as na
import os
import glob
import sys

#==========フォルダパスを取得する================================================================
def finder():
    root = tk.Tk()
    root.withdraw()

    fpath = fd.askdirectory(
        initialdir=sys.argv[0]
    )
    root.quit()
    return fpath

#=============================================================================================

#====フォルダ生成===============================================================================
def maker(folder_dir, name):
    fname = os.path.join(folder_dir, f"{name}")
    os.makedirs(fname)
    return fname
#=============================================================================================

#====フォルダ内のファイルパス一括取得==============================================================

def lister(folder_dir, extention):
    files = folder_dir + f"/*.{extention}"
    get_paths = [glob.glob(files)]
    file_path = na.natsorted(get_paths[0])
    return file_path

#==============================================================================================