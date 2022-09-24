#使用するモジュール
import tkinter as tk
from tkinter import filedialog as fd
import sys

#==========pdbファイルパスを呼び出す==========

def pdb():

    root = tk.Tk()
    root.withdraw()

    fType = [("PDBファイル", "*.pdb")]

    pdb_file = fd.askopenfilename(
        filetype=fType,
        initialdir=sys.argv[0]
    )

    root.quit()
    return pdb_file

#==========================================

#==========dcdファイルパスを呼び出す==========

def dcd():
    root = tk.Tk()
    root.withdraw()

    fType = [("DCDファイル", "*.dcd")]

    pdb_file = fd.askopenfilename(
        filetype=fType,
        initialdir=sys.argv[0]
    )

    root.quit()
    return pdb_file

#==========================================

#==========csvファイルパスを呼び出す==========

def csv():
    root = tk.Tk()
    root.withdraw()

    fType = [("CSVファイル", "*.csv")]

    pdb_file = fd.askopenfilename(
        filetype=fType,
        initialdir=sys.argv[0]
    )

    root.quit()
    return pdb_file

#==========================================

#==========TEXTファイルパスを呼び出す==========

def txt():
    root = tk.Tk()
    root.withdraw()

    fType = [("TEXTファイル", "*.txt")]

    pdb_file = fd.askopenfilename(
        filetype=fType,
        initialdir=sys.argv[0]
    )

    root.quit()
    return pdb_file

#==========================================
