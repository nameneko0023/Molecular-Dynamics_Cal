#使用モジュール
import numpy as np
from tqdm import tqdm
import math
import os

#使用した自己関数
import Path_pull as PP
import folder_processer as folder
import file_processer as file


ps = 10**-12
GHz = 10**9


print("解析するファイルを選択してください")
AC_path = PP.csv()  #時間窓で分割するファイルを選択
AC = file.csv_reader(AC_path)  #選択したファイルのデータを読み込む
save_dir = os.path.dirname(AC_path)  #選択ファイルのあるディレクトリを取得
print(f"読み込んだファイル: {AC_path}")

all_time = float(input("時間範囲[ps]: "))
dt = int(input("時間窓の刻み[ps]: "))
window_path = folder.maker(save_dir, "time_window")  #選択ファイルと同じディレクトリにtime_windowフォルダを生成

i = 0
time = 0
print("", end="")
with tqdm() as pber:
    while time < all_time:
        div_time = np.where(AC[:, 0] <= i * ps)
        n = int(len(div_time[0]))
        data = AC[:n, 0:]
        wt = "{:.2f}".format(i * dt)
        name = f"data{wt}ps.csv"
        fname = os.path.join(window_path, name)
        np.savetxt(
            fname,
            data,
            delimiter=","
        )
        time = i * dt
        pber.update(dt)
        i += 1
i = 0
wt = "{:.2f}".format(i * dt)
os.remove(f"{window_path}/data{wt}ps.csv")

#==========ラプラス・フーリエ変換を連続実行==========
files = folder.lister(window_path, extention="csv")  #time_windowフォルダ内のcsvファイルをリスト化
LFT_path = folder.maker(window_path, name="MD_relaxation")  #time_window内にMD_relaxationフォルダを生成
n = int(len(files))

print("=====Start　LFT=====")
s_fre = float(input("開始周波数[GHz]: "))
e_fre = float(input("終端周波数[GHz]: "))
df = float(input("周波数刻み[GHz]: "))

print("", end="")
for k in tqdm(range(n)):
    file_path = files[k]
    csv = file.csv_reader(f"{file_path}")
    i = 0
    relax_data = []
    fre = s_fre
    m = int(len(csv))

    while fre <= e_fre:
        w = 2 * math.pi * fre * GHz
        fre = s_fre + i * df

        re = np.sum(csv[:, 1] * np.cos(w * csv[:, 0]) * csv[1][0])
        im = np.sum(csv[:, 1] * np.sin(w * csv[:, 0]) * csv[1][0])

        relax_data.append([np.log10(fre * GHz), re, im])
        i += 1

    wt = "{:.2f}".format((k + 1) * dt)
    fname = f"relax{wt}ps.csv"
    relax_path = os.path.join(LFT_path, fname)
    np.savetxt(relax_path, relax_data, delimiter=",")
print("", end="")
print("=====End　LFT=====")

#==========時間窓に対する静的誘電率, 誘電損失のピーク周波数==========
print("=====傾向分析開始=====")
files = folder.lister(LFT_path, extention="csv")
folder_path = folder.maker(LFT_path, name="analyze_data")
n = int(len(files))
analys_data = []

print("", end="")
for i in tqdm(range(n)):
    in_file = files[i]
    csv = file.csv_reader(f"{in_file}")
    m = int(len(csv))
    es = csv[0][1]
    imax = np.max(csv[:, 2])
    j = 0
    while j < m:
        if csv[j][2] == imax:
            freq = csv[j][0]
            analys_data.append([(i + 1) * dt, es, freq])
            pass
        j += 1

fname = os.path.join(folder_path, "result.csv")
np.savetxt(fname,
            analys_data,
            delimiter=","
           )
print("", end="")