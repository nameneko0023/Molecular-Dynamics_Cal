import os
import math
from tqdm import tqdm
import numpy as np

import folder_processer as fol


#生成するフォルダを格納するフォルダ名を指定する
print("保存先を指定してください")
save_dir = fol.finder()

#指定したフォルダ内にフォルダを生成する
print("格納フォルダ名を入力してください")
fname = input("フォルダ名: ")
save_path = fol.maker(save_dir, fname)


#==========モデル計算==========
cal_time = float(input("計算時間[ps]: "))
dt = float(input("時間刻み[ps]: "))
tau = float(input("緩和時間[ps]: "))
DC = float(input("DC成分: "))

data = []
i = 0
time = 0
print(f"到達数値{cal_time / dt}[ps]")

with tqdm() as pber:
    while time <= cal_time:
        time = i * dt
        fu = math.exp(-time / tau) + DC
        data.append([time, fu])
        pber.update(i)
        i += 1

name = (f"{cal_time}ps_dt{dt}ps_tau{tau}ps_DC+{DC}")
name_path = os.path.join(save_path, name)
np.savetxt(
    name_path,
    data,
    delimiter=","
)

#=============================