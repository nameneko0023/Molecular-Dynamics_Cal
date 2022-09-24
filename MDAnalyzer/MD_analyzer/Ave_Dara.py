import numpy as np
import MDAnalysis as md
import os
from tqdm import tqdm
import Path_pull as PP

a = np.ones(100)
print(a)

b = 2 * np.ones(100)
print(b)

c = np.zeros(len(b))
print(c)


for i in range(3):
    c += a
    print(c)

e = 1.6 * 10**-19
debye = 3.3 * 10** -23
ps = 10 ** -12
ong = 10 ** -10

#solvation時の水分子の双極子モーメント計算================================================================================
def wat_dip(PDB, DCD):
    u = md.Universe(f"{PDB}", f"{DCD}")
    print(f"原子配列確認: {u.atoms}")
    wat_dipole = []

    #水分子の分極電荷
    H_charge = 0.4174 * e
    O_charge = -0.8340 * 10 * e

    try:
        #solvation由来のHOHと名の付く水分子
        for ts in tqdm(u.trajectory, desc="water_trajectory"):
            wat = u.select_atoms("resname HOH*")
            O = wat.select_atoms("name OH2").positions
            H1 = wat.select_atoms("name H1").positions
            H2 = wat.select_atoms("name H2").positions

            OH1 = H_charge * np.array(O - H1)
            OH2 = H_charge * np.array(O - H2)
            dipole =(ong / debye) * (OH1 + OH2)
            wat_dipole.append(dipole)

    except ValueError:
        wat_dipole = np.zeros(len(u.tranjectory))

    return wat_dipole  #水分子の双極子モーメントの3次元成分を出力


print("PDBファイルを選択してください")
PDB = PP.pdb()
print(f"選択ファイル{PDB}")
print("DCDファイルを選択してください")
DCD = PP.dcd()
print(f"選択ファイル{DCD}")
ts = float(input("シミュレーション時間刻み[ps]: "))

sub = wat_dip(PDB, DCD)
sub_name = "water"

t = len(sub)
s = len(sub[0])
scdom = s**2 - s

for x in range


