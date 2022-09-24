#座標から双極子モーメントを計算する
import os
import numpy as np
import MDAnalysis as md
from tqdm import tqdm

import Path_pull as PP



ong = 10**-10
debye = 3.33564 * 10**-30
e = 1.602176634 * 10**-19
ps = 10**-12

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
        wat_dipole = []

    return wat_dipole  #水分子の双極子モーメントの3次元成分を出力

#メタノール分子の双極子モーメント計算=====================================================================================
def met_dip(PDB, DCD):
    u = md.Universe(f"{PDB}", f"{DCD}")
    print(f"原子配列確認: {u.atoms}")
    print(u.atoms)
    met_dipole = []

    # methanolの分極電荷
    H1_charge = 0.4000 * e  # C原子周りのH原子の電荷
    O_charge = -0.6500 * e
    H46_charge = 0.0500 * e
    C_charge = 0.1000 * e

    try:
        #UNK属性のメタノール分子を抽出
        for ts in tqdm(u.trajectory, desc="methanol_trajectory"):
            met = u.select_atoms("resname UNK")
            C = met.select_atoms("name C3").positions
            O = met.select_atoms("name O2").positions
            H1 = met.select_atoms("name H1").positions
            H2 = met.select_atoms("name H4").positions
            H3 = met.select_atoms("name H5").positions
            H4 = met.select_atoms("name H6").positions

            CH = H1_charge * (np.array(C - H1) + np.array(C - H2) + np.array(C - H3))
            CO = C_charge * np.array(O - C) + H46_charge * np.array(O - H4)

            dipole = (ong / debye) * (CH + CO)
            met_dipole.append(dipole)

    except ValueError:
        met_dipole = []

    return met_dipole

 #1,2-Dichloroethane====================================================================================================
def dce12_dip(PDB, DCD):
    u = md.Universe(f"{PDB}", f"{DCD}")
    # print(f"原子配列確認: {u.atoms}")
    dce_dipole = []

    C_charge = 0.1 * e
    H_charge = 0.05 * e
    Cl_charge = -0.2 * e

    try:
        for ts in tqdm(u.trajectory, desc="1,2-dce_trajectory"):
            dce = u.select_atoms("resname UNK")
            C1 = dce.select_atoms("name C1").positions
            C2 = dce.select_atoms("name C2").positions
            H1 = dce.select_atoms("name H1").positions
            H2 = dce.select_atoms("name H2").positions
            H3 = dce.select_atoms("name H3").positions
            H4 = dce.select_atoms("name H4").positions
            Cl1 = dce.select_atoms("name Cl2").positions
            Cl2 = dce.select_atoms("name CL6").positions

            r1 = Cl_charge * np.array(C1 - Cl1) + H_charge * (np.array(C1 - H1) + np.array(C1 - H2))
            r2 = Cl_charge * np.array(C2 - Cl2) + H_charge * (np.array(C2 - H3) + np.array(C2 - H4))
            dipole = (ong / debye) * (r1 + r2)
            dce_dipole.append(dipole)

    except ValueError:
        dce_dipole = []

    return dce_dipole


#1,1,2-Trichloroethane==================================================================================================
def tce112_dip(PDB, DCD):
    u = md.Universe(f"{PDB}", f"{DCD}")
    # print(f"原子配列確認: {u.atoms}")
    tce_dipole = []

    C1_charge = 0.1 * e
    C2_charge = 0.35 * e
    H_charge = 0.05 * e
    Cl_charge = -0.2 * e

    try:
        for ts in tqdm(u.trajectory, desc="1,1,2-tce_trajecroty"):
            tce = u.select_atoms("resname UNK")
            C1 = tce.select_atoms("name C1").positions
            C2 = tce.select_atoms("name C2").positions
            H1 = tce.select_atoms("name H1").positions
            H2 = tce.select_atoms("name H2").positions
            H3 = tce.select_atoms("name H3").positions
            H4 = tce.select_atoms("name H4").positions
            Cl1 = tce.select_atoms("name Cl1").positions
            Cl2 = tce.select_atoms("name CL2").positions

            # 仮置き
            r1 = np.array(C1 - Cl1) + np.array(C1 - H1) + np.array(C1 - H2)
            r2 = np.array(C2 - Cl2) + np.array(C2 - H3) + np.array(C2 - H4)
            dipole = (ong / debye) * (r1 + r2)
            tce_dipole.append(dipole)

    except ValueError:
        tce_dipole = []

    return tce_dipole

#urea===================================================================================================================
def urea_dip(PDB, DCD):
    u = md.Universe(f"{PDB}", f"{DCD}")
    # print(f"原子配列確認: {u.atoms}")
    urea_dipole = []

    #尿素の分極電荷
    N_charge = -0.55 * e
    C_charge = 0.65 * e
    H_charge = 0.25 * e
    O_charge = -0.55 * e

    try:
        for ts in tqdm(u.trajectory, desc="urea_trajecrory"):
            urea = u.select_atoms("name UNK*")
            C = urea.select_atoms("name C2").positions
            O = urea.select_atoms("name O1").positions
            N1 = urea.select_atoms("name N3").positions
            N2 = urea.select_atoms("name N4").positions
            H1 = urea.select_atoms("name H1").positions
            H2 = urea.select_atoms("name H2").positions
            H3 = urea.select_atoms("name H3").positions
            H4 = urea.select_atoms("name H4").positions

            N1H = H_charge * (np.array(N1 - H1) + np.array(N1 - H2))
            N2H = H_charge * (np.array(N2 - H3) + np.array(N2 - H4))
            ON = N_charge * (np.array(O - N1) + np.array(O - N2))
            OC = C_charge * np.array(O - C)

            dipole = (ong / debye) * (N1H + N2H + ON + OC)
            urea_dipole.append(dipole)

    except ValueError:
        urea_dipole = []

    return urea_dipole

def gly_dip(PDB, DCD):
    u = md.Universe(f"{PDB}", f"{DCD}")
    # print(f"原子配列確認: {u.atoms}")
    gly_dipole = []

    # 尿素の分極電荷
    N_charge = -0.9000 * e
    C2_charge = 0.1000 * e
    C3_charge = 0.7000 * e
    O4_charge = -0.5500 * e
    O5_charge = -0.4000 * e
    H12_charge = 0.3500 * e
    H34_charge = 0.0500 * e
    H5_charge = 0.2500 * e

    for ts in tqdm(u.trajectory, desc="gry_trajecrory"):
        gly = u.select_atoms("resname UNK")
        N1 = gly.select_atoms("name N1").positions
        C2 = gly.select_atoms("name C2").positions
        C3 = gly.select_atoms("name C3").positions
        O4 = gly.select_atoms("name O4").positions
        O5 = gly.select_atoms("name O5").positions
        H1 = gly.select_atoms("name H1").positions
        H2 = gly.select_atoms("name H2").positions
        H3 = gly.select_atoms("name H3").positions
        H4 = gly.select_atoms("name H4").positions
        H5 = gly.select_atoms("name H5").positions

        HN = H12_charge * (np.array(N1 - H1) + np.array(N1 - H2))
        CN1 = C2_charge * np.array(N1 - C2)
        CH = H34_charge * (np.array(C2 - H3) + np.array(C2 - H4))
        CO = C3_charge * np.array(O4 - C3)
        OC = C3_charge * np.array(O5 - C3)
        HO = H5_charge * np.array(O5 - H5)

        dipole = (ong / debye) * (HN + CN1 + CH + CO + OC + HO)
        gly_dipole.append(dipole)

    return gly_dipole

def b_ala_dip(PDB, DCD):
    u = md.Universe(f"{PDB}", f"{DCD}")
    #print(f"原子配列確認: {u.atoms}")
    b_ala_dipole = []

    # 尿素の分極電荷
    N_charge = -0.9000 * e
    C2_charge = 0.1000 * e
    C3_charge = -0.1000 * e
    C4_charge = 0.7000 * e
    O5_charge = -0.5500 * e
    O6_charge = 0.4000 * e
    H12_charge = 0.3500 * e
    H36_charge = 0.0500 * e
    H7_charge = 0.2500 * e

    for ts in tqdm(u.trajectory, desc="b_ala_trajecrory"):
        gly = u.select_atoms("resname UNK")
        N1 = gly.select_atoms("name N1").positions
        C2 = gly.select_atoms("name C2").positions
        C3 = gly.select_atoms("name C3").positions
        C4 = gly.select_atoms("name C4").positions
        O5 = gly.select_atoms("name O5").positions
        O6 = gly.select_atoms("name O6").positions
        H1 = gly.select_atoms("name H1").positions
        H2 = gly.select_atoms("name H2").positions
        H3 = gly.select_atoms("name H3").positions
        H4 = gly.select_atoms("name H4").positions
        H5 = gly.select_atoms("name H5").positions
        H6 = gly.select_atoms("name H6").positions
        H7 = gly.select_atoms("name H7").positions

        HN = H12_charge * (np.array(N1 - H1) + np.array(N1 - H2))
        CN = C2_charge * np.array(N1 - C2)
        CH1 = H36_charge * (np.array(C2 - H3) + np.array(C2 - H4))
        CH2 = H36_charge * (np.array(C3 - H5) + np.array(C3 - H6))
        OC1 = C4_charge * np.array(O5 - C4)
        OC2 = C4_charge * np.array(O6 - C4)
        HO = H7_charge * np.array(O6 - H7)

        dipole = (ong / debye) * (HN + CN + CH1 + CH2 + OC1 + OC2 + HO)
        b_ala_dipole.append(dipole)

    return b_ala_dipole

if __name__=="__main__":
    print("物質名を指定してください")
    print("<<<物質: 操作番号>>>")
    print("Water: 0")
    print("Methanol: 1")
    print("1,2-Dichloroethane: 2")
    print("1,1,2-Trichloroethane: 3")
    print("urea: 4")
    print("glycine: 5")
    print("β-alanine: 6")
    ope = input("操作番号: ")

    print("PDBファイルを選択してください")
    PDB = PP.pdb()
    print(f"選択ファイル{PDB}")
    print("DCDファイルを選択してください")
    DCD = PP.dcd()
    print(f"選択ファイル{DCD}")
    ts = float(input("シミュレーション時間刻み[ps]: "))

    if ope == str(0):
        sub = wat_dip(PDB, DCD)
        sub_name = "water"
    elif ope == str(1):
        sub = met_dip(PDB, DCD)
        sub_name = "mathanol"
    elif ope == str(2):
        sub = dce12_dip(PDB, DCD)
        sub_name = "1,2-Dichloroethane"
    elif ope == str(3):
        sub = tce112_dip(PDB, DCD)
        sub_name = "1,1,2-Trichloroethane"
    elif ope == str(4):
        sub = urea_dip(PDB, DCD)
        sub_name = "urea"
    elif ope == str(5):
        sub = gly_dip(PDB, DCD)
        sub_name = "glycine"
    elif ope == str(6):
        sub = b_ala_dip(PDB, DCD)
        sub_name = "β-alanine"


    t = len(sub)
    s = len(sub[0])
    scdom = s**2 - s

    # 分母(t = 0)
    sub_czero = 0
    sub_azero = 0
    for i in tqdm(range(int(s)), desc=f"計算準備中({sub_name})"):
        for j in range(int(s)):
            if i < j:
                sub_czero += np.dot(sub[0][i], sub[0][j])  # Static auto correlation

            elif i == j:
                sub_azero += np.dot(sub[0][i], sub[0][j])  # Static cross correlation

    # 分子(t = t')
    time = []
    sub_cross = []
    sub_auto = []
    for k in tqdm(range(int(t)), desc=f"計算中({sub_name}_A&C)"):
        sub_cnum = 0
        sub_anum = 0
        for i in range(int(s)):
            for j in range(int(s)):
                if i < j:
                    sub_cnum += np.dot(sub[k][i], sub[0][j])  # cross correlation

                elif i == j:
                    sub_anum += np.dot(sub[k][i], sub[0][j])  # auto correlation

        time.append(k * ts * ps)
        sub_cross.append(sub_cnum)
        sub_auto.append(sub_anum)

    sub_num = np.array(sub_auto) / s + 2 * np.array(sub_cross) / scdom
    sub_dom = np.array(sub_azero) / s + 2 * np.array(sub_czero) / scdom
    sub_fu = np.divide(sub_num, sub_dom)

    func = np.c_[np.array(time), np.array(sub_fu)]
    fpath = os.path.dirname(PDB)

    print("", end="")
    print("計算完了")

    name = os.path.join(fpath, f"{sub_name}.csv")
    np.savetxt(
        name,
        func,
        delimiter=","
    )

    print("保存が完了しました。")
    print(f"保存ファイルパス: {name}")