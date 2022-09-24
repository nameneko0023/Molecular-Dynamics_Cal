import tkinter as tk
from tkinter import filedialog as fd
import sys
from tqdm import tqdm
import os
import math
import numpy as np
import MDAnalysis as md
from time import sleep


# =====トポロジーファイル(pdb, crd)ファイルパスを呼び出す=====
def top():
    root = tk.Tk()
    root.withdraw()

    file_name = [("PDBファイル", "*.pdb"), ("CRDファイル", "*.crd")]

    top_file = fd.askopenfilename(
        filetype=file_name,
        # initialdir=sys.argv[0]
    )

    root.quit()
    return top_file


# =====dcdファイルパスを呼び出す=====
def dcd():
    root = tk.Tk()
    root.withdraw()

    file_name = [("DCDファイル", "*.dcd")]

    dcd_file = fd.askopenfilename(
        filetype=file_name,
        # initialdir=sys.argv[0]
    )
    root.quit()
    return dcd_file


# =====単位定義=====
ong = 10 ** -10
D = 3.33564 * 10 ** -30
e = 1.602176634 * 10 ** -19
ps = 10 ** -12


# =====DS使用時ファイル解析==============================================================================================
# =====Water_molecule=====
def WAT_dip(u):
    global OH2
    wat_dipole = []

    # 水分子の分極電荷
    H_charge = 0.4174 * e
    O_charge = -0.8340 * e

    # =====solvation由来のHOHと名の付く水分子構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="water_trajectory"):
        wat = u.select_atoms("resname HOH*")
        O = wat.select_atoms("name OH2").positions
        H1 = wat.select_atoms("name H1").positions
        H2 = wat.select_atoms("name H2").positions

        OH1 = H_charge * np.array(O - H1) * (ong / D)
        OH2 = H_charge * np.array(O - H2) * (ong / D)

        wat_dip = OH1 + OH2
        wat_dipole.append(wat_dip)
    sleep(0.05)

    print("=====" + "水分子の双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"O-H1: {math.sqrt(OH1[0][0] * OH1[0][0] + OH1[0][1] * OH1[0][1] + OH1[0][2] * OH1[0][2])}")
    print(f"O-H2: {math.sqrt(OH2[0][0] * OH2[0][0] + OH2[0][1] * OH2[0][1] + OH2[0][2] * OH2[0][2])}")
    print(
        f"water: {math.sqrt(wat_dipole[0][0][0] * wat_dipole[0][0][0] + wat_dipole[0][0][1] * wat_dipole[0][0][1] + wat_dipole[0][0][2] * wat_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return wat_dipole  # 水分子の双極子モーメントの3次元成分を出力


# =====Methanol_molecule=====
def MeOH_dip(u):
    met_dipole = []

    # =====methanolの部分電荷=====
    H1_charge = 0.4000 * e
    O_charge = -0.6500 * e
    H46_charge = 0.0500 * e
    C_charge = 0.1000 * e

    # =====UNK属性のメタノール構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="methanol_trajectory"):
        met = u.select_atoms("resname UNK*")
        C = met.select_atoms("name C3").positions
        O = met.select_atoms("name O2").positions
        H1 = met.select_atoms("name H1").positions
        H4 = met.select_atoms("name H4").positions
        H5 = met.select_atoms("name H5").positions
        H6 = met.select_atoms("name H6").positions

        CH_CO = (H46_charge * (np.array(C - H4) + np.array(C - H5) + np.array(C - H6)) + C_charge * np.array(O - C)) * (
                ong / D)
        CH = H1_charge * np.array(O - H1) * (ong / D)

        dipole = CH_CO + CH
        met_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "methanoalの双極子モーメント(永久双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"CH_CO: {math.sqrt(CH_CO[0][0] * CH_CO[0][0] + CH_CO[0][1] * CH_CO[0][1] + CH_CO[0][2] * CH_CO[0][2])}")
    print(f"CH: {math.sqrt(CH[0][0] * CH[0][0] + CH[0][1] * CH[0][1] + CH[0][2] * CH[0][2])}")
    print(
        f"methanol: {math.sqrt(met_dipole[0][0][0] * met_dipole[0][0][0] + met_dipole[0][0][1] * met_dipole[0][0][1] + met_dipole[0][0][2] * met_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return met_dipole


# =====Ethanol_molecule=====
def EtOH_dip(u):
    eto_dipole = []

    # =====エタノールの分極電荷=====
    C1_charge = -0.1500 * e  # C原子周りのH原子の電荷
    C2_charge = 0.1500 * e
    O3_charge = -0.6500 * e
    H15_charge = 0.0500 * e
    H6_charge = 0.4000 * e

    # =====UNK属性のエタノール構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="ethanol_trajectory"):
        eto = u.select_atoms("resname UNK")
        C1 = eto.select_atoms("name C1").positions
        C2 = eto.select_atoms("name C2").positions
        O3 = eto.select_atoms("name O3").positions
        H1 = eto.select_atoms("name H1").positions
        H2 = eto.select_atoms("name H2").positions
        H3 = eto.select_atoms("name H3").positions
        H4 = eto.select_atoms("name H4").positions
        H5 = eto.select_atoms("name H5").positions
        H6 = eto.select_atoms("name H6").positions

        C1_O = (H15_charge * (np.array(C1 - H1) + np.array(C1 - H2) + np.array(C1 - H3))) * (ong / D)
        CH_OH = (H15_charge * (np.array(C2 - H4) + np.array(C2 - H5)) + H6_charge * np.array(O3 - H6)) * (ong / D)

        dipole = C1_O + CH_OH
        eto_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "methanoalの双極子モーメント(永久双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"C1_O: {math.sqrt(C1_O[0][0] * C1_O[0][0] + C1_O[0][1] * C1_O[0][1] + C1_O[0][2] * C1_O[0][2])}")
    print(f"CH: {math.sqrt(CH_OH[0][0] * CH_OH[0][0] + CH_OH[0][1] * CH_OH[0][1] + CH_OH[0][2] * CH_OH[0][2])}")
    print(
        f"methanol: {math.sqrt(eto_dipole[0][0][0] * eto_dipole[0][0][0] + eto_dipole[0][0][1] * eto_dipole[0][0][1] + eto_dipole[0][0][2] * eto_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return eto_dipole


# =====Discovery Studio データの解析=====
def DS(u):
    print("==========物質名一覧==========")
    print("Water")
    print("Ethanol")
    print("Methanol")
    print("Ethan")
    print("1,2-DCE")
    print("1,1,2-TCE")

    print("", end="\n")
    print("物質名を入力してください")
    ope = input("物質名: ")

    if ope == "Water":
        sol = WAT_dip(u)

    return ope, sol


# ====計算サーバー使用時ファイル解析======================================================================================
# =====Water_molecule=====
def WAT_dipS(u):
    wat_dipole = []
    # 水分子の分極電荷
    H_charge = 0.4174 * e
    O_charge = -0.8340 * e
    # solvation由来のHOHと名の付く水分子
    for ts in tqdm(u.trajectory, desc="water_trajectory"):
        wat = u.select_atoms("resname TIP*")
        O = wat.select_atoms("name OH2").positions
        H1 = wat.select_atoms("name H1").positions
        H2 = wat.select_atoms("name H2").positions
        OH1 = H_charge * np.array(O - H1) * (ong / D)
        OH2 = H_charge * np.array(O - H2) * (ong / D)
        wat_dip = OH1 + OH2
        wat_dipole.append(wat_dip)
    sleep(0.05)

    print("=====" + "水分子の双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"O-H1: {math.sqrt(OH1[0][0] * OH1[0][0] + OH1[0][1] * OH1[0][1] + OH1[0][2] * OH1[0][2])}")
    print(f"O-H2: {math.sqrt(OH2[0][0] * OH2[0][0] + OH2[0][1] * OH2[0][1] + OH2[0][2] * OH2[0][2])}")
    print(
        f"water: {math.sqrt(wat_dipole[0][0][0] * wat_dipole[0][0][0] + wat_dipole[0][0][1] * wat_dipole[0][0][1] + wat_dipole[0][0][2] * wat_dipole[0][0][2])}")
    print("=" * 117)
    return wat_dipole  # 水分子の双極子モーメントの3次元成分を出力


# =====Methanol_molecule=====
def MeOH_dipS(u):
    met_dipole = []

    # =====methanolの部分電荷=====
    H1_charge = 0.4000 * e  # C原子周りのH原子の電荷
    O_charge = -0.6500 * e
    H46_charge = 0.0500 * e
    C_charge = 0.1000 * e

    # =====R01属性のメタノール分子の構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="methanol_trajectory"):
        met = u.select_atoms("resname R01")
        C = met.select_atoms("name C3").positions
        O = met.select_atoms("name O2").positions
        H1 = met.select_atoms("name H1").positions
        H4 = met.select_atoms("name H4").positions
        H5 = met.select_atoms("name H5").positions
        H6 = met.select_atoms("name H6").positions

        CH_CO = (H46_charge *
                 (np.array(C - H4) + np.array(C - H5) + np.array(C - H6))
                 + C_charge * np.array(O - C)) * (ong / D)

        CH = H1_charge * np.array(O - H1) * (ong / D)

        dipole = CH_CO + CH
        met_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "methanoalの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"CH_CO: {math.sqrt(CH_CO[0][0] * CH_CO[0][0] + CH_CO[0][1] * CH_CO[0][1] + CH_CO[0][2] * CH_CO[0][2])}")
    print(f"CH: {math.sqrt(CH[0][0] * CH[0][0] + CH[0][1] * CH[0][1] + CH[0][2] * CH[0][2])}")
    print(
        f"methanol: {math.sqrt(met_dipole[0][0][0] * met_dipole[0][0][0] + met_dipole[0][0][1] * met_dipole[0][0][1] + met_dipole[0][0][2] * met_dipole[0][0][2])}")
    print("=" * 117)
    return met_dipole


# =====Ethanol_molecule=====
def EtOH_dipS(u):
    eto_dipole = []

    # =====methanolの部分電荷=====
    C1_charge = -0.1500 * e  # C原子周りのH原子の電荷
    C2_charge = 0.1500 * e
    O3_charge = -0.6500 * e
    H15_charge = 0.0500 * e
    H6_charge = 0.4000 * e

    # =====R01属性のエタノール分子の構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="methanol_trajectory"):
        eto = u.select_atoms("resname R01")
        C1 = eto.select_atoms("name C1").positions
        C2 = eto.select_atoms("name C2").positions
        O3 = eto.select_atoms("name O3").positions
        H1 = eto.select_atoms("name H1").positions
        H2 = eto.select_atoms("name H2").positions
        H3 = eto.select_atoms("name H3").positions
        H4 = eto.select_atoms("name H4").positions
        H5 = eto.select_atoms("name H5").positions
        H6 = eto.select_atoms("name H6").positions

        C1_O = (H15_charge * (np.array(C1 - H1) + np.array(C1 - H2) + np.array(C1 - H3))) * (ong / D)
        CH_OH = (H15_charge * (np.array(C2 - H4) + np.array(C2 - H5)) + H6_charge * np.array(O3 - H6)) * (ong / D)

        dipole = C1_O + CH_OH
        eto_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "methanoalの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"C1_O: {math.sqrt(C1_O[0][0] * C1_O[0][0] + C1_O[0][1] * C1_O[0][1] + C1_O[0][2] * C1_O[0][2])}")
    print(f"CH: {math.sqrt(CH_OH[0][0] * CH_OH[0][0] + CH_OH[0][1] * CH_OH[0][1] + CH_OH[0][2] * CH_OH[0][2])}")
    print(
        f"methanol: {math.sqrt(eto_dipole[0][0][0] * eto_dipole[0][0][0] + eto_dipole[0][0][1] * eto_dipole[0][0][1] + eto_dipole[0][0][2] * eto_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return eto_dipole


# =====Ethane_molecule=====
def Et_dipS(u):
    eth_dipole = []

    # =====ethanの部分電荷=====
    C_charge = -0.15 * e
    H_charge = 0.05 * e

    # =====R01属性のエタン分子の構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="ethane_trajectory"):
        eth = u.select_atoms("resname R01")
        C1 = eth.select_atoms("name C1").positions
        C2 = eth.select_atoms("name C2").positions
        H1 = eth.select_atoms("name H1").positions
        H2 = eth.select_atoms("name H2").positions
        H3 = eth.select_atoms("name H3").positions
        H4 = eth.select_atoms("name H4").positions
        H5 = eth.select_atoms("name H5").positions
        H6 = eth.select_atoms("name H6").positions

        C1_dip = H_charge * (np.array(C1 - H4) + np.array(C1 - H5) + np.array(C1 - H6)) * (ong / D)
        C2_dip = H_charge * (np.array(C2 - H1) + np.array(C1 - H2) + np.array(C1 - H3)) * (ong / D)
        dipole = C1_dip + C2_dip

        eth_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "ethanの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"C1-H: {math.sqrt(C1_dip[0][0] * C1_dip[0][0] + C1_dip[0][1] * C1_dip[0][1] + C1_dip[0][2] * C1_dip[0][2])}")
    print(f"C2-H: {math.sqrt(C2_dip[0][0] * C2_dip[0][0] + C2_dip[0][1] * C2_dip[0][1] + C2_dip[0][2] * C2_dip[0][2])}")
    print(
        f"ethanoal: {math.sqrt(eth_dipole[0][0][0] * eth_dipole[0][0][0] + eth_dipole[0][0][1] * eth_dipole[0][0][1] + eth_dipole[0][0][2] * eth_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return eth_dipole


# 1,2-Dichloroethane=====
def DCE12_dipS(u):
    dce12_dipole = []

    # =====1,2-Dichloroethaneの部分電荷=====
    C_charge = 0.1 * e
    H_charge = 0.05 * e
    Cl_charge = -0.2 * e

    # =====R01属性の1,2-ジクロロエタン分子の構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="1,2-dce_trajectory"):
        dce = u.select_atoms("resname R01")
        C1 = dce.select_atoms("name C1").positions
        C2 = dce.select_atoms("name C2").positions
        H1 = dce.select_atoms("name H1").positions
        Cl2 = dce.select_atoms("name CL2").positions
        H3 = dce.select_atoms("name H3").positions
        Cl4 = dce.select_atoms("name CL4").positions
        H5 = dce.select_atoms("name H5").positions
        H6 = dce.select_atoms("name H6").positions

        C1_Cl1 = (Cl_charge * np.array(C1 - Cl2) + H_charge * (np.array(C1 - H1) + np.array(C1 - H3))) * (ong / D)
        C2_Cl2 = (Cl_charge * np.array(C2 - Cl4) + H_charge * (np.array(C2 - H5) + np.array(C2 - H6))) * (ong / D)
        dipole = C1_Cl1 + C2_Cl2
        dce12_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "ethanの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(
        f"C1-Cl1: {math.sqrt(C1_Cl1[0][0] * C1_Cl1[0][0] + C1_Cl1[0][1] * C1_Cl1[0][1] + C1_Cl1[0][2] * C1_Cl1[0][2])}")
    print(
        f"C2-Cl2: {math.sqrt(C2_Cl2[0][0] * C2_Cl2[0][0] + C2_Cl2[0][1] * C2_Cl2[0][1] + C2_Cl2[0][2] * C2_Cl2[0][2])}")
    print(
        f"dce: {math.sqrt(dce12_dipole[0][0][0] * dce12_dipole[0][0][0] + dce12_dipole[0][0][1] * dce12_dipole[0][0][1] + dce12_dipole[0][0][2] * dce12_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return dce12_dipole


# 1,1,2-Trichloroethane=====
def TCE112_dipS(u):
    tce112_dipole = []

    # =====1,1,2-Trichloroethaneの部分電荷=====
    C1_charge = 0.35 * e
    C2_charge = 0.1 * e
    H_charge = 0.05 * e
    Cl_charge = -0.2 * e

    # =====R01属性の1,1,2-トリクロロエタン分子の構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="1,1,2-tce_trajecroty"):
        tce = u.select_atoms("resname R01")
        C1 = tce.select_atoms("name C1").positions
        C2 = tce.select_atoms("name C2").positions
        H1 = tce.select_atoms("name H1").positions
        CL2 = tce.select_atoms("name CL2").positions
        H3 = tce.select_atoms("name H3").positions
        CL4 = tce.select_atoms("name CL4").positions
        CL5 = tce.select_atoms("name CL5").positions
        H6 = tce.select_atoms("name H6").positions

        C1_Cl = (Cl_charge * (np.array(C1 - CL4) + np.array(C1 - CL5)) + H_charge * np.array(C1 - H6)) * (ong / D)
        C2_Cl = (Cl_charge * np.array(C2 - CL2) + H_charge * (np.array(C2 - H1) + np.array(C2 - H3))) * (ong / D)
        dipole = C1_Cl + C2_Cl
        tce112_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "1,1,2-Trichloroethaneの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"C1-Cl1: {math.sqrt(C1_Cl[0][0] * C1_Cl[0][0] + C1_Cl[0][1] * C1_Cl[0][1] + C1_Cl[0][2] * C1_Cl[0][2])}")
    print(f"C2-Cl2: {math.sqrt(C2_Cl[0][0] * C2_Cl[0][0] + C2_Cl[0][1] * C2_Cl[0][1] + C2_Cl[0][2] * C2_Cl[0][2])}")
    print(
        f"tce: {math.sqrt(tce112_dipole[0][0][0] * tce112_dipole[0][0][0] + tce112_dipole[0][0][1] * tce112_dipole[0][0][1] + tce112_dipole[0][0][2] * tce112_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return tce112_dipole


# =====Hexachloroethane_molecule=====
def HCE_dipS(u):
    hce_dipole = []

    # =====Hexachloroethane_moleculeの部分電荷=====
    C_charge = 0.6 * e
    Cl_charge = -0.2 * e

    # =====R01属性のヘキサクロロエタン分子の構成原子を抽出=====
    for ts in tqdm(u.trajectory, desc="Hexachloroethane_trajecroty"):
        tce = u.select_atoms("resname R01")
        C1 = tce.select_atoms("name C1").positions
        C2 = tce.select_atoms("name C2").positions
        CL3 = tce.select_atoms("name CL3").positions
        CL4 = tce.select_atoms("name CL4").positions
        CL5 = tce.select_atoms("name CL5").positions
        CL6 = tce.select_atoms("name CL6").positions
        CL7 = tce.select_atoms("name CL7").positions
        CL8 = tce.select_atoms("name CL8").positions

        C1_Cl = Cl_charge * (np.array(C1 - CL6) + np.array(C1 - CL7) + np.array(C1 - CL8)) * (ong / D)
        C2_Cl = Cl_charge * (np.array(C2 - CL3) + np.array(C2 - CL4) + np.array(C2 - CL5)) * (ong / D)
        dipole = C1_Cl + C2_Cl
        hce_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "1,1,2-Trichloroethaneの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"C1-Cl1: {math.sqrt(C1_Cl[0][0] * C1_Cl[0][0] + C1_Cl[0][1] * C1_Cl[0][1] + C1_Cl[0][2] * C1_Cl[0][2])}")
    print(f"C2-Cl2: {math.sqrt(C2_Cl[0][0] * C2_Cl[0][0] + C2_Cl[0][1] * C2_Cl[0][1] + C2_Cl[0][2] * C2_Cl[0][2])}")
    print(
        f"tce: {math.sqrt(hce_dipole[0][0][0] * hce_dipole[0][0][0] + hce_dipole[0][0][1] * hce_dipole[0][0][1] + hce_dipole[0][0][2] * hce_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return hce_dipole


# =====acethon_molecule=====
def ACE_dipS(u):
    ace_dipole = []

    # =====acethonの部分電荷=====
    O_charge = 0.55 * e
    C13_charge = -0.15 * e
    C2_charge = 0.55 * e
    H_charge = 0.05 * e

    # =====R01属性のアセトン分子を抽出=====
    for ts in tqdm(u.trajectory, desc="acethon_trajectory"):
        ace = u.select_atoms("resname R01")
        C1 = ace.select_atoms("name C1").positions
        C2 = ace.select_atoms("name C2").positions
        C3 = ace.select_atoms("name C3").positions
        O4 = ace.select_atoms("name O4").positions
        H1 = ace.select_atoms("name H1").positions
        H2 = ace.select_atoms("name H2").positions
        H3 = ace.select_atoms("name H3").positions
        H4 = ace.select_atoms("name H4").positions
        H5 = ace.select_atoms("name H5").positions
        H6 = ace.select_atoms("name H6").positions

        H_C1 = H_charge * (np.array(C1 - H1) + np.array(C1 - H2) + np.array(C1 - H3)) * (ong / D)
        H_C2 = H_charge * (np.array(C1 - H4) + np.array(C1 - H5) + np.array(C1 - H6)) * (ong / D)
        O_H = O_charge * np.array(C1 - O4) * (ong / D)

        dipole = H_C1 + H_C2 + O_H
        ace_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "acetonの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"H_C1: {math.sqrt(H_C1[0][0] * H_C1[0][0] + H_C1[0][1] * H_C1[0][1] + H_C1[0][2] * H_C1[0][2])}")
    print(f"H_C2: {math.sqrt(H_C2[0][0] * H_C2[0][0] + H_C2[0][1] * H_C2[0][1] + H_C2[0][2] * H_C2[0][2])}")
    print(f"O_H: {math.sqrt(O_H[0][0] * O_H[0][0] + O_H[0][1] * O_H[0][1] + O_H[0][2] * O_H[0][2])}")
    print(
        f"aceton: {math.sqrt(ace_dipole[0][0][0] * ace_dipole[0][0][0] + ace_dipole[0][0][1] * ace_dipole[0][0][1] + ace_dipole[0][0][2] * ace_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return ace_dipole


# =====DMSO_molecule=====
def DMSO_dipS(u):
    dim_dipole = []

    # =====DMSOの部分電荷=====
    C_charge = -0.0858 * e
    S_charge = 0.7647 * e
    H_charge = 0.0077 * e
    O_charge = -0.6213 * e

    # =====R01属性のDMSO分子を抽出=====
    for ts in tqdm(u.trajectory, desc="DMSO_trajectory"):
        dim = u.select_atoms("resname R01")
        C1 = dim.select_atoms("name C1").positions
        S2 = dim.select_atoms("name S2").positions
        C3 = dim.select_atoms("name C3").positions
        H1 = dim.select_atoms("name H1").positions
        H2 = dim.select_atoms("name H2").positions
        H3 = dim.select_atoms("name H3").positions
        H4 = dim.select_atoms("name H4").positions
        H5 = dim.select_atoms("name H5").positions
        H6 = dim.select_atoms("name H6").positions
        O10 = dim.select_atoms("name O10").positions

        H_C1 = H_charge * (np.array(C1 - H1) + np.array(C1 - H2) + np.array(C1 - H3)) * (ong / D)
        H_C3 = H_charge * (np.array(C3 - H4) + np.array(C3 - H5) + np.array(C3 - H6)) * (ong / D)
        O_H = S_charge * np.array(S2 - O10) * (ong / D)

        dipole = H_C1 + H_C3 + O_H
        dim_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "DMSOの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"H_C1: {math.sqrt(H_C1[0][0] * H_C1[0][0] + H_C1[0][1] * H_C1[0][1] + H_C1[0][2] * H_C1[0][2])}")
    print(f"H_C3: {math.sqrt(H_C3[0][0] * H_C3[0][0] + H_C3[0][1] * H_C3[0][1] + H_C3[0][2] * H_C3[0][2])}")
    print(f"O_H: {math.sqrt(O_H[0][0] * O_H[0][0] + O_H[0][1] * O_H[0][1] + O_H[0][2] * O_H[0][2])}")
    print(
        f"DNSO: {math.sqrt(dim_dipole[0][0][0] * dim_dipole[0][0][0] + dim_dipole[0][0][1] * dim_dipole[0][0][1] + dim_dipole[0][0][2] * dim_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return dim_dipole


# =====Urea=====
def urea_dipS(u):
    urea_dipole = []

    # Ureaの分極電荷
    O_charge = -0.55 * e
    C_charge = 0.65 * e
    N_charge = 0.55 * e
    H_charge = 0.25 * e

    for ts in tqdm(u.trajectory, desc="Urea_trajecrory"):
        ure = u.select_atoms("resname R01")
        O1 = ure.select_atoms("name O1").positions
        C1 = ure.select_atoms("name C1").positions
        N3 = ure.select_atoms("name N3").positions
        N4 = ure.select_atoms("name N4").positions
        H1 = ure.select_atoms("name H1").positions
        H2 = ure.select_atoms("name H2").positions
        H3 = ure.select_atoms("name H3").positions
        H4 = ure.select_atoms("name H4").positions

        C_N3 = (C_charge * np.array(N3 - C1) + H_charge * (np.array(N3 - H1) + np.array(N3 - H2))) * (ong / D)
        C_N4 = (C_charge * np.array(N4 - C1) + H_charge * (np.array(N4 - H3) + np.array(N4 - H4))) * (ong / D)
        C_O = (C_charge * np.array(O1 - C1)) * (ong / D)
        dipole = C_N3 + C_N4 + C_O

        urea_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "ureaの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print(f"C-N3: {math.sqrt(C_N3[0][0] * C_N3[0][0] + C_N3[0][1] * C_N3[0][1] + C_N3[0][2] * C_N3[0][2])}")
    print(f"C-O4: {math.sqrt(C_N4[0][0] * C_N4[0][0] + C_N4[0][1] * C_N4[0][1] + C_N4[0][2] * C_N4[0][2])}")
    print(f"C-O4: {math.sqrt(C_O[0][0] * C_O[0][0] + C_O[0][1] * C_O[0][1] + C_O[0][2] * C_O[0][2])}")
    print(
        f"urea: {math.sqrt(urea_dipole[0][0][0] * urea_dipole[0][0][0] + urea_dipole[0][0][1] * urea_dipole[0][0][1] + urea_dipole[0][0][2] * urea_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return urea_dipole


# 1,4-dioxane=まだ=====
def DIO14_dipS(u):
    dio14_dipole = []

    # 1,4-ジオキサンの分極電荷
    C_charge = 0.08 * e
    H_charge = 0.05 * e
    O_charge = -0.36 * e

    for ts in tqdm(u.trajectory, desc="1,4dioxane_trajecrory"):
        dio14 = u.select_atoms("name R01")
        O1 = dio14.select_atoms("name O1").positions
        C2 = dio14.select_atoms("name C2").positions
        C3 = dio14.select_atoms("name C3").positions
        O4 = dio14.select_atoms("name O4").positions
        C5 = dio14.select_atoms("name C5").positions
        C6 = dio14.select_atoms("name C6").positions
        H1 = dio14.select_atoms("name H1").positions
        H2 = dio14.select_atoms("name H2").positions
        H3 = dio14.select_atoms("name H3").positions
        H4 = dio14.select_atoms("name H4").positions
        H5 = dio14.select_atoms("name H5").positions
        H6 = dio14.select_atoms("name H6").positions
        H7 = dio14.select_atoms("name H7").positions
        H8 = dio14.select_atoms("name H8").positions

        CO1 = C_charge * (np.array(O1 - C2) + np.array(O1 - C6)) * (ong / D)
        CO4 = C_charge * (np.array(O4 - C3) + np.array(O4 - C5)) * (ong / D)
        C1H = H_charge * (np.array(C2 - H1) + np.array(C2 - H2)) * (ong / D)
        C3H = H_charge * (np.array(C3 - H3) + np.array(C3 - H4)) * (ong / D)
        C5H = H_charge * (np.array(C5 - H5) + np.array(C5 - H6)) * (ong / D)
        C6H = H_charge * (np.array(C6 - H7) + np.array(C6 - H8)) * (ong / D)

        dipole = CO1 + CO4 + C1H + C3H + C5H + C6H
        dio14_dipole.append(dipole)
    sleep(0.05)

    print("=====" + "1,4-dioの双極子モーメント(永久双極子モーメント　+ 誘起双極子モーメント)" + "=" * 50)
    print("※数値を使用する場合はメモをする")
    print(f"C-O1: {math.sqrt(CO1[0][0] * CO1[0][0] + CO1[0][1] * CO1[0][1] + CO1[0][2] * CO1[0][2])}")
    print(f"C-O4: {math.sqrt(CO4[0][0] * CO4[0][0] + CO4[0][1] * CO4[0][1] + CO4[0][2] * CO4[0][2])}")
    print(
        f"1,4-dio: {math.sqrt(dio14_dipole[0][0][0] * dio14_dipole[0][0][0] + dio14_dipole[0][0][1] * dio14_dipole[0][0][1] + dio14_dipole[0][0][2] * dio14_dipole[0][0][2])}")
    print("=" * 117, end="\n\n")
    return dio14_dipole


def server(u):
    print("==========物質名一覧==========")
    print("Water")
    print("Ethanol")
    print("Methanol")
    print("Ethan")
    print("1,2-DCE")
    print("1,1,2-TCE")
    print("Hexachloroethane")
    print("Acethon")
    print("Urea")
    print("DMSO")
    print()

    print("", end="\n")
    print("物質名を入力してください")
    ope = input("物質名: ")

    if ope == "Water":
        sol = WAT_dipS(u)
    elif ope == "Methanol":
        sol = MeOH_dipS(u)
    elif ope == "Ethanol":
        sol = Et_dipS(u)
    elif ope == "Ethan":
        sol = Et_dipS(u)
    elif ope == "1,2-DCE":
        sol = DCE12_dipS(u)
    elif ope == "1,1,2-TCE":
        sol = TCE112_dipS(u)
    elif ope == "Hexachloroethane":
        sol = HCE_dipS(u)
    elif ope == "Urea":
        sol = urea_dipS(u)
    elif ope == "DMSO":
        sol = DMSO_dipS(u)
    elif ope == "Acethon":
        sol = ACE_dipS(u)
    elif ope == "1,4-dioxane":
        sol = DIO14_dipS(u)

    return ope, sol


# =====純物質溶液系の自己・相互相関関数=====================================================================================
def ac(top, ope, mate):
    t = len(mate)
    s = len(mate[0])
    scdom = s ** 2 - s

    ts = float(input("Save Result Interval[ps]: "))
    print("時間の積分範囲を指定してください[ps]")
    scope = float(input("積分範囲0ps ~ :"))

    cut = scope // ts  # 分割した次のスタート地点をコマ数に換算
    n = t // cut  # 全体時間に対して何回積算できるのかを計算

    div = input("積算計算する(yes or no)：　")

    if div == "yes":
        print(f"Number of average {n}回で計算を実行します", end="\n\n")

        # 時間軸作成
        time_axis = np.arange(0, cut * ts, ts)

        # 積算用ゼロ配列
        sfunc = np.zeros(int(cut))  # A&Ccorrelation積算用0行列
        sfunc_auto = np.zeros(int(cut))  # Auto-correlation積算用0行列
        sfunc_cross = np.zeros(int(cut))  # Cross-correlation積算用0行列

        # 積算ループ
        for x in range(int(n)):
            y = x + 1
            print(f"{y}回目")
            sleep(0.05)

            # 計算範囲の指定
            cut_sub = mate[x * int(cut):(x + 1) * int(cut)]

            # 自己・相互相関関数の分母計算
            czero = 0
            azero = 0
            for i in tqdm(range(s), desc="計算準備中"):
                for j in range(s):

                    if i == j:
                        azero += np.dot(cut_sub[0][i], cut_sub[0][j])

                    if i < j:
                        czero += np.dot(cut_sub[0][i], cut_sub[0][j])
            sleep(0.05)

            # 自己・相互相関関数の分子計算
            sauto = []
            scross = []
            for t in tqdm(range(int(cut)), desc="計算中"):
                anum = 0
                cnum = 0
                for i in range(s):
                    for j in range(s):

                        if i == j:
                            anum += np.dot(cut_sub[t][i], cut_sub[0][j])
                        if i < j:
                            cnum += np.dot(cut_sub[t][i], cut_sub[0][j])

                sauto.append(anum)
                scross.append(cnum)
            sleep(0.05)

            # 自己・相互相関関数
            fnum = np.array(sauto) / s + (2 * np.array(scross)) / scdom
            fdom = azero / s + (2 * czero) / scdom
            func = fnum / fdom

            # 自己相関関数
            fnum_auto = np.array(sauto) / s
            fdom_auto = azero / s
            func_auto = fnum_auto / fdom_auto

            # 相互相関関数
            fnum_cross = (2 * np.array(scross)) / scdom
            fdom_cross = (2 * czero) / scdom
            func_cross = fnum_cross / fdom_cross

            # 計算結果を積算する
            sfunc += np.array(func)
            sfunc_auto += np.array(func_auto)
            sfunc_cross += np.array(func_cross)
        sleep(0.05)
        print("=" * 100, end="\n")

        # 積算したAuto- and Cross- correlation の平均値を計算
        ave_func = np.array(sfunc) / n
        ave_auto = np.array(sfunc_auto) / n
        ave_cross = np.array(sfunc_cross) / n

        # 平均値を時間軸を統合する
        ave = np.c_[time_axis, ave_func]
        ave_A = np.c_[time_axis, ave_auto]
        ave_C = np.c_[time_axis, ave_cross]

        # topファイルの格納先ディレクトリを取得
        dir = os.path.dirname(top)

        # topファイルと同じフォルダへ格納する
        fn = f"{ope}_{n}average.csv"
        name = os.path.join(dir, fn)

        np.savetxt(
            name,
            ave,
            delimiter=","
        )

        print("完了")
        print(name)

        fn_a = f"{ope}_Auto_{n}average.csv"
        name = os.path.join(dir, fn_a)

        np.savetxt(
            name,
            ave_A,
            delimiter=","
        )

        print("完了")
        print(name)

        fn_c = f"{ope}_Cross_{n}average.csv"
        name = os.path.join(dir, fn_c)

        np.savetxt(
            name,
            ave_C,
            delimiter=","
        )

        print("完了")
        print(name)

    elif div == "no":

        # 時間軸作成
        time_axis = np.arange(0, cut * ts, ts)

        cut_sub = mate[0:int(cut)]

        # 自己・相互相関関数の分母計算
        czero = 0
        azero = 0
        for i in tqdm(range(s), desc="計算準備中"):
            for j in range(s):

                if i == j:
                    azero += np.dot(cut_sub[0][i], cut_sub[0][j])

                if i < j:
                    czero += np.dot(cut_sub[0][i], cut_sub[0][j])
        sleep(0.05)

        # 自己・相互相関関数の分子計算
        sauto = []
        scross = []
        for t in tqdm(range(int(cut)), desc="計算中"):
            anum = 0
            cnum = 0
            for i in range(s):
                for j in range(s):

                    if i == j:
                        anum += np.dot(cut_sub[t][i], cut_sub[0][j])
                    if i < j:
                        cnum += np.dot(cut_sub[t][i], cut_sub[0][j])

            sauto.append(anum)
            scross.append(cnum)
        sleep(0.05)

        # 自己・相互相関関数
        fnum = np.array(sauto) / s + (2 * np.array(scross)) / scdom
        fdom = azero / s + (2 * czero) / scdom
        func = fnum / fdom
        data1 = np.c_[time_axis, func]

        # 自己相関関数
        fnum_auto = np.array(sauto) / s
        fdom_auto = azero / s
        func_auto = fnum_auto / fdom_auto
        data2 = np.c_[time_axis, func_auto]

        # 相互相関関数
        fnum_cross = (2 * np.array(scross)) / scdom
        fdom_cross = (2 * czero) / scdom
        func_cross = fnum_cross / fdom_cross
        data3 = np.c_[time_axis, func_cross]

        dir = os.path.dirname(top)

        fn = f"{ope}.csv"
        name = os.path.join(dir, fn)

        np.savetxt(
            name,
            data1,
            delimiter=","
        )

        print("完了")
        print(name)

        fn_a = f"{ope}_Auto.csv"
        name = os.path.join(dir, fn_a)

        np.savetxt(
            name,
            data2,
            delimiter=","
        )

        print("完了")
        print(name)

        fn_c = f"{ope}_Cross.csv"
        name = os.path.join(dir, fn_c)

        np.savetxt(
            name,
            data3,
            delimiter=","
        )

        print("完了")
        print(name)


# =====混合物質溶液系の自己・相互相関関数===================================================================================
def acaq(top, ope, wat, sub):
    WAT = 18.0152

    if ope == "Ethanol":
        Mw = 46.0684
    elif ope == "Methanol":
        Mw = 16.0424
    elif ope == "Ethan":
        Mw = 30.0690
    elif ope == "1,2-DCE":
        Mw = 98.9585
    elif ope == "1,1,2-TCE":
        Mw = 133.4033
    elif ope == "Urea":
        Mw = 60.055
    elif ope == "DMSO":
        Mw = 78.1334
    elif ope == "1,4-Dioxane":
        Mw = 88.1051
    elif ope == "Acethon":
        Mw = 96.5085
    elif ope == "Glycine":
        Mw = 75.0666
    elif ope == "β-Alanine":
        Mw = 89.0932
    else:
        print("""
        <<=====================>>
        分子名が間違っています。
        もう一度やり直してください。
        <<=====================>>
        """)
        pass

    t = len(wat)
    w = len(wat[0])
    wcdom = w * w - w
    s = len(sub[0])
    scdom = s * s - s

    ts = float(input("Save Result Interval[ps]: "))
    print("時間の積分範囲を指定してください[ps]")
    scope = float(input("積分範囲0ps ~ :"))

    cut = scope // ts  # 分割した次のスタート地点をコマ数に換算
    n = t // cut  # 全体時間に対して何回積算できるのかを計算

    div = input("積算計算する(yes or no)：　")

    if div == "yes":
        print(f"Number of average {n}回で計算を実行します", end="\n\n")

        # 時間軸作成
        time_axis = np.arange(0, cut * ts, ts)

        mix = np.zeros(int(cut))  # A&Ccorrelation積算用0行列
        solvent = np.zeros(int(cut))  # Auto-correlation積算用0行列
        solute = np.zeros(int(cut))  # Cross-correlation積算用0行列

        for x in range(int(n)):
            y = x + 1
            print(f"{y}回目")
            sleep(0.05)

            cut_wat = wat[x * int(cut):(x + 1) * int(cut)]

            # 自己・相互相関関数の分母計算
            wat_czero = 0
            wat_azero = 0
            for i in tqdm(range(int(w)), desc=f"計算準備中(water)"):
                for j in range(int(w)):
                    if i < j:
                        wat_czero += np.dot(cut_wat[0][i], cut_wat[0][j])  # Static auto correlation

                    elif i == j:
                        wat_azero += np.dot(cut_wat[0][i], cut_wat[0][j])  # Static cross correlation
            sleep(0.05)

            # 自己・相互相関関数の分母計算
            wat_cross = []
            wat_auto = []
            for k in tqdm(range(int(cut)), desc="計算中(water_A&C)"):
                wat_cnum = 0
                wat_anum = 0
                for i in range(int(w)):
                    for j in range(int(w)):
                        if i < j:
                            wat_cnum += np.dot(cut_wat[k][i], cut_wat[0][j])  # cross correlation

                        elif i == j:
                            wat_anum += np.dot(cut_wat[k][i], cut_wat[0][j])  # auto correlation

                wat_cross.append(wat_cnum)
                wat_auto.append(wat_anum)
                sleep(0.05)
            sleep(0.05)

            # ==========溶質の自己・相互相関==========

            cut_sol = sub[x * int(cut):(x + 1) * int(cut)]

            # 自己・相互相関関数の分母計算
            sol_czero = 0
            sol_azero = 0
            for i in tqdm(range(int(s)), desc=f"計算準備中({ope})"):
                for j in range(int(s)):
                    if i < j:
                        sol_czero += np.dot(cut_sol[0][i], cut_sol[0][j])  # Static auto correlation

                    elif i == j:
                        sol_azero += np.dot(cut_sol[0][i], cut_sol[0][j])  # Static cross correlation
                sleep(0.05)
            sleep(0.05)

            # 自己・相互相関関数の分子計算
            sol_cross = []
            sol_auto = []
            for k in tqdm(range(int(cut)), desc=f"計算中({ope}_A&C)"):
                # 積分範囲を更新するごとに分母を0にリセットする
                sol_cnum = 0
                sol_anum = 0
                for i in range(int(s)):
                    for j in range(int(s)):
                        if i < j:
                            sol_cnum += np.dot(cut_sol[k][i], cut_sol[0][j])  # cross correlation

                        elif i == j:
                            sol_anum += np.dot(cut_sol[k][i], cut_sol[0][j])  # auto correlation

                sol_cross.append(sol_cnum)
                sol_auto.append(sol_anum)
            sleep(0.05)

            # ==========水分子と溶質分子の相互相関==========
            # 分母(t = 0)
            mix_czero = 0
            for i in tqdm(range(int(s)), desc=f"計算準備中({ope}aq)"):
                for j in range(int(w)):
                    mix_czero += np.dot(cut_sol[0][i], wat[0][j])  # Static auto correlation
            sleep(0.05)

            # 分子(t = t')
            mix_cross = []
            for k in tqdm(range(int(cut)), desc=f"計算中({ope}aq_A&C)"):
                mix_cnum = 0
                for i in range(int(s)):
                    for j in range(int(w)):
                        mix_cnum += np.dot(cut_sol[k][i], wat[0][j])  # cross correlation
                mix_cross.append(mix_cnum)
            sleep(0.05)

            # ==========永久双極子モーメントの自己・相互相関関数==========
            # ==水==
            wat_num = np.array(wat_auto) / w + 2 * np.array(wat_cross) / wcdom
            wat_dom = np.array(wat_azero) / w + 2 * np.array(wat_czero) / wcdom

            # ==溶質==
            sol_num = np.array(sol_auto) / s + 2 * np.array(sol_cross) / scdom
            sol_dom = np.array(sol_azero) / s + 2 * np.array(sol_czero) / scdom

            # ==水-溶質==
            mix_num = np.array(mix_cross) / (w * s)
            mix_dom = np.array(mix_czero) / (w * s)

            mix += (wat_num + sol_num + mix_num) / (wat_dom + sol_dom + mix_dom)
            solvent += wat_num / wat_dom
            solute += sol_num / sol_dom

        # 積算した混合水溶液、水、溶質のAuto- and Cross- correlation の平均値を計算
        mix_ave = mix / n
        solvent_ave = solvent / n
        solute_ave = solute / n

        func_wat = np.c_[np.array(time_axis), np.array(solvent_ave)]
        func_sol = np.c_[np.array(time_axis), np.array(solute_ave)]
        func = np.c_[np.array(time_axis), np.array(mix_ave)]

        file_path = os.path.dirname(top)
        print("計算完了")

        # 重量パーセント濃度
        wt = (Mw * s) / (Mw * s + WAT * w) * 100

        # 水溶液の水 Auto- and Cross-Correlation
        name = os.path.join(file_path, f"{ope}aq_solvent_{n}ave.csv")
        np.savetxt(
            name,
            func_wat,
            delimiter=","
        )

        # 水溶液の溶質A&C
        name = os.path.join(file_path, f"{ope}aq_substans_{n}ave.csv")
        np.savetxt(
            name,
            func_sol,
            delimiter=","
        )

        # 水溶液のA&C
        name = os.path.join(file_path, f"{ope}{wt:.2f}wt%_{n}ave.csv")
        np.savetxt(
            name,
            func,
            delimiter=","
        )

        print("保存が完了しました。")
        print(f"ファイルパス: {name}")

        # モル比率濃度
        mol = s / (s + w) * 100

        print("=" * 100)
        print("数値を使用する場合はメモをする")
        # 小数点第2位まで表示する
        print(f"{ope}の重量濃度は{wt:.2f}wt%です")
        print(f"{ope}のモル比率は{mol:.2f}mol%です")

    elif div == "no":

        # 時間軸作成
        time_axis = np.arange(0, cut * ts, ts)

        # 時間積分を指定
        cut_wat = wat[0:int(cut)]

        # 相互相関関数の分母計算
        wat_czero = 0
        wat_azero = 0
        for i in tqdm(range(int(w)), desc=f"計算準備中(water)"):
            for j in range(int(w)):
                if i < j:
                    wat_czero += np.dot(cut_wat[0][i], cut_wat[0][j])  # Static auto correlation
                elif i == j:
                    wat_azero += np.dot(cut_wat[0][i], cut_wat[0][j])  # Static cross correlation
        sleep(0.05)

        # 相互相関関数の分子計算
        wat_cross = []
        wat_auto = []
        for k in tqdm(range(int(cut)), desc="計算中(water_A&C)"):
            wat_cnum = 0
            wat_anum = 0
            for i in range(int(w)):
                for j in range(int(w)):
                    if i < j:
                        wat_cnum += np.dot(cut_wat[k][i], cut_wat[0][j])  # cross correlation
                    elif i == j:
                        wat_anum += np.dot(cut_wat[k][i], cut_wat[0][j])  # auto correlation
            wat_cross.append(wat_cnum)
            wat_auto.append(wat_anum)
        sleep(0.05)

        # ==========溶質の自己・相互相関==========
        # 時間積分範囲を指定
        cut_sol = sub[0:int(cut)]

        # 相互相関関数の分母計算
        sol_czero = 0
        sol_azero = 0
        for i in tqdm(range(int(s)), desc=f"計算準備中({ope})"):
            for j in range(int(s)):
                if i < j:
                    sol_czero += np.dot(cut_sol[0][i], cut_sol[0][j])  # Static auto correlation
                elif i == j:
                    sol_azero += np.dot(cut_sol[0][i], cut_sol[0][j])  # Static cross correlation
        sleep(0.05)

        # 相互相関関数の分子計算
        sol_cross = []
        sol_auto = []
        for k in tqdm(range(int(cut)), desc=f"計算中({ope}_A&C)"):
            sol_cnum = 0
            sol_anum = 0
            for i in range(int(s)):
                for j in range(int(s)):
                    if i < j:
                        sol_cnum += np.dot(cut_sol[k][i], cut_sol[0][j])  # cross correlation
                    elif i == j:
                        sol_anum += np.dot(cut_sol[k][i], cut_sol[0][j])  # auto correlation
            sol_cross.append(sol_cnum)
            sol_auto.append(sol_anum)
        sleep(0.05)

        # ==========水分子と溶質分子の相互相関==========
        # 分母(t = 0)
        mix_czero = 0
        for i in tqdm(range(int(s)), desc=f"計算準備中({ope}aq)"):
            for j in range(int(w)):
                mix_czero += np.dot(cut_sol[0][i], wat[0][j])  # Static auto correlation
        sleep(0.05)

        # 分子(t = t')
        mix_cross = []
        for k in tqdm(range(int(cut)), desc=f"計算中({ope}aq_A&C)"):
            mix_cnum = 0
            for i in range(int(s)):
                for j in range(int(w)):
                    mix_cnum += np.dot(cut_sol[k][i], wat[0][j])  # cross correlation
            mix_cross.append(mix_cnum)
        sleep(0.05)

        # ==========初期位置の双極子モーメントの自己・相互相関関数で正規化==========
        # ==水==
        wat_num = np.array(wat_auto) / w + 2 * np.array(wat_cross) / wcdom
        wat_dom = np.array(wat_azero) / w + 2 * np.array(wat_czero) / wcdom
        water = wat_num / wat_dom
        # ==溶質==
        sol_num = np.array(sol_auto) / s + 2 * np.array(sol_cross) / scdom
        sol_dom = np.array(sol_azero) / s + 2 * np.array(sol_czero) / scdom
        substance = sol_num / sol_dom
        # ==水-溶質==
        mix_num = np.array(mix_cross) / (w * s)
        mix_dom = np.array(mix_czero) / (w * s)

        mix = (wat_num + sol_num + mix_num) / (wat_dom + sol_dom + mix_dom)
        solute = sol_num / sol_dom

        func = np.c_[np.array(time_axis), np.array(mix)]
        func_wat = np.c_[np.array(time_axis), np.array(water)]
        func_sol = np.c_[np.array(time_axis), np.array(substance)]

        file_path = os.path.dirname(top)
        print("計算完了")

        # 重量パーセント濃度
        wt = (Mw * s) / (Mw * s + WAT * w) * 100

        # 水溶液の水A&C
        name = os.path.join(file_path, f"{ope}aq_solvent.csv")
        np.savetxt(
            name,
            func_wat,
            delimiter=","
        )

        # 水溶液の溶質A&C
        name = os.path.join(file_path, f"{ope}aq_substansave.csv")
        np.savetxt(
            name,
            func_sol,
            delimiter=","
        )

        # 水溶液のA&C
        name = os.path.join(file_path, f"{ope}{wt:.2f}wt%.csv")
        np.savetxt(
            name,
            func,
            delimiter=","
        )

        print("保存が完了しました。")
        print(f"ファイルパス: {name}")

        # モル比率濃度
        mol = s / (s + w) * 100

        print("=" * 100)
        print("数値を使用する場合はメモをする")
        print(f"{ope}の重量濃度は{wt:.2f}wt%です")
        print(f"{ope}のモル比率は{mol:.2f}mol%です")


# =====本計算============================================================================================================
print("pdbファイルかcrdファイルを選択してください")
TOP = top()  # pdb, crdファイルパスをよみこむ
TOP_path = os.path.splitext(TOP)  # TOPのパスTOP_path[0]と拡張子TOP_oath[1]を分割する
top_dom = TOP_path[1]  # TOPファイルの拡張子
ABS_top = os.path.abspath(TOP)  # TOP絶対パス

# =====DS上で作成したデータを用いた計算=================================================================================
if top_dom == ".pdb":
    top = TOP
    print(f"pdb file:[{ABS_top}]")
    print("dcdファイルを選んでください")
    DCD = dcd()
    ABS_dcd = os.path.abspath(DCD)
    print(f"dcd file:[{ABS_dcd}]", end="\n\n")
    u = md.Universe(TOP, DCD)
    print(f"原子名/原子グループ：[{u.atoms}]", end="\n\n")

    choose = input("純物質？？：　")
    if choose == "":
        meta = DS(u)

        ac(top, ope=meta[0], mate=meta[1])
    else:
        wat = WAT_dip(u)
        meta = DS(u)
        acaq(top, ope=meta[0], wat=wat, sub=meta[1])

    # =====計算サーバー上で作成したデータを用いた計算=========================================================================
elif top_dom == ".crd":
    top = TOP
    print(f"crd file:[{ABS_top}]")
    print("dcdファイルを選んでください")
    DCD = dcd()
    ABS_dcd = os.path.abspath(DCD)
    print(f"dcd file:[{ABS_dcd}]", end="\n\n")
    u = md.Universe(top, DCD)
    print(f"原子名/原子グループ：[{u.atoms}]", end="\n\n")

    choose = input("純物質？？：　")
    if choose == "":
        meta = server(u)
        ac(top, ope=meta[0], mate=meta[1])
    else:
        wat = WAT_dipS(u)
        meta = server(u)
        acaq(top, ope=meta[0], wat=wat, sub=meta[1])
