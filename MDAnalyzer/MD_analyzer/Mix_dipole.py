import os
import sys
import numpy as np
from tqdm import tqdm
import Path_pull as PP
import Pure_dipole as Pd

ps = 10**-12

if __name__=="__main__":
    print("溶質名を指定してください")
    print("<<<溶質: 操作番号>>>")
    print("methanol: 0")
    print("1,2-Dichloroethane: 1")
    print("1,1,2-Trichloroethane: 2")
    print("urea: 3")
    print("glycine: 4")
    print("β-alanine: 5")
    ope = input("操作番号: ")

    print("PDBファイルを選択してください")
    PDB = PP.pdb()
    print(f"選択ファイル{PDB}")
    print("DCDファイルを選択してください")
    DCD = PP.dcd()
    print(f"選択ファイル{DCD}")
    ts = float(input("シミュレーション時間刻み[ps]: "))

    wat = Pd.wat_dip(PDB, DCD)
    t = len(wat)
    w = len(wat[0])
    wcdom = w ** 2 - w
    WAT = 18.0152

    if ope == str(0):
        sub = "Methanol"
        Mw = 16.04240
        sol = Pd.met_dip(PDB, DCD)
    elif ope == str(1):
        sub = "1,2-Dichloroethane"
        Mw = 98.9585
        sol = Pd.dce12_dip(PDB, DCD)
    elif ope == str(2):
        sub = "1,1,2-Trichloroethane"
        Mw = 133.4033
        sol = Pd.tce112_dip(PDB, DCD)
    elif ope == str(3):
        sub = "Urea"
        Mw = 60.055
        sol = Pd.urea_dip(PDB, DCD)
    elif ope == str(4):
        sub = "Glycine"
        Mw = 75.0666
        sol = Pd.gly_dip(PDB, DCD)
    elif ope == str(5):
        sub = "β-Alanine"
        Mw = 89.0932
        sol = Pd.b_ala_dip(PDB, DCD)

    # 混合液系の計算===========================================================
    # ==========水の自己・相互相関==========
    w = len(wat[0])
    wcdom = w ** 2 - w

    # 分母(t = 0)
    wat_czero = 0
    wat_azero = 0
    for i in tqdm(range(int(w)), desc=f"計算準備中(water)"):
        for j in range(int(w)):
            if i < j:
                wat_czero += np.dot(wat[0][i], wat[0][j])  # Static auto correlation

            elif i == j:
                wat_azero += np.dot(wat[0][i], wat[0][j])  # Static cross correlation

    # 分子(t = t')
    time = []
    wat_cross = []
    wat_auto = []
    for k in tqdm(range(int(t)), desc="計算中(water_A&C)"):
        wat_cnum = 0
        wat_anum = 0
        for i in range(int(w)):
            for j in range(int(w)):
                if i < j:
                    wat_cnum += np.dot(wat[k][i], wat[0][j])  # cross correlation

                elif i == j:
                    wat_anum += np.dot(wat[k][i], wat[0][j])  # auto correlation

        time.append(k * ts * ps)
        wat_cross.append(wat_cnum)
        wat_auto.append(wat_anum)

    # ==========溶質の自己・相互相関==========
    # 分母(t = 0)
    s = len(sol[0])
    scdom = s ** 2 - s

    sol_czero = 0
    sol_azero = 0
    for i in tqdm(range(int(s)), desc=f"計算準備中({sub})"):
        for j in range(int(s)):
            if i < j:
                sol_czero += np.dot(sol[0][i], sol[0][j])  # Static auto correlation

            elif i == j:
                sol_azero += np.dot(sol[0][i], sol[0][j])  # Static cross correlation

    # 分子(t = t')
    sol_cross = []
    sol_auto = []
    for k in tqdm(range(int(t)), desc=f"計算中({sub}_A&C)"):
        sol_cnum = 0
        sol_anum = 0
        for i in range(int(s)):
            for j in range(int(s)):
                if i < j:
                    sol_cnum += np.dot(sol[k][i], sol[0][j])  # cross correlation

                elif i == j:
                    sol_anum += np.dot(sol[k][i], sol[0][j])  # auto correlation

        sol_cross.append(sol_cnum)
        sol_auto.append(sol_anum)

    # ==========水と溶質の相互相関==========
    # 分母(t = 0)
    mix_czero = 0
    for i in tqdm(range(int(s)), desc=f"計算準備中({sub}aq)"):
        for j in range(int(w)):
            mix_czero += np.dot(sol[0][i], wat[0][j])  # Static auto correlation

    # 分子(t = t')
    mix_cross = []
    for k in tqdm(range(int(t)), desc=f"計算中({sub}aq_A&C)"):
        mix_cnum = 0
        for i in range(int(s)):
            for j in range(int(w)):
                mix_cnum += np.dot(sol[k][i], wat[0][j])  # cross correlation

        mix_cross.append(mix_cnum)

    # ==========初期位置の双極子モーメントの自己・相互相関関数で正規化==========
    # ==水==
    wat_num = np.array(wat_auto) / w + 2 * np.array(wat_cross) / wcdom
    wat_dom = np.array(wat_azero) / w + 2 * np.array(wat_czero) / wcdom
    # ==溶質==
    sol_num = np.array(sol_auto) / s + 2 * np.array(sol_cross) / scdom
    sol_dom = np.array(sol_azero) / s + 2 * np.array(sol_czero) / scdom
    # ==水-溶質==
    mix_num = 2 * np.array(mix_cross) / (w * s)
    mix_dom = 2 * np.array(mix_czero) / (w * s)

    mix = (wat_num + sol_num + mix_num) / (wat_dom + sol_dom + mix_dom)

    func = np.c_[np.array(time), np.array(mix)]
    fpath = os.path.dirname(PDB)

    print("計算完了")

    wt = (Mw * s) / (Mw * s + WAT * w) * 100
    name = os.path.join(fpath, f"{sub}{wt:.2f}wt%.csv")
    np.savetxt(
        name,
        func,
        delimiter=","
    )

    print("保存が完了しました。")
    print(f"ファイルパス: {name}")
    mol = s / (s + w) * 100
    print(f"{sub}の重量濃度は{wt:.2f}wt%です")
    print(f"{sub}のモル比率は{mol:.2f}%です")