from tqdm import tqdm
import subprocess
import numpy as np
import os

import folder_processer as folder
import file_processer as file
import Path_pull as PP

gver1 = "GHz"
GHz = 10**9

gvar2 = "ps"
ps = 10 ** -12

#==========モデル関数生成==========
def model():
    print("保存先のフォルダを選択してください")
    main_path = folder.finder()

    print("モデル関数を格納するフォルダ名を入力してください")
    fname = input("フォルダ名: ")
    save_path = folder.maker(main_path, fname)

    end_time = float(input("計算時間[ps]: "))
    tau = float(input("緩和時間[ps]: "))
    dc = float(input("DC成分: "))
    dt = float(input("時間刻み[ps]: "))

    data = []
    i = 0
    time = float(0)
    while time < end_time:
        time = i * dt
        fu = np.exp(- time / tau) + dc
        data.append([time * ps, fu])
        i += 1

    fname = os.path.join(save_path, f"{fname}.csv")
    np.savetxt(f"{fname}",
               data,
               delimiter=",")

    print("")
    print("モデル関数が生成されました")
#==========ラプラス・フーリエ変換単独実行==========

def single_LFT():
    file_path = PP.csv()
    folder_path = os.path.dirname(file_path)
    save_name = input("保存名: ")
    relax_path = folder.maker(folder_path, save_name)



    s_fre = float(input("開始周波数[GHz]: "))
    e_fre = float(input("終端周波数[GHz]: "))
    df = float(input("周波数刻み[GHz]: "))

    print("=====Start　LFT=====")

    i = 0
    relax_data = []
    fre = s_fre
    csv = file.csv_reader(file_path)
    print(f"目標周波数{e_fre}")
    #名前　計算実行周波数
    with tqdm() as pber:
        while fre <= e_fre:
            fre = s_fre + i * df
            w = 2 * np.pi * fre * GHz

            re = np.sum(csv[:, 1] * np.cos(w * csv[:, 0]) * csv[1][0])
            im = np.sum(csv[:, 1] * np.sin(w * csv[:, 0]) * csv[1][0])

            relax_data.append([np.log10(fre * GHz), re, im])
            pber.update(df)
            i += 1
    relax = os.path.join(relax_path, f"{save_name}.txt")
    np.savetxt(relax, relax_data, delimiter="\t")


# メイン部分
while True:
    print("= " * 40)
    print("モデル関数作成: 1")
    print("auto-and cross-correlation function 生成: 2")
    print("ラプラス・フーリエ変換: 3")
    print("精度分析: 4")
    print("終了コマンド: Push Enter_key")
    ope_1 = input("操作番号を入力してください。(半角英数): ")
    print("")

    if (ope_1 == "1"):
        model()
    elif(ope_1 == "2"):
        print("=====溶液系を選んでください=====")
        print("Pure: 0")
        print("Mixture: 1")
        ope_2 = input("system: ")
        if(ope_2 == "0"):
            subprocess.run(["Python", "Pure_dipole.py"], check=True)
        elif(ope_2 == "1"):
            subprocess.run(["Python", "Mix_dipole.py"], check=True)
    elif(ope_1 == "3"):
        single_LFT()
    elif(ope_1 == "4"):
        subprocess.run(["Python", "Result_Analyzer.py"], check=True)
    elif ope_1 == "":  # エンターキーを押した場合、MD_relaxation_maker.pyを終了する.
        print("<<<操作を終了します>>>")
        break
    else:  # その他の入力は再入力をようきゅされるようにwhileループさせている.
        print("コマンドが間違っています。最初からやり直してください\n")