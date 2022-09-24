import tkinter as tk
from tkinter import filedialog as fd
import sys
import os
import shutil
import numpy as np
import natsort as na
import glob
from tqdm import tqdm
from time import sleep
import math

# =====define unit=====
ps = 10 ** -12


# ==========フォルダパスを取得する================================================================


def finder():
    root = tk.Tk()
    root.withdraw()

    folder_path = fd.askdirectory(
        initialdir=sys.argv[0]
    )
    root.quit()
    return folder_path


# =====call file path(.csv file)=====


def csv():
    root = tk.Tk()
    root.withdraw()
    filename = [("CSVファイル", "*.csv")]

    csv_file = fd.askopenfilename(
        filetype=filename,
        initialdir=sys.argv[0]
    )
    root.quit()
    return csv_file


# =====read file path=====


def csv_reader():
    # =====jadged file path=====
    while True:
        csv_path = csv()
        jad = bool(csv_path)
        if jad:
            break
        else:
            print('<<<ERROR>>>')
            print("ファイルを選択してください。")

    csv_data = np.loadtxt(
        fname=csv_path,
        dtype=float,  # float型として読み込む
        delimiter=","  # コンマ区切りのデータとして読み込む
    )
    return csv_path, csv_data

# =====load file path=====


def csv_loader(csv_path):

    # ===== load csv data
    csv_data = np.loadtxt(
        fname=csv_path,
        dtype=float,  # float型として読み込む
        delimiter=","  # コンマ区切りのデータとして読み込む
    )
    return csv_data


# =====make folder with csv directory=====


def maker(folder_dir, name):
    # 指定したディレクトリ下にnameのフォルダ名の入ったパスを生成する
    directory = os.path.join(folder_dir, f"{name}")
    # ファルダ生成を行う
    try:
        os.makedirs(directory)
    except FileExistsError:
        print("<<<Caution>>>")
        print(f"{directory}は既に存在しています。")
        sle = str(input("消去しますか(yes or no): "))
        while True:
            if sle == "yes":
                shutil.rmtree(directory)
                os.makedirs(directory)
                break
            elif sle == "no":
                print("別名ファイルを作成します")
                new_name = input("ファイル名: ")
                directory = os.path.join(folder_dir, f"{new_name}")
                os.makedirs(directory)
                break
            else:
                print("非対応文字です。もう一度やり直して", end="\n\n")

    return directory


# ===== get all csv file path in the folder=====


def lister(folder_dir):
    files = folder_dir + "/*.csv"
    get_paths = [glob.glob(files)]
    file_path = na.natsorted(get_paths[0])
    return file_path


# ===== dividing time range=====


def dividing_time(csv_path, csv_data):
    csv_dir = os.path.dirname(csv_path)
    dt = csv_data[1][0]
    time = int((int(len(csv_data)) + 1) * dt)
    print(f"全時間: {time}ps")
    print(f"下限時間刻み: {dt}ps")
    print("")

    while True:
        try:
            all_time = float(input('時間範囲[ps]: '))
            time_range = float(input('時間窓の刻み[ps]: '))
            break
        except ValueError:
            print("非対応文字が含まれています")

    koma = int(all_time // time_range)

    folder_path = maker(csv_dir, "time_windows")

    # 指定範囲まで別ファイルとして保存する

    for div_time in tqdm(range(koma), desc="分割中"):
        last_line = int((div_time + 1) * time_range // dt) + 2
        div_csv = csv_data[:last_line, :]
        wt = "{:.2f}".format((div_time + 1) * time_range)
        file_name = os.path.join(folder_path, f"time_range0.00ps~{wt}ps.csv")

        np.savetxt(
            file_name,
            div_csv,
            delimiter=","
        )
    sleep(0.05)
    return folder_path, time_range


# =====Laplace fourier translation =====


def continue_LFT(window_path, time_range):
    LFT_path = maker(window_path, name="MD_relaxation")
    files = lister(window_path)

    while True:
        try:
            start_fre = float(input("開始周波数[Hz]: "))
            end_fre = float(input("終端周波数[Hz]: "))
            df = float(input("周波数刻み[Hz]: "))
            break
        except ValueError:
            print("非対応文字が含まれています",end="\n\nsss")

    # =====指定周波数帯域を生成する=====
    frequency = np.arange(start_fre, end_fre, df)

    # =====指定周波数帯域を生成する(等間隔)=====
    # frequency = []
    # fre_range = int(math.log10(end_fre // start_fre))
    # for i in range(fre_range):
    #    area = np.arange(start_fre * 10 ** i, start_fre * 10 ** (i + 1), df * 10 ** i)
    #   frequency.extend(area)
    # frequency.extend([end_fre])
    # sleep(0.05)

    # 工学屋さんの対数で加工する(物理屋さんはln)
    frequency10 = np.log10(np.array(frequency))
    print("=====Start LFT=====")

    dloss_peack = []
    relax_time = []
    timerange_axis = []
    count = 1

    sleep(0.05)
    for number in tqdm(files, desc="ラプラス・フーリエ変換中"):
        csv = csv_loader(number)
        all_time = int(len(csv))
        dt = csv[1][0]

        # =====前進差分法(中心差分法がよかったな・・・)=====

        dfunc = np.diff(csv[:, 1])
        # 自己・相互相関関数の微分関数
        function = -1 * np.array(dfunc) / dt
        # 時間軸の要素数合わせのために時間軸要素1減らす
        time_axis = csv[:all_time - 1, 0]

        Re = []
        Im = []
        # =====start Laplace fourier translation===
        for fre in frequency:
            w = 2 * math.pi * fre

            # =====DFT=====
            real = np.sum(np.array(function) * np.cos(w * np.array(time_axis) * ps) * csv[1][0] * ps)
            image = np.sum(np.array(function) * np.sin(w * np.array(time_axis) * ps) * csv[1][0] * ps)

            # =====FFT=====

            # =====開始周波数の整数倍の結果のみを保存する=====
            # if fre == stock_number * start_fre:
            Re.append(real)
            Im.append(image)

        # ===== dielectric relaxation curves =====
        es = Re[0]
        wt = (count * time_range)
        Real_norm = np.array(Re) / es
        Imag_norm = np.array(Im) / es
        relaxation_curves = np.c_[frequency10, Real_norm, Imag_norm]
        fpath = os.path.join(LFT_path, f"relax{wt}ps.txt")
        np.savetxt(
            fpath,
            relaxation_curves,
            delimiter="\t"
        )

        # ===== peack of the dielectoric loss =====
        imax = np.max(Imag_norm)
        # ===== relaxation times =====
        # =====誘電損失がピークになる配列数を探す=====
        imax_number = np.argmax(Im)
        # =====誘電損失ピークに対応する周波数を探す=====
        peack_fre = frequency[imax_number]
        # =====緩和時間を計算=====
        rtime = 1 / (2 * math.pi * peack_fre) / ps

        dloss_peack.append(imax)
        relax_time.append(rtime)
        timerange_axis.append(wt)
        count += 1
    sleep(0.05)

    relax_path = maker(LFT_path, name="analyze_data")

    dloss_data = np.c_[timerange_axis, dloss_peack]
    relax_filename = os.path.join(relax_path, "Dielectric_loss.csv")
    np.savetxt(
        relax_filename,
        dloss_data,
        delimiter=","
    )
    relax_time_data = np.c_[timerange_axis, relax_time]
    relax_filename = os.path.join(relax_path, "Relaxation_time.csv")
    np.savetxt(
        relax_filename,
        relax_time_data,
        delimiter=","
    )


# =====handle=====

csv_info = csv_reader()
object = dividing_time(csv_path=csv_info[0], csv_data=csv_info[1])
continue_LFT(window_path=object[0], time_range=object[1])
