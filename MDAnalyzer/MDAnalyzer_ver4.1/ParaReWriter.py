import tkinter as tk
from tkinter import filedialog as fd
import sys
import re
import os
import shutil
import numpy as np


# ==========実行内容==========


def main():
    Folder_Path = Find_Folder()
    print(f"選択されてたフォルダパス : {Folder_Path}")

    Para_Folders = (["1_Minimization", "minimization.inp"],
                    ["2_Minimization", "minimization.inp"],
                    ["3_Heating", "heating.inp"],
                    ["4_Equilibration", "equilibration.inp"],
                    ["5_Production", "production.inp"]
                    )

    select = input("ACCELRYS_LICENSE_USERNAMEを書き換えますか(yes or no) : ")
    if select == "yes":
        my_num = input("自分の学籍番号 : ")
        print("", end="\n\n")
        for name in Para_Folders:
            cstudent_number(Folder_Path, Folder_name=name[0], New_USERNAME=my_num)
    else:
        pass

    # =====inpファイルのパスの作成====
    Paths = Parameter_File_Paths(Folder_Path, Para_Folders)

    # =====Minimization, Minimization2=====
    for counter1 in range(2):
        print(f"[{Para_Folders[counter1][0]}>>>{Para_Folders[counter1][1]}]")
        print(Paths[counter1])
        Minimization(Paths[counter1], counter1)
        print("", end="\n")

    # =====Heating=====
    print(f"[{Para_Folders[2][0]}>>>{Para_Folders[2][1]}]")
    print(Paths[2])
    FINAL = Heating(Paths[2])
    print("", end="\n")

    # =====Equilibration, Production=====
    for counter2 in range(2):
        number = int(counter2 + 3)
        print(f"[{Para_Folders[number][0]}>>>{Para_Folders[number][1]}]")
        print(Paths[number])
        the_others(Paths[number], FINAL)
        print("", end="\n")

    # =====内容確認=====
    print("=====<<Parameters>>=====")
    i_number = 0
    for counter in Paths:
        name = os.path.basename(counter)
        Step_Name = Para_Folders[i_number][0]
        with open(counter) as f:
            para = f.readlines()

        if name == "minimization.inp":
            step = Step_Name
            print(f"[{step}]")

            Para_line = para[40].split()
            X = Para_line[3]
            print(f"MaxStep : {X}", end="\n\n")

        elif name == "heating.inp":
            step = Step_Name
            print(f"[{step}]")
            Time_line = para[40].split()
            Time_Step = float(Time_line[1])
            Now_NSTEP = float(Time_line[3])
            Now_Simu_Time = int(Now_NSTEP * Time_Step)
            print(f"Simulation Time : {Now_Simu_Time}ps")
            Temp_line = para[41].split()
            rFIRST = Temp_line[1]
            rFINAL = Temp_line[3]
            print(f"Initial Temperature : {rFIRST}K")
            print(f"Target Temperature : {rFINAL}K", end="\n\n")


        elif name == "equilibration.inp":
            step = Step_Name
            print(f"[{step}]")
            Para_line = para[40].split()
            Time_Step = float(Para_line[1])
            Now_NSTEP = float(Para_line[3])
            Now_Simu_Time = int(Now_NSTEP * Time_Step)
            print(f"Simulation Time : {Now_Simu_Time}ps")

            Temp_line = para[41].split()
            rFIRST = Temp_line[1]
            rFINAL = Temp_line[3]
            print(f"Initial Temperature : {rFIRST}K")
            print(f"Target Temperature : {rFINAL}K", end="\n\n")

        elif name == "production.inp":
            step = Step_Name
            print(f"[{step}]")
            Para_line = para[40].split()
            Time_Step = float(Para_line[1])
            Now_NSTEP = float(Para_line[3])
            Now_Simu_Time = int(Now_NSTEP * Time_Step)
            print(f"Simulation Time : {Now_Simu_Time}ps")

            Temp_line = para[41].split()
            rFIRST = Temp_line[1]
            rFINAL = Temp_line[3]
            print(f"Initial Temperature : {rFIRST}K")
            print(f"Target Temperature : {rFINAL}K", end="\n\n")

        i_number += 1


# ==========フォルダパスを取得する==========


def Find_Folder():
    root = tk.Tk()
    root.withdraw()

    fpath = fd.askdirectory(
        initialdir=sys.argv[0]
    )
    root.quit()
    return fpath


# ==========CHARMm.shファイルのACCELRYS_LICENSE_USERNAMEを書き換える。==========

def cstudent_number(MDResult_file_path, Folder_name, New_USERNAME):
    Intermediate_Folder_Path = os.path.join(MDResult_file_path, "Intermediate")
    File_path = os.path.join(Intermediate_Folder_Path, f"{Folder_name}")
    CHARMM_path = os.path.join(File_path, "RunCHARMm.sh")
    print(f'==========[{Folder_name}>>>RunCHARMm.sh]==========')
    print(f"Path :{CHARMM_path}")

    # =====変更前のファイル(バックアップ)を作成=====
    Back_Up_Name = CHARMM_path + ".bak"
    shutil.copy(CHARMM_path, Back_Up_Name)

    # ==========ファイルを読み込む==========
    with open(CHARMM_path) as f:
        text = f.read()

    USERNAME_line = re.findall('ACCELRYS_LICENSE_USERNAME.*"', text)
    # print(USERNAME_line)
    USERNAME = re.findall('"(.*)"', f'{USERNAME_line[0]}')
    # print(USERNAME)

    # =====文字列の置換=====
    text = text.replace(f'{USERNAME[0]}', f'{New_USERNAME}')
    print(text)
    print("書き換え完了", end="\n\n")

    # =====変更内容を保存=====
    with open(CHARMM_path, mode="w") as f:
        f.writelines(text)

    return


# ==========フォルダ内の~.inpファイルを探す==========
def Parameter_File_Paths(MDResult_file_path, Para_Folders):
    Intermediate_Folder_Path = os.path.join(MDResult_file_path, "Intermediate")
    # print(Intermediate_Folder_Path) #Intermediateフォルダのカレントディレクトの取得

    Para_Files = []

    # 1_Minimization, ...5_Productionのパスを作成
    for i in Para_Folders:
        Folder_Path = os.path.join(Intermediate_Folder_Path, f"{i[0]}")
        Para_Path = os.path.join(Folder_Path, f"{i[1]}")

        # 1_Minimization, ...5_Productionの有無をboolで処理
        Path = os.path.exists(Para_Path)

        if Path:
            Para_Files.append(Para_Path)
            print(f"{i[0]}>>>{i[1]} : Found")
        else:
            print("ファイルが見つかりませんでした。")

    return Para_Files


def Minimization(Para_Path, Folder_counter):
    # ==========ファイルを読み込む==========
    with open(Para_Path) as f:
        texts = f.readlines()
    Para_line = texts[40].split()
    Dfo_Max_Step = str(Para_line[3])

    print("MaxStep数を入力してください(Max : 1e9)")
    print("変更の無い場合はエンターキーを押してください。")
    print(f"現在のMaxStep : {Dfo_Max_Step}")
    Max_Step = input("変更後のMaxStep : ")

    if Max_Step == "":
        # 操作をスキップ
        pass

    elif int(Max_Step) <= (Folder_counter + 1) * 1000000000:

        # =====変更前のファイル(バックアップ)を保存=====
        Back_Up_Name = Para_Path + ".bak"
        shutil.copy(Para_Path, Back_Up_Name)

        # =====NSTEPパラメータを変更=====
        Dfo_Max_Step = str(Max_Step)
        texts[40] = texts[40] = "  " + ' '.join(Para_line) + " \n"
        # print(texts[40])

        # =====変更内容を保存=====
        with open(Para_Path, mode="w") as f:
            f.writelines(texts)

    else:
        print("値が大きすぎます")
        print("最初からやり直してください")


def Heating(Para_Path):
    # ==========ファイルを読み込む==========
    with open(Para_Path) as f:
        texts = f.readlines()
    Para_line = texts[40].split()
    Time_Step = float(Para_line[1])
    Now_NSTEP = Para_line[3]
    Now_Simu_Time = int(float(Now_NSTEP) * Time_Step)
    print("Simulation Timeを入力してください(推奨 : 4ps)")
    print(f"現在のSimulation Time : {Now_Simu_Time}[ps]")
    Time = input("変更後のSimulation Time[ps] : ")

    # Simulation Time　の上限判定 and 数字かどうか判定
    if int(Time) <= 1000000000:
        FIRST = input("Initial_Temperature[K] : ")
        FINAL = input("Target_Temperature[K] : ")
        # Initial_Temperature = int(Target_Temperature) + 100
        # =====変更前のファイル(バックアップ)を保存=====
        Back_Up_Name = Para_Path + ".bak"
        shutil.copy(Para_Path, Back_Up_Name)

        # =====NSTEPパラメータを変更=====
        # print(texts[41])
        MaxStep = round(int(Time) // float(Time_Step), -1)

        Para_line[3] = str(MaxStep)
        # print(Para_line[3])
        texts[40] = "  " + ' '.join(Para_line) + " \n"
        # print(texts[40])

        # ==========改装中=========
        Temp_Para_line = texts[41].split()
        # print(Temp_Para_line)

        # =====Initial Temp=====
        Temp_Para_line[1] = FIRST
        # print(Temp_Para_line[1])

        # =====Target Temp=====
        Temp_Para_line[3] = FINAL
        # print(Temp_Para_line[3])

        temp = "  " + ' '.join(Temp_Para_line) + " \n"
        # print(temp)
        texts[41] = temp
        print(np.array(texts[40 : 42]))

        # ==========================


        # =====変更内容を保存=====
        with open(Para_Path, mode="w") as f:
            f.writelines(texts)

    else:
        print("値が大きすぎます")
        print("最初からやり直してください")

    return FINAL


# =====Equilibration, Production=====
def the_others(Para_Path, Target_Temperature):
    # ==========ファイルを読み込む==========
    with open(Para_Path) as f:
        texts = f.readlines()
    Para_line = texts[40].split()
    # print(Para_line)
    Time_Step = float(Para_line[1])
    Now_NSTEP = float(Para_line[3])
    Now_Simu_Time = round(Time_Step * Now_NSTEP, -1)
    print("Simulation Timeを入力してください")
    print(f"現在のSimulation Time : {Now_Simu_Time}[ps]")
    Time = input("変更後のSimulation Time[ps] : ")

    # Simulation Time　の上限
    if int(Time) <= 1000000000:

        # =====変更前のファイル(バックアップ)を保存=====
        Back_Up_Name = Para_Path + ".bak"
        shutil.copy(Para_Path, Back_Up_Name)
        # =====NSTEPパラメータを変更=====
        MaxStep = round(int(Time) // Time_Step, -1)
        Para_line[3] = str(MaxStep)
        print(Para_line[3])
        replased_line = "  " + ' '.join(Para_line) + " \n"
        # print(replased_line)
        texts[40] = replased_line
        # =====Temperature=====
        Temp_Para_line = texts[41].split()
        FIRST = texts[41].replace(f"{Temp_Para_line[1]}", f"{Target_Temperature}")
        # FINAL = FIRST.replace(f"{Temp_Para_line[3]}", f"{Target_Temperature}")
        texts[41] = FIRST

        print(np.array(texts[40:42]))

        # =====変更内容を保存=====
        with open(Para_Path, mode="w") as f:
            f.writelines(texts)

    else:
        print("値が大きすぎます")
        print("最初からやり直してください")


main()
