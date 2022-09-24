import Pure_dipole as PD
import Path_pull as PP

print("PDBファイルを選択してください")
PDB = PP.pdb()
print(PDB)
print("DCDファイルを選択してください")
DCD = PP.dcd()
print(DCD)
ts = float(input("シミュレーション時間刻み[ps]: "))
ps = 10**-12

wat = PD.wat_dip(PDB, DCD)

