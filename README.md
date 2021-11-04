# LAMMPS-utils
python scripts for LAMMPS

poscar_lammps.pyの使い方：

## Setup ##
  filename = 'POSCAR.vasp'    # VASP POSCARのファイル名を入力
  ntype  = 2                  # 原子種の数を入力
  nbtype = 1                  # 結合の種類の数を入力
  natype = 1                  # 結合角の種類の数を入力
  ndtype = 0                  # 二面角の種類の数を入力
  nitype = 0                  # improper角の種類の数を入力

  astyle = 'full'             # atomic, charge, molecular, または fullを入力 

  q = [0 for i in range(ntype)]; mass = [0 for i in range(ntype)]           
  q[0]     = -0.82000         # 原子種0の電荷を入力（only for astyle = 'charge' or 'full'）
  q[1]     = 0.41000          # 原子種1の電荷を入力（only for astyle = 'charge' or 'full'）
  mass[0]  = 16.000           # 原子種0の質量を入力
  mass[1]  = 1.010            # 原子種1の質量を入力

  btype = [0 for i in range(nbtype)]
  btype[0] = [0, 1]           # 結合0を構成する原子種ペアを入力（番号は0始まり）        

  atype = [0 for i in range(natype)]
  atype[0] = [1, 0, 1]        # 結合角0を構成する原子種ペアを入力（番号は0始まり）

  dtype = [0 for i in range(ndtype)]
#  dtype[0] = [0, 1, 1, 1]     # 二面角0を構成する原子種ペアを入力（番号は0始まり）

  itype = [0 for i in range(nitype)]
#  itype[0] = [1, 1, 1, 1]     # improper角0を構成する原子種ペアを入力（番号は0始まり）

  mul = [2, 2, 2]             # ユニットセルをmul[0]*mul[1]*mul[2]倍したスーパーセルを作成

  bcoeff = ' 1 22.965000 1.0120000'          # Bond Coeffsに記入する文字列
  acoeff = ' 1 harmonic 1.6456800 113.24000' # Angle Coeffsに記入する文字列
  dcoeff = ''                                # Dihedral Coeffsに記入する文字列
  icoeff = ''                                # Improper Coeffsに記入する文字列

  cutoff = [[0 for i in range(ntype)] for j in range (ntype)]
  cutoff[0][0] = 0.0          # 原子種0-原子種0間の結合カットオフ距離(Å)を入力（番号は0始まり）
  cutoff[0][1] = 1.2          # 原子種0-原子種1間の結合カットオフ距離(Å)を入力（番号は0始まり）
  cutoff[1][1] = 0.0          # 原子種1-原子種1間の結合カットオフ距離(Å)を入力（番号は0始まり）

  # Note ##########################################################################
  # ある原子種のペアに対して，原子間距離が結合カットオフ距離以下である場合に結合していると判定．
  # 結合カットオフ距離を0に設定した場合，その原子種のペアの間には結合は生じないと見なす．
  # 本プログラムでは，この結合カットオフ距離に基づいてNeighbor listを作り，これによって 
  # Molecule-ID, Bonds, Angles, Dihedrals, Impropers のリストを作成する．
  #################################################################################
