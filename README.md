# LAMMPS-utils
python scripts for LAMMPS

# poscar_lammps.pyの使い方

Setupの項目に設定を入力する：

  filename       # VASP POSCARのファイル名を入力
  ntype          # 原子種の数を入力
  nbtype         # 結合の種類の数を入力
  natype         # 結合角の種類の数を入力
  ndtype         # 二面角の種類の数を入力
  nitype         # improper角の種類の数を入力

  astyle         # atomic, charge, molecular, または fullを入力 
        
  q[X]           # 原子タイプXの電荷を入力（タイプ番号Xは0始まりインデックス，以下同様）
  mass[X]        # 原子タイプXの質量を入力

  btype[X]       # 結合タイプXを構成する原子種ペアを入力 ※nbtype=0の場合はコメントアウト 

  atype[X]       # 結合角タイプXを構成する原子種ペアを入力 ※natype=0の場合はコメントアウト

  dtype[X]       # 二面角タイプXを構成する原子種ペアを入力 ※ndtype=0の場合はコメントアウト

  itype[X]       # improper角タイプXを構成する原子種ペアを入力 ※nitype=0の場合はコメントアウト

  mul            # ユニットセルをmul[0]*mul[1]*mul[2]倍したスーパーセルを作成

  bcoeff         # Bond Coeffsに記入する文字列を入力
  acoeff         # Angle Coeffsに記入する文字列を入力
  dcoeff         # Dihedral Coeffsに記入する文字列を入力
  icoeff         # Improper Coeffsに記入する文字列を入力

  cutoff[X][X]   # 原子種X-原子種X間の結合カットオフ距離(Å)を入力

  ある原子種のペアに対して，原子間距離が結合カットオフ距離以下である場合に結合していると判定．
  結合カットオフ距離を0に設定した場合，その原子種のペアの間には結合は生じないと見なす．
  本プログラムでは，この結合カットオフ距離に基づいてNeighbor listを作り，これによって 
  Molecule-ID, Bonds, Angles, Dihedrals, Impropers のリストを作成する．
