# improved_rDock
分子シミュレーションソフトウェアrDockに、疾患標的タンパク質の機能を阻害/活性化する低分子化合物やその結合様式を計算により高速・高精度に予測するための機能を追加した改良版

## 動作環境
- AmberTools, Open Babelがインストールされている必要がある。
- 本スクリプトで使用する下記実行形式へのパスが通っている必要がある。
  - MMPBSA.py
  - pymdpbsa
  - tleap
  - pytleap
  - antechamber
  - parmchk
  - babe
- /~/amberXX/binの中にあるtleapの内容の一部を下記のように変更して, 名前をtleap_rDockとする （AmberTools17の場合、tleapの内容を下記のように修正する必要がある）。
  - 変更前 `export AMBERHOME="$(dirname "$(cd "$(dirname "$0")" && pwd)")"`
  - 変更後 `export AMBERHOME=/~/amberXX`
- 以下のPythonスクリプトを作業ディレクトリ内に配置する必要がある.
  - to_after_Rdock_argv_5_br171114.py
  - conv.py
## 改良版rDockの機能

### 側鎖構造最適化
rbdockのコマンドオプションで指定 `-sideChain`
### リガンド分割並列計算
mpiexec, あるいはmpirunコマンドをrbdockコマンドに付加して実行すると, 機能が追加される. 

分割数は, 並列数と島数から自動的に割り当てられる.
### 非極性水素付加
rbdockのコマンドオプションで指定 `-h`あるいは`-attachH`
### 島モデルGA
mpiexec, あるいはmpirunコマンドをrbdockコマンドに付加して実行すると, 機能が追加される. 

島数はrbdockのコマンドオプションで指定 `-l iIsland`

移住機能は `-m`
### MM-GB/PBSA
rbdockのコマンドオプションで指定 `-g mmgbpbsa_starting_point`

オプション`-g`の後に, 何回目のGAからMM-GB/PBSAを使用するかを指定する.
- `-g 3` : ３回目のGAからMM-GB/PBSAが実行され, 1,2回目のGAはrDockのスコアを使用する.
- `-g 0` : 全てのスコアがMM-GB/PBSAで計算される.
### タンパク質座標出力
rbdockのコマンドオプションで指定 `-outputRec`
### スコア計算分割並列計算
rbdockのコマンドオプションで指定 `-M`
