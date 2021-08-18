# NewBdG1
本プログラムはCuO2面のバルク層と表面層から成る2層モデルのみ計算可能。

### mwf1
カーパリネロ法によりエネルギー最小の波動関数(多価)を計算する。

### single_shot
LAPACK の対角化サブルーチンを一度だけ用いてxiの因子を取り出した波動関数を計算する。

### mwf2
chi の計算し, 一価波動関数を得る。また、電流分布を計算する。

## コンパイル
Intel® Fortran Compiler を使用。
```
make mwf1
make single_shot
make mwf2
```

OPT=0をつけると最適化なしでデバック用のコンパイル           \
OPT=1をつけると最適化(-O3)ありでコンパイル                       \
OPT=2をつけると最適化(-O3)ありで並列(OpenMP)もありのコンパイル   \
※指定のない場合はOPT=1                                         \
※オプション変更して実行する際は`make clean`を実行してから行う。
## 実行
./config/vconf を用意し,
mwf1, single_shot, mwf2の順番で実行する。
### 実行コマンド例
```
./mwf1.exe > 1.out &
```
```
./single_shot.exe > 2.out &
```
```
./mwf2.exe > 3.out &
```

## 入力ファイル
### mwf1
./config/vconf                                            

### single_shot
./config/vconf                                            \
fort.100                                                \
fort.101                                                

### mwf2
./config/vconf                                            \
fort.50 - fort.58                                         



## 出力ファイル
拡張子plt, dat はgnuplot等での表示用データ。gnuplotでの使用例については後述。 \
fort.*** は別で実行するプログラムで使用する変数のバイナリデータ。

### mwf1
| file名 | 説明 |
| ---- | ---- |
| xi.plt | スピン等をgnuplotで表示するためのフォーマット |
| init_xi.dat | 初期値に使用したスピン角$\xi$  |
| cp_xi.dat | car-parrinello計算後のスピン角$\xi$  |
| cp_spin.dat | car-parrinello計算後のスピン$S$  |
| cp_spin_n.dat | car-parrinello計算後のスピン$S$(規格化)  |
| cp_eg.dat | car-parrinello計算後のエネルギー固有値$E_n$  |
| fort.100 | エネルギー固有値$E_n$ real(8), dimension(4\*n_acc)  |
| fort.101 | 波動関数$u_{j\sigma},v_{j\sigma}$ complex(8), dimension(4\*n_acc, 4\*n_acc)  |

### single_shot
| file名 | 説明 |
| ---- | ---- |
| ss_spin.dat | single_shot計算後のスピン$S$  |
| ss_spin_n.dat | single_shot計算後のスピン$S$(規格化)  |
| ss_eg.dat | single_shot計算後のエネルギー固有値$E_n$  |
| fort.50 | (mwf2で使用)波動関数\tilde{u}_{j\sigma},\tilde{v}_{j\sigma}$ complex(8), dimension(4\*n_acc, 4\*n_acc)  |
| fort.51 | (mwf2で使用)エネルギー固有値$E_n$ real(8), dimension(4\*n_acc)  |
| fort.52 | (mwf2で使用)スピン角$\xi$ real(8), dimension(n)  |
| fort.54 | (mwf2で使用)$\uparrow$スピン電子数 real(8), dimension(n)  |
| fort.55 | (mwf2で使用)$\downarrow$スピン電子数 real(8), dimension(n)  |
| fort.56 | (mwf2で使用)$\Delta$ complex(8), dimension(n)  |
| fort.57 | (mwf2で使用)スピン Sx Sy Sz real(8), dimension(3\*n)  |
| fort.58 | (mwf2で使用)化学ポテンシャル$\mu$ real(8), dimension(2)  |
| fort.110 | エネルギー固有値$E_n$ real(8), dimension(4\*n_acc)  |
| fort.111 | 波動関数$u_{j\sigma},v_{j\sigma}$ complex(8), dimension(4\*n_acc, 4\*n_acc)  |

### mwf2
| file名 | 説明 | 
| ---- | ---- |
| chi.plt | chiや電流のgnuplotで表示するためのフォーマット  (巻き数はバルクのスモールポーラロンのだけ表示) |
| chi_sb_init.plt | chi.plt に表面の巻数を追加(初期値で使用した巻き数) |
| chi_sb_conv.plt | chi.plt に表面の巻数を追加(fchi-optimizerで収束したchiの巻き数) |
| eta.dat | $\eta$ ($\xi$から反強磁性成分を除いた角度) |
| wp_eta.dat | $\eta$の巻数(+1) |
| wm_eta.dat | $\eta$の巻数(-1) |
| wp_chi.dat | $\chi$の巻数(+1) |
| wm_chi.dat | $\chi$の巻数(-1) |
| init_chi.dat | 初期値に使用した$\chi$ |
| fchi.dat | fchi計算後の$\chi$ |
| current.dat | fchi計算後の電流分布 |
| fchi_eg.dat | fchi計算後のエネルギー固有値$E_n$ |
| fchi_spin.dat | fchi計算後のスピン$S$ |
| fchi_spin_n.dat | fchi計算後のスピン$S$(規格化) |
| fort.70 | (モンテカルロで使用)波動関数$\tilde{u}_{j\sigma},\tilde{v}_{j\sigma}$, complex(8), dimension(4\*n_acc, 4\*n_acc)  |
| fort.71 | (モンテカルロで使用)エネルギー固有値$E_n$ real(8), dimension(4\*n_acc)  |
| fort.72 | (モンテカルロで使用)スピン角$\xi$ real(8), dimension(n)  |
| fort.73 | (モンテカルロで使用)ベリー位相$\chi$ real(8), dimension(n)  |
| fort.74 | (モンテカルロで使用)$\uparrow$スピン電子数 real(8), dimension(n)  |
| fort.75 | (モンテカルロで使用)$\downarrow$スピン電子数 real(8), dimension(n)  |
| fort.76 | (モンテカルロで使用)表面層にある渦の位置 real(8), dimension(2, n_pole)  |
| fort.77 | (モンテカルロで使用)表面層にある渦の巻数 integer, dimension(n_pole)  |
| fort.78 | (モンテカルロで使用)表面層にある渦の数`n_pole`  integer |
| fort.79 | (モンテカルロで使用)化学ポテンシャル$\mu$ real(8), dimension(2)  |
| fort.120 | エネルギー固有値 real(8), dimension(4\*n_acc)  || fort.121 | 波動関数 complex(8), dimension(4\*n_acc, 4\*n_acc)  |
| fort.121 | 波動関数$u_{j\sigma},v_{j\sigma}$, complex(8), dimension(4\*n_acc, 4\*n_acc)  |
| fort.122 | $\chi$を含まない波動関数$e^{\frac{i}{2}\xi_j}\tilde{u}_{j\sigma},e^{\frac{i}{2}\xi_j}\tilde{v}_{j\sigma}$, complex(8), dimension(4\*n_acc, 4\*n_acc)  |


## gnuplot
### スピンの表示例
```
gnuplot> load "xi.plt"
gnuplot> splot "cp_spin.dat" with vector
```

### 電流の表示例
```
gnuplot> load "chi.plt"
gnuplot> splot "current.dat" with vector
```

### 巻数のあるループの表示例
```
gnuplot> load "chi_sb_conv.plt"
gnuplot> splot "wp_chi.dat" with vector, "wm_chi.dat" with vector
```

## その他
スタックメモリのサイズ制限によりエラーが出ることがある。その場合は
```
ulimit -s unlimited
```
を行う。

