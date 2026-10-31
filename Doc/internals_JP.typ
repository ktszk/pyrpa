//文章フォーマット設定
#set text(font: ("Times New Roman", "Liberation Serif"), size: 11pt)
#show regex("[　-龥＀-￯]"): set text(font: ("Yu Gothic", "Noto Sans CJK JP"), lang: "ja")
#set page(paper: "a4", margin: (x: 2cm, y: 2cm))
#set math.equation(numbering: numbering.with("(1)"), supplement: "式")
#set heading(numbering: "1.")
#import "callgraph.typ": module-graph-plibs, module-graph-kernels, option-graphs

#align(center)[
  #text(size: 18pt, weight: "bold")[pyrpa 開発者向け資料]
  #v(0.5em)
  #text(size: 11pt)[内部関数・実装・検証・変更履歴]
]

#v(1em)

#outline(title: "目次", depth: 2, indent: auto)

#pagebreak()

= この資料について

利用者向けの説明は `Doc/manual_JP.typ` にあります。この資料には, 次のものをまとめています。

- `main.py` からは設定せず, `libs.plibs` / `libs.flibs` の関数を直接呼ぶときに使う引数
- ビルドと実装の内部事情
- ファイル形式の互換性
- 精度の検証結果
- テストの内容と変更履歴

マニュアルの各節に対応する項目は, 同じ節の名前で参照しています。

= コードの構成

- `main.py`: 設定（上部の変数）と実行部。実行部は `setup()`（入力の検査, 模型の読み込み, 格子・相互作用・化学ポテンシャルの準備）と, `option` ごとの `run_*` 関数（対応表 `RUNNERS`）からなります。
- `libs/plibs`: 計算の Python 側（模型の読み込み, 対称性, ドライバ）。図は描きません。
- `libs/flibs`: Fortran カーネル `libfmod.so` の ctypes ラッパー。
- `libs/plots.py`: matplotlib による図。
- `libs/src/*.f90`: Fortran カーネル（$chi$, FLEX, Eliashberg, 伝導, Eilenberger など）。

#table(
  columns: (auto, 1fr),
  table.header([*モジュール*], [*内容*]),
  [`plibs._lattice`], [格子ベクトル, ホッピングの読み込み（全 `ftype`）, $bold(k)$・$bold(r)$ メッシュ, BZ],
  [`plibs._symmetry`], [点群操作, 対称性の検出, 既約 $bold(k)$ 点],
  [`plibs._bands`], [固有値, 化学ポテンシャル, フェルミ面, 断面積],
  [`plibs._response`], [スピン・ペアリング感受率のドライバ, ギャップの形状因子, 超伝導状態の $chi$],
  [`plibs._calc`], [FLEX・Eliashberg・スペクトル・CPA のドライバ, `sigma.bin` / ギャップの入出力],
  [`plibs._eliash_nosym`], [対称性を仮定しない線形 Eliashberg（Zeeman, FFLO）],
  [`plibs._freq_grid`], [IR 基底の周波数格子],
  [`plibs._conductivity`], [Boltzmann・久保公式の伝導, ホール係数, 緩和時間モデル],
  [`plibs._optical`, `_color`], [光学定数と物質の色],
  [`plibs._fermiology`], [サイクロトロン質量, dHvA, 経路上の有効質量],
  [`plibs._impurity`], [不純物スーパーセルと unfold したスペクトル],
  [`plibs._nmr`], [超伝導状態の NMR],
  [`plibs._phonon_pairing`], [EPA からのフォノン媒介の対相互作用],
  [`plibs._wannier_io`], [EPA の読み込み, U/J パラメータセット, Wannier 実空間形式の入出力],
  [`plibs._maxent`], [MaxEnt 解析接続（`ana_cont`）],
  [`plibs._eilenberger*`], [準古典 Eilenberger（一様系・モデル FS, 表面, 渦・渦格子, スピン行列版）],
  [`plibs._units`], [物理定数（`sw_unit`）],
  [`flibs._loader`], [`libfmod*.so` の選択と読み込み, OpenMP の設定],
  [`flibs._hamiltonian`, `_green`], [$H(bold(k))$, 固有値, グリーン関数],
  [`flibs._susceptibility`, `_interaction`], [$chi_0$, RPA, 相互作用頂点],
  [`flibs._selfenergy`, `_eliashberg`], [FLEX の自己エネルギー, 線形・非線形 Eliashberg],
  [`flibs._transport`], [伝導のカーネル],
  [`flibs._freq`, `_grid`, `_impurity`, `_eilenberger`], [周波数変換, 渦のグリッド, 不純物, Riccati のカーネル],
)

= `main.py` の `option` と呼び出す関数

どの `option` でも, 最初に `setup()` が模型の読み込み（`plibs.import_hoppings` / `import_MLO_hoppings` / `read_HamRsMLO`）, 格子（`plibs.get_ptv`, `get_bvec`）, 必要なら化学ポテンシャル（`plibs.get_mu`）, 反転操作（`plibs.inversion_op`）, $chi$ の軌道対（`plibs.get_chi_orb_list`）, 相互作用行列（`flibs.gen_SCmatrix`, `gen_SCmatrix_orb`）を用意します。その後に `option` ごとの `run_*` 関数が呼ぶ主な関数は次の通りです。

#table(
  columns: (auto, 1fr),
  table.header([*option*], [*主な関数*]),
  [0 `BAND`], [`plibs.mk_klist`, `get_eigs` → `plots.plot_band`],
  [1 `DOS`], [`flibs.gen_dos` → `plots.plot_dos`],
  [2 `FERMI_2D`], [`plibs.get_eigs_2d`, `get_kf_points`, `get_colors` → `plots.plot_FS`],
  [3 `FERMI_3D`], [`plibs.gen_3d_surf_points` → `plots.plot_3d_surf`],
  [4 `SPECTRUM`], [`plibs.spectrum_on_path` → `plots.plot_path_map`],
  [5 `CONDUCTIVITY_BT`], [`plibs.calc_conductivity_Boltzmann`, `get_hall_coe`],
  [6 `CONDUCTIVITY_PT`], [`plibs.calc_conductivity_lrt` → `plots.plot_optical`],
  [7, 10], [`plibs.mk_qlist`, `calc_path_spectrum`],
  [8], [`plibs.chis_q_point` → `plots.plot_chis_qpoint`],
  [9, 11], [`plibs.chis_qmap`, `phi_qmap` → `plots.plot_qmap`],
  [12, 13], [`plibs.sc_gap_on_mesh`, `chis_spectrum_sc`, `chis_q_point_sc`, `bdg_bands_path`],
  [14 `FLEX`], [`plibs.calc_flex` / `calc_flex_soc`],
  [15 `LIN_ELIASHBERG`], [`plibs.calc_lin_eliashberg_eq` / `calc_lin_eliash_soc` / `calc_lin_eliash_nosym`],
  [16 `GAP_FUNCTION`], [`plibs.output_Fk`],
  [17 `BAND_FILLING`], [`plibs.report_band_filling`],
  [18 `CYCLOTRON_MASS`], [`plibs.get_mass`],
  [19 `DHVA`], [`plibs.get_dhva_band` → `plots.plot_dhva`],
  [20 `ELECTRON_MASS`], [`plibs.band_mass_path`, `report_band_mass`],
  [21 `SPECTRUM_IMPURITY`], [`plibs.random_impurity_sites`, `calc_impurity_spectrum`],
  [22 `SIGMA_CPA`], [`plibs.calc_cpa`],
  [23 `NONLIN_ELIASHBERG`], [`plibs.calc_eliashberg_eq`],
  [24 `EILENBERGER`], [`plibs.calc_eilenberger`, `calc_penetration_depth`, `calc_fs_penetration`, `calc_pauli_limit`, `calc_spin_pauli`, `calc_free_energy`],
  [25 `EILENBERGER_SURFACE`], [`plibs.calc_surface`, `calc_surface_dvector`],
  [26 `EILENBERGER_VORTEX`], [`plibs.calc_vortex`, `calc_vortex_lattice_periodic`, `calc_vortex_lattice_sc`, `calc_vortex_current`, `calc_vortex_maxwell`, `calc_vortex_chiral`, `calc_vortex_dvector`],
  [27 `NMR_SC`], [`plibs.nmr_qmesh`, `nmr_sweep`, `write_nmr_dat` → `plots.plot_nmr`],
)

= 関数の相関図

`Doc/gen_internals.py` がソースを静的に解析して作る呼び出しグラフ（doxygen の call graph に相当）です。描画は Typst パッケージ diagraph（Graphviz を内蔵）で行うので, `dot` コマンドは要りません。

- 色: #box(fill: rgb("#e8e8e8"), inset: 2pt)[`main.py`], #box(fill: rgb("#dbe9f7"), inset: 2pt)[`plibs`], #box(fill: rgb("#dcf0dc"), inset: 2pt)[`flibs`], #box(fill: rgb("#fbe3c8"), inset: 2pt, radius: 3pt)[Fortran（角丸, 括弧内はソースファイル）]
- 矢印は「呼び出す側 → 呼び出される側」です。モジュールの図では, 線の太さが関数同士の呼び出しの数を表します。
- `option` ごとの図は `main.py` の `run_*` 関数から 5 段まで, 40 個までの関数を描きます。それより先がある関数は点線の「…」で示します。
- `flibs._loader` の引数変換（`i64`, `dbl`, `c_array` など）は全ての呼び出しに現れるので省いています。
- 静的解析なので, 関数を変数に入れて呼ぶ場合などは拾えません。

== `plibs` のモジュール間の依存

#module-graph-plibs

== `plibs` → `flibs` → Fortran

#module-graph-kernels

== `option` ごとの呼び出し

#option-graphs

= 関数リファレンス

`libs.plibs` と `libs.flibs` の公開関数（各サブモジュールの `__all__`）の一覧です。説明は各関数の docstring の最初の 1 文（`@brief` があればそれ）で, `Doc/gen_internals.py` が自動生成します（コードと一致させるため英語のままです）。引数の型と詳しい意味は docstring を参照してください。API を変えたら `python Doc/gen_internals.py` で作り直してください。

#include "api_reference.typ"

= ライブラリ関数の引数（`main.py` からは設定しないもの）

== 線形 Eliashberg の固有値ソルバー（`arnoldi_m`）

`flibs.linearized_eliashberg` と `flibs.linearized_eliashberg_soc` は `arnoldi_m`（既定 10）を取ります。`arnoldi_m > 0` で Arnoldi 法, `0` で shift + deflation 付きのべき乗法です。`main.py`（`plibs.calc_lin_eliashberg_eq`, `calc_lin_eliash_soc`）は既定値を使うので, 常に Arnoldi 法で解いています。

SOC がある場合は Arnoldi 法を使う必要があります。べき乗法にはブランチ追跡が無く, パリティセクターの支配的な状態（奇周波数の相手になり得る）に収束します。Arnoldi 法は `gap_sym` のスピン性格に一致するブランチを返します。

対称性を仮定しない版（`flibs.linearized_eliashberg_nosym`）の Krylov 次元は `main.py` の `nosym_arnoldi_m` です。

== MaxEnt による解析接続（`sw_maxent`）

`plibs.spectrum_on_path`（option=4）と `plibs.output_Fk`（option=16）は `sw_maxent=True`（既定 False）で Padé の代わりに MaxEnt を使います。`ana_cont` パッケージが必要です。`main.py` は Padé だけを使います。IR 基底の $Sigma$ から接続する場合, MaxEnt は `ir_eps` によらず安定です。

== FLEX の静的成分の減算（`sub_sigma`）

`flibs.mkself`, `mkself_soc` と非線形 Eliashberg のラッパーは `sub_sigma` を取ります。各反復で $Sigma$ から静的成分を差し引き, バンドの静的なシフトを化学ポテンシャル側に吸収させてフェルミ面を安定させます。

- `1`: $Sigma(bold(k), i omega_0)$ の Hermite 部分（既定。`main.py` はこれを使います）
- `2`: 周波数平均（一様格子専用。IR 基底では使えません）
- `0`: 減算なし

厳密な Hartree--Fock 項だけの減算ではないので, 絶対的なバンドシフトを議論するなら `sub_sigma=0` を使ってください。$chi$ の軌道対を相関軌道に限った模型（`chiolist`）では $Sigma$ の静的部分が相関軌道にだけ入り, $epsilon_d - epsilon_p$ を実質的に動かすので, 二重計上の扱いに注意してください。AiiDA インターフェースでは `options` の `sub_sigma` で指定できます。

== `sigma.bin` の再グリッド（`plibs.regrid_sigma_bin`）

`sw_in_self=True` で読み込む `sigma.bin` の `Nw` が今回の計算と違う場合, `plibs._auto_regrid_sigma_seed` が `regrid_sigma_bin(temp_new=temp, Nw_new=Nw, w_scale=temp_old/temp)` を呼び, ファイルをその場で書き換えます。手で呼ぶ場合は次の通りです。

```python
plibs.regrid_sigma_bin(temp_new=..., Nw_new=..., w_scale=...)
```

- `w_scale=temp_old/temp_new` は, 同じ `Nw` の `sigma.bin` をそのまま読んだときの周波数軸の圧縮を再現します（臨界点の近くで冷やす場合）。`w_scale=1`（関数の既定）は忠実な補間です。
- 補間は共役対称性 $Sigma_(l m)(-i omega) = Sigma_(m l)(i omega)^*$ で負の周波数に拡張した 3 次スプラインです。循環畳み込みの折返しを含む最終点は除き, 窓の外は内側でフィットした $c_0 + c_1\/(i omega)$ で補います。
- SOC 版（`mkself_soc`）の `sigma.bin` も同じ形式で, 同じ関数が使えます。
- 元の格子は `sigma.bin` 末尾のメタデータ（温度・`Nw`・軌道数）から取得します。明示した値がこれと矛盾するとエラーになります。メタデータの無いファイルについては「ファイル形式と互換性」の節を参照してください。

== CPA の濃度とポテンシャル

`plibs.calc_cpa(..., x_cpa=0.5, VB_shift=1.0)` は 2 元合金の CPA です。`x_cpa` は B 種の濃度, `VB_shift` は B 種のオンサイト準位のずれ [eV] です。`main.py`（option=22）は上の固定値を渡しています。

== ギャップの形状因子（`set_harmonic_lattice` ほか）

`main.py` は形状因子のための格子を自動で登録します（`plibs.set_harmonic_lattice(avec, rvec)`; $s^plus.minus$ のシェルは `s_rcut=` で, `main.py` の `gap_rcut` がこれに渡ります）。`plibs.lattice_harmonic()` で単独の調和関数を直接構成でき, `gap_symms(..., avec=, rvec=)` で格子を明示的に渡すこともできます。実際に使われたシェルの $|bold(R)|$ は `get_initial_gap` がログに出力します。

直交系の 1 サイト格子では, マニュアルの表の全ての指定・任意の軸比について, 以前の簡約座標の式（$cos(2 pi k_x) - cos(2 pi k_y)$ など）をビット単位で再現します。

== 渦格子の対称性（`calc_vortex_lattice_symmetry`）

`plibs.calc_vortex_lattice_symmetry` は `main.py` からは呼ばれないドライバです。Ichioka–Machida の格子自由エネルギーをセルの頂角とギャップ対格子方位 $theta_0$ について最小化し, *安定な渦格子対称性* とその磁場発展（例: $H_(c 2)$ 近傍での $d$ 波の三角→正方転移）を決めます。Wannier FS を与えると $theta_0$ は結晶全体（FS ＋ギャップ）を剛体回転させるため, フェルミ速度の異方性も対称性の選択に寄与します。

= ビルドの内部

- *`SL` の自動判定.* 読み込まれているモジュール（`MKLROOT`, `AOCL_ROOT`）に合わせて判定しています。root が設定されていないライブラリを選ぶと, BLAS のシンボルが未解決のまま import 時に初めて落ちる `libfmod.so` ができてしまうためです。
- *gfortran と BLAS の整数型.* gfortran は同じファイル内で同じ外部手続きを異なる整数の型で呼ぶとエラーにするため, BLAS/LAPACK の整数引数はファイルごとに型を揃えています（`fmod.f90` では `c_int64_t`）。
- *CMake.* 参照実装は makefile 側です。C++ コンパイラが必要なのは pocketfft 経路だけで, MKL 経路では CXX 言語自体を有効にしないため, `ifx` はあるが `icpx` が無い oneAPI 環境でも configure できます。`flang` と AOCC のモジュールを併用する場合, makefile は子プロセスの `CPATH` と `LIBRARY_PATH` の掃除も行いますが, CMake 側は行いません。

= FLEX と Eliashberg の実装

- *`sw_rescale_flex` のクランプ先.* $chi_0$ を $(1 - 10^(-4))\/S$ 倍と 1 の直下にしているのは意図的で, クランプ直後の $Sigma$ の過大応答が偽の固定点から抜け出す機構になります。0.98 程度に緩めると, $S < 1$ の正しい固定点がある近臨界の系でも偽の固定点に速く「収束」してしまうことを確認しています。
- *温度アニーリングの初期値.* 同じ `Nw` の `sigma.bin` をそのまま読むと周波数軸を $T_"old"\/T_"new"$ 倍に圧縮した初期値になります。臨界点の近くでは, これが忠実な補間よりむしろ速く収束することを確認しています。
- *`sigma.bin` の照合.* 読み込み側で `Nw`・軌道数・$bold(k)$ 点数を照合するのは, Fortran の read がレコードの大きいファイルを黙って誤った並びで読む事故を防ぐためです。
- *非線形 Eliashberg の振幅方向 Newton 加速.* 詳細は `libs/src/ffeliash.f90` の `sw_amp_newton` 関連のコメントを参照してください。
- *`sw_chi0_tail` の検証.* 厳密な Lindhard 関数に対し, $chi_0$ の収束次数が 1.0 から 2.0 に上がることを確認しています。
- *IR 基底での $Sigma$ の書き出し.* `self_en_wannier.npz` を一様な $omega_n$ の形式で書き出すとき, 基底で表せない定数 $Sigma(i omega -> infinity)$ は最も高い 3 つのサンプリング点の Hermite 部分から $1\/omega^2$, $1\/omega^4$ 項とともに推定して除いてから展開します。
- *`eil_gap_file` の Fourier 変換の符号.* $Delta_"orb"(bold(k)) = sum_bold(R) e^(-i 2 pi bold(k) dot bold(R)) Delta(bold(R))$ の負符号は, 書き出し側の `np.fft.ifftn` の厳密な逆変換です。書き出し→読み戻しの往復で検証しています。

= ファイル形式と互換性

- *`sigma.bin`.* 末尾に温度・`Nw`・軌道数のメタデータを追記しています（追記前のコードの読み込みには影響しません）。メタデータの無いファイルは自動の再グリッドの対象にならないので, `regrid_sigma_bin` に `temp_old`/`Norb`/`Nw_old` を指定してください。SOC 版の旧・要素別レコード形式は読めません。
- *`self_en_wannier.npz` の $bold(k)$ メッシュ.* `tau_mode='flex'` は FLEX の $bold(k)$ メッシュをファイルから読みます。記録のない古いファイルでは $bold(R)$ の範囲から推定するので, FLEX を `sw_out_self=True` で解き直してください。
- *ギャップファイルの軌道順序（`orb_conv`）.* `gap_wannier.npz` の軌道ペアは物理的な順序 $Delta_(a b)$ で格納します（書き出し時に Fortran 配列の反転ビューを一度だけ転置し, npz に `orb_conv` を刻印します）。刻印が無い旧形式も正しく読めますが, その旨が 1 行表示されるので, 気になる場合は書き出し直してください。この順序は無害な付け替えではありません。$Delta$ はエルミートですが対称ではなく, 転置で反対称成分の符号が反転します（000AsP の RPA 解では $max |Delta|$ の 63%。完全ギャップの $s_(plus.minus)$ ポケットがノードを持つように見えてしまいます）。`eil_gap_file` と `nmr_gap_file` で同じ扱いです。
- *U, J のパラメータセット.* TOML（`.toml`）のほか JSON（`.json`）も読めます。同名の `.toml` があればそちらを優先します。

= 精度の検証

== IR 基底

- *`ir_wmax` の既定値.* 000AsP, $T = 0.02$ eV で, $omega_"max" = 2.88$ eV では $mu$ が 0.27 meV ずれ, 6.0 eV で一様格子の外挿と $6 times 10^(-6)$ eV で一致しました。既定値を $2 times max|epsilon_bold(k) - mu|$ にしているのはこのためです。
- *一様格子の外挿との一致.* $chi_0$ は厳密な Lindhard 和と $10^(-8)$, FLEX の $mu$ と Stoner 因子は $10^(-5)$ 以下, 線形 Eliashberg の $lambda$ は $10^(-7)$〜$10^(-8)$, 非線形 Eliashberg の収束したギャップは $10^(-4)$ 以下で一致します。
- *一様格子の誤差.* 1 軌道 $16 times 16$, $T = 0.02$ eV, $lambda = 1.34$ の例では, $N_w = 4096$ でもギャップが 3% 小さく出ます。FLEX では $chi_0$ の tail 補正を入れても $Sigma$ の打ち切りで $O(1\/N_w)$ の誤差が残り, $N_w = 1024$ でも $mu$ が 1 meV 程度ずれることがあります。
- *Padé をサンプリング点から取る理由.* Padé は悪条件で, 入力に相対 $10^(-10)$ の乱雑な誤差を加えただけでスペクトルが数十 % 変わります。IR の値を一様な $omega_n$ で評価し直すと, 松原軸の上では一様格子より正確です（`square_2orb.hop`, $16 times 16$, $T = 0.05$ eV: 外挿との差は $G$ で $9 times 10^(-6)$, 一様格子 $N_w = 4096$ では $2 times 10^(-4)$）。しかし残る基底の打ち切り誤差は滑らかで系統的なため, Padé はそれを実軸の偽の構造として接続してしまいます。精度は `ir_eps` で決まり, FLEX の収束判定にはほぼよりません。
- *Padé の精度*（Richardson 外挿した一様格子のスペクトルとの相対 L2 差, $bold(k)$ 点の中央値 / 最大）:
  - `square_2orb.hop`: `ir_eps=1e-8` で 12% / 45%, `1e-12` で 4% / 34%, `1e-14` で 3% / 6%（一様格子 $N_w = 1024$: 1.5% / 3%）
  - 000AsP（$16 times 16 times 4$, $T = 0.05$ eV）: `1e-12`, `1e-14` で 4% / 13%（$N_w = 256$: 3% / 14%）

== MLO 基底

- テイル補正なしの $N_w = 128$（窓 40 eV）では, Cu の MLO 模型で $lambda$ が数倍ずれた例があります（バンドが約 $+30$ eV まで延びるため）。
- 規格化しない場合, ecalj の MLO はノルムが 1 ではなく, Cu の s 軌道では $1\/S_(a a)(0)^2$ が約 4000 になります。

== その他

- *`sw_dec_axis`.* 物理的な $bold(k)$ における $E(bold(k))$ が全ての `brav` について $10^(-14)$ で保存されることを確認しています。
- *dHvA（option=19）.* 幾何と単位変換は, 正方晶格子での $F(theta)=F(0)/cos(theta)$ を含む解析的な円筒フェルミ面テストで検証しています。物質模型では, 実験格子定数を用いた 10 軌道 FeS Wannier 模型で $B parallel c$ のとき $0.489$–$3.123$ kT となり, T. Terashima _et al._, Phys. Rev. B *99*, 134501 (2019) の Fig. 3(b) の GGA 振動数（およそ $0.5$–$3$ kT）と整合します。

== 参考: DFT のバンドと実験の比較

上の FeS の論文の実験値は $0.15$–$1.87$ kT で, 計算値より小さくなっています。この差は計算の不備ではなく想定されるもので, 鉄系に限った話でもありません。強相関金属では DFT（LDA/GGA）のバンド構造が電子相関によるバンドの narrowing と軌道依存のバンドシフトを含まないため, 一般にフェルミ面のポケットを系統的に大きく, 質量繰り込み（したがって $E_F$ での状態密度）を小さく見積もります。この偏りは本コードの対相関の計算にもそのまま伝播します。繰り込みを受けていない DFT バンドの上で評価したクーパー不安定性は弱く出るため, RPA/Eliashberg の固有値や $T_c$ は過小評価されます。

鉄系以外ではずれがさらに大きくなる例もあります。

- 重い電子系反強磁性体 $"CeRhIn"_5$ では, $4f$ 電子を遍歴として扱うバンド計算は large Fermi surface を与えますが, dHvA で観測されるフェルミ面は $4f$ を持たない参照物質 $"LaRhIn"_5$ とほぼ同じで, サイクロトロン質量だけが約 1 桁重くなります（H. Shishido _et al._, J. Phys. Soc. Jpn. *71*, 162 (2002)）。
- $"Na"_x "CoO"_2$ では, LDA が予言する 6 個の小さな $e'_g$ ホールポケット（D. J. Singh, Phys. Rev. B *61*, 13397 (2000)）が光電子分光では観測されず, $a_(1g)$ 由来のシートのみが見えます（H.-B. Yang _et al._, Phys. Rev. Lett. *92*, 246403 (2004)）。これを局所相関だけで説明できるかについては, 現在も議論が続いています。

= テスト

実行方法は `tests/README_JP.md` とマニュアルの付録にあります。全体で約 7 分（ifx 版, 24 スレッド）かかり, その大半は `test_eliash_nosym.py`（全ブリルアンゾーンで解く）と `test_eilenberger.py`（2 次元グリッド上の自己無撞着場の問題）です。`test_aiida.py` は `aiida-core` が, `test_freq_grid.py` は `sparse-ir` が無い環境ではスキップされます。

#table(
  columns: (auto, 1fr),
  align: (left, left),
  table.header([*ファイル*], [*確認している内容*]),
  [`test_bands.py`], [バンドとフェルミ面の幾何: Wannier90 `_hr.dat` の読み込みと `gen_ham` の位相, 3 次元フェルミ面のための第一ブリルアンゾーン（面, 折り畳み）, 非直交（MLO）基底の群速度（overlap 項, エルミート性）, 直交・MLO 基底での逆有効質量（バンド間項, overlap 項, 縮退）, サイクロトロン軌道と dHvA 極値軌道（断面の幾何, 開軌道の棄却, ゾーン還元）],
  [`test_transport.py`], [直流極限での Kubo$arrow.l.r$ボルツマン一致, Wiedemann--Franz 則, 対称化されたバンド間熱流頂点, バンド縮退, $bold(k)$ 点の重み, 弱磁場ホール応答, 緩和時間モデル（EPA, FLEX）],
  [`test_optical.py`], [光学定数と色: 反射率から色への変換（CIE XYZ, D65, sRGB）, 屈折率, 反射率, 吸収係数, 透過, Kubelka--Munk],
  [`test_impurity.py`], [不純物スーパーセルと unfold したスペクトル（位相の符号, 総和則）, 非直交基底での unfold と CPA],
  [`test_rpa_flex.py`], [RPA / FLEX / Eliashberg の構成要素: 相互作用頂点, RPA 行列代数, 厳密 Lindhard 和・対バブルとの照合による $chi_0$ と $phi_0$, tail 補正 $chi_0$, `mkself` / `mkself_soc` と二重計上の相殺, 非線形 Eliashberg の Nambu 解, 自己エネルギーの再グリッド, 時間反転対称性下の既約 $bold(k)$ メッシュ; 非直交（MLO）基底の Green 関数・FLEX・Eliashberg（行列の Dyson 方程式, スケールに対する共変性, 軌道の規格化）; EPA からのフォノン媒介の対相互作用],
  [`test_eliash_nosym.py`], [対称性を仮定しない線形 Eliashberg（Zeeman 場, FFLO）, SOC 版との整合（裸の項 $-V\/2$, triplet の初期値とパリティ射影, バンド基底の出力）, 多軌道 triplet の厳密カーネルとの照合, SOC 下のブランチ追跡],
  [`test_freq_grid.py`], [IR 基底: 基底の再現性, 厳密な $chi_0$, FLEX・線形 / 非線形 Eliashberg・SOC・`sw_nosym`・フォノン頂点が一様格子の $N_w -> infinity$ 極限と一致すること, IR での tail 補正の無視, サンプリング点からの Padé 接続],
  [`test_nmr.py`], [超伝導状態の動的感受率（非直交基底の BdG を含む）と NMR: BCS ギャップ内挿, $bold(q)$ メッシュの折り畳み, 自己無撞着 Eliashberg ギャップの読み込み, Knight shift（Yosida と triplet）, Hebel--Slichter ピーク, ライン節の $T^3$ 則],
  [`test_eilenberger.py`], [準古典 Eilenberger / Riccati ソルバー: 均一系の極限, Fortran カーネルと Python 参照実装の一致, 表面 Andreev 束縛状態, 渦糸および自己無撞着渦糸格子, 凝縮自由エネルギー, model / Wannier フェルミ面, Pauli 制限, triplet $d$-vector texture; 渦ソルバーの Fortran グリッドカーネルと Riccati の chord カーネル],
  [`test_main.py`], [`main.py` の Eilenberger モードでの化学ポテンシャルの扱い（Wannier フェルミ面と固定 $mu$）],
  [`test_aiida.py`], [AiiDA インターフェース（`examples/aiida`）: 各ステップと workflow が直接の `plibs` / `flibs` 呼び出しと一致すること（RPA / FLEX / Eliashberg, MLO, SOC, Eilenberger）, `options`, 格子の記録, 来歴の記録],
  [`test_tools.py`], [他ファイルが照合先とする参照実装（Fermi / Bose 恒等式, 閉形式のモデル, 厳密 $chi_0$）, CPU に応じたライブラリの選択（AMD, Intel）, gfortran 版の実行環境（`OMP_WAIT_POLICY`, 標準出力のバッファリング）の設定],
)

= 変更履歴

以下は過去のバージョンで得た結果を変えた修正です。該当する計算は再計算してください（それ以前の変更は git の履歴を参照してください）。

- *非直交 MLO 基底（ftype=3,4）の FLEX・Eliashberg.* 自己エネルギーを含む計算で $G^(-1)$ に重なり行列 $S$ が入っておらず（$G^(-1) = i omega + mu - H - Sigma$）, 粒子数の補正項とも整合していませんでした。Cu の MLO 模型では FLEX の $mu$ が裸の値から約 $+2.9$ eV ずれ, Stoner 因子が負になります。現在は $G^(-1) = (i omega + mu) S - H - Sigma$, 粒子数は $"tr"(S G)$ です。*以前の MLO の FLEX, FLEX を使った Eliashberg, 自己エネルギー付きスペクトルは全て誤りです。* RPA（$Sigma = 0$）は固有ベクトルから $G_0$ を作っていたので正しく, 変わりません。
- *MLO 軌道の規格化.* ecalj の MLO はノルムが 1 ではなく（Cu: $S_(a a)(0)$ が s で 0.016, d で 0.565）, 生の軌道に置いた U は実効的に $U\/S_(a a)(0)^2$ 程度で効いていました（Cu の RPA, U = 1 で Stoner 因子 85）。現在は読み込み時に各軌道をノルム 1 に規格化します（バンドは不変）。U, J, 軌道の重み, 軌道表示の $Sigma$ や $Delta$ の値は変わり, 規格化前に保存したファイルとは規約が合いません。バンド間の速度行列要素（Löwdin 直交化を通るもの）は数 % 変わり, 光学伝導度・反射率もそれに応じて変わります（Cu で $Delta R lt.eq 0.02$）。群速度とそれだけを使う量（Boltzmann 輸送, 有効質量, dHvA）は変わりません。
- *MLO での `sw_chi0_tail`.* $tau = 0$ での $G_0$ の跳びを単位行列と仮定していたため, 非直交基底では打ち切り誤差が $O(1\/N_w^2)$ にならず $O(1\/N_w)$ のままでした。現在は跳び $S^(-1)(bold(k))$ と密度行列の $bold(k)$ 畳み込みで扱います。
- *非線形 Eliashberg 方程式（option=23）.* 超伝導状態の正常グリーン関数を, 前の反復の $F$ を使った $G = G_0 (1 + Delta F^dagger)$ の部分反復ではなく, Nambu の Dyson 方程式 $G = [G_0^(-1) + Delta G_0^dagger Delta^dagger]^(-1)$ から直接求めるようにしました（旧方式は $Delta$ がおよそ $pi T$ を超えると発散）。`sw_self=True` では各反復で $Sigma$ の静的部分を差し引き（`sub_sigma`, $mu$ を固定しているため）, $Sigma$ も混合率 `pp` で混合します。$Delta$ が NaN/Inf になった場合は収束扱いにせずエラーで止まります。収束したギャップと $Sigma$ は以前と変わります。
- *Fortran 側の論理型.* Python から 1 バイトの `c_bool` で渡しているのに Fortran 側が既定の 4 バイトの `logical` で受けていた引数がありました（超伝導状態の $chi$ のスピンチャンネル `sw_spsym`, $chi$ の先頭固有値とトレースの切り替え `sw_eig` など）。残りの 3 バイトは不定なので, 実行環境によっては誤った値として読まれていた可能性があります。option=12, 13, 27 の結果に疑問がある場合は再計算してください。

同じ時期に追加した機能（以前の結果は変えません）: IR 基底の周波数格子（`freq_grid='ir'`）, 対称性を仮定しない線形 Eliashberg（`sw_nosym`）, AiiDA インターフェース, gfortran でビルドしたライブラリの自動設定。
