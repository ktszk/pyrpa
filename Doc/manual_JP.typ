//文章フォーマット設定
#set text(font: ("Times New Roman", "Liberation Serif"), size: 11pt)
#show regex("[　-龥＀-￯]"): set text(font: ("Yu Gothic", "Noto Sans CJK JP"), lang: "ja")
#set page(paper: "a4", margin: (x: 2cm, y: 2cm))
#set math.equation(numbering: numbering.with("(1)"), supplement: "式")
#set heading(numbering: "1.")
//外部パッケージ
#import "@preview/physica:0.9.3": *

// 読み飛ばしてよい補足（実装・精度検証・細かな注意）
#let detail(body) = block(fill: luma(246), stroke: (left: 2pt + luma(170)), inset: 8pt, width: 100%)[
  #text(size: 9.5pt)[#text(weight: "bold")[詳細] #h(0.4em) #body]
]
// 重要な注意
#let caution(body) = block(fill: luma(235), inset: 8pt, radius: 4pt, width: 100%)[#body]

#align(center)[
  #text(size: 18pt, weight: "bold")[pyrpa マニュアル]
  #v(0.5em)
  #text(size: 11pt)[2026年]
]

#v(1em)

#outline(title: "目次", depth: 2, indent: auto)

#pagebreak()

= pyrpaについて

pyrpa は, タイトバインディング模型（Wannier90 や ecalj の出力など）を入力として, バンド構造から輸送係数, 磁気揺らぎ, 超伝導までを計算・可視化する Python + Fortran のコードです。主に第一原理計算の結果のポスト処理を想定しています。

主な機能は次の通りです。

- *電子構造*: バンド, 状態密度（DOS）, フェルミ面（2D / 3D）, スペクトル関数 $A(bold(k), omega)$, 有効質量・サイクロトロン質量・dHvA 振動数, バンドフィリング
- *輸送・光学*: ボルツマン理論・久保公式による伝導係数（電気伝導率, 熱伝導率, ゼーベック係数など）, 弱磁場ホール係数（バンド分解つき）, 反射率と物質の色
- *揺らぎと超伝導*: RPA のスピン・電荷・ペアリング感受率, FLEX 自己エネルギー, 線形・非線形 Eliashberg 方程式, 対称性を仮定しない Eliashberg（Zeeman 場, FFLO）, 超伝導状態の動的スピン感受率と NMR
- *準古典理論*: Eilenberger/Riccati ソルバによる $T_c$, 磁場侵入長, 表面 Andreev 束縛状態, 渦・渦格子
- *その他*: 不純物・CPA, 非直交 MLO 基底（ecalj）, IR 基底の周波数格子, AiiDA による来歴付き実行

設定はすべて `main.py` の上部（ライブラリの import 行より上）にある変数を書き換えて行います。

== このマニュアルの読み方

初めて使う場合は, 次の順に読んでください。

+ 「クイックスタート」で, ビルドから最初の計算までを試す
+ 「計算モード」の冒頭の一覧表で, 目的の計算に対応する `option` を探す
+ その `option` の節と, 「典型的な計算フローの例」を読む

灰色の帯がついた *詳細* の枠は, 計算の仕組みや特殊な場合の注意です。初めは読み飛ばして構いません。パラメータの完全な一覧は「各種パラメータの詳細説明」にあります。

内部関数・`main.py` から設定しないライブラリ引数・実装の詳細・精度の検証結果・変更履歴は, 開発者向けの資料 `Doc/internals_JP.typ` にまとめてあります。

== クイックスタート

+ Fortran カーネルをビルドします（コンパイラ・数値ライブラリ・CPU は自動判定されます。詳しくは「Fortran カーネルのビルド」）。
  ```bash
  cd libs && make && cd ..
  ```
+ `main.py` の上部で, 模型・計算モード・条件を設定します。例えば同梱の 1 軌道正方格子模型のバンドを描く場合:
  ```python
  fname, ftype, brav, sw_soc = 'inputs/square.hop', 1, 0, False
                      # 模型ファイル, 形式, 格子, SOC
  option = 0            # 計算モード: 0 = バンド図
  abc = [1., 1., 1.]    # 格子定数 [Å]
  fill = 0.45           # フィリング（SOC なしではスピンあたりの電子数）
  tempK = 200           # 温度 [K]
  ```
+ 実行します。結果は図として表示され, モードによっては `*.png`, `*.dat` などのファイルにも出力されます。
  ```bash
  python main.py
  ```

続けて `option = 2`（フェルミ面）, `option = 7`（スピン感受率）, `option = 15`（線形 Eliashberg）などに変えてみてください。必要なパラメータは各 `option` の節に書いてあります。

== 用語

#table(
  columns: (auto, 1fr),
  [*用語*], [*意味*],
  [タイトバインディング模型], [ホッピング $H(bold(R))$ で表した一電子ハミルトニアン。Wannier90 などで第一原理計算から作ります],
  [RPA], [乱雑位相近似。裸の感受率 $chi_0$ に相互作用 $U$, $J$ を足し上げてスピン・電荷感受率を求めます],
  [FLEX], [揺らぎ交換近似。RPA の揺らぎによる自己エネルギー $Sigma$ を自己無撞着に求めます],
  [Eliashberg 方程式], [揺らぎを媒介とする超伝導のギャップ方程式。線形化した方程式の最大固有値 $lambda$ が 1 になる温度が $T_c$ です],
  [Stoner 因子], [$chi_0 U$ の最大固有値。1 に達すると磁気秩序の不安定性（RPA の発散）を意味します],
  [松原周波数, `Nw`], [温度グリーン関数の虚周波数 $omega_n = (2n+1) pi T$ と, その点数],
  [IR 基底], [松原周波数を少数（20〜30 点程度）のサンプリング点で表す方法。`freq_grid='ir'` で使います],
  [MLO], [ecalj の muffin-tin 軌道に基づく非直交な局在基底。重なり行列 $S$ を伴います],
)


= 動作環境

== 必要なもの

- Python 3.11 以降（`U`, `J` のパラメータセットの読み込みに標準の `tomllib` を使います）
- Fortran コンパイラ（ifx / ifort, flang / amdflang, gfortran のいずれか）と BLAS/LAPACK（MKL, AOCL, OpenBLAS のいずれか）
- Python パッケージ:
  - `numpy`, `scipy`, `matplotlib`
  - `skimage`（フェルミ面 option=2,3, サイクロトロン質量 option=18, dHvA option=19 で使用）
  - `sparse-ir`（任意。IR 基底 `freq_grid='ir'` で使用。1.x, 2.x; 2.1.6 で確認）
  - `aiida-core`（任意。AiiDA インターフェース `examples/aiida/` でのみ使用）

計算の重い部分（$chi$, FLEX, Eliashberg など）は OpenMP 並列の Fortran カーネル `libs/libfmod.so` で行います。`libs/flibs`, `libs/plibs` はこれを呼び出す Python パッケージです。

== Fortran カーネルのビルド

=== 手元のマシンで使う場合

`libs` ディレクトリで `make` を実行するだけです。

```bash
cd libs
make                          # すべて自動判定
make FC=gfortran SL=OpenBLAS  # コンパイラと数値ライブラリを指定する場合
make ARCH=znver4              # CPU ターゲットを指定する場合
```

自動判定の内容は次の通りで, コマンドラインで指定した値が常に優先されます。

#table(
  columns: (auto, auto, 1fr),
  [*変数*], [*意味*], [*自動判定*],
  [`FC`], [コンパイラ], [`PATH` 上の `ifx`, `ifort`, `flang`, `amdflang`, `gfortran` の最初に見つかったもの],
  [`SL`], [BLAS/LAPACK], [`MKLROOT` が設定されていれば MKL, `AOCL_ROOT` なら AOCL, どちらも無ければ OpenBLAS],
  [`ARCH`], [CPU ターゲット], [AMD では `znver3`, Intel では `x86-64-v3`（どちらも AVX2 を使う可搬な設定）。AVX-512 を持つ CPU 専用にするなら `znver4` / `x86-64-v4`],
)

=== クラスタで使う場合

ログインノードと計算ノードで CPU が違うことがあるため, ビルドしたマシンの CPU に合わせると計算ノードで動かない（SIGILL で落ちる）ことがあります。その場合は `make auto` を使ってください。

#table(
  columns: (auto, 1fr),
  [*コマンド*], [*内容*],
  [`make auto`], [*クラスタでの推奨.* 計算ノードの CPU を調べ, CPU の種類ごとに `libfmod_<target>.so` をビルドします。実行時には `libs/flibs/_loader.py` が, そのノードで動く最良のライブラリを自動で選びます],
  [`make dispatch`], [CPU を調べずに, AVX2 用と AVX-512 用の 2 つを作ります（AMD: `znver3` + `znver4`, Intel: `x86-64-v3` + `x86-64-v4`; `DISPATCH_ARCHS` で変更可）],
  [`make`], [ライブラリを 1 つだけ作ります（上の「手元のマシンで使う場合」）],
)

CPU の一覧は `libs/detect_arch.sh` で確認できます。

```
$ sh libs/detect_arch.sh --table
NODE         CPU                                      -march
salmon       AMD Ryzen 7 5700G with Radeon Graphics   znver3
unagi        AMD Ryzen 9 7900 12-Core Processor       znver4
yakisoba     AMD Ryzen 9 7900X 12-Core Processor      znver4
```

結果は `libs/.arch_cache` にキャッシュされます（`--refresh` で再取得, `--local` で自機のみ）。`make auto ARCH_LIST="znver3 znver4"` のように一覧を直接与えることもできます。

AMD 上で MKL を使う場合は `MKL_ENABLE_INSTRUCTIONS=AVX512`（Zen 4）または `AVX2`（Zen 3）を export してください（Intel CPU では MKL が自動で選びます）。`OMP_NUM_THREADS` は決め打ちせず, スケジューラの値（`$SLURM_CPUS_PER_TASK` など）から設定してください。

#detail[
  *CPU の調べ方.* `detect_arch.sh` はスケジューラの台帳（Slurm の `sinfo`, PBS の `pbsnodes`, どちらも無ければ自機のみ）にある各ノードで `gcc -march=native` がその CPU を何と呼ぶかを問い合わせ, 重複を除いたターゲットを出力します。応答しないノードは無視せず, 可搬な既定値（AMD では `znver3`, Intel では `x86-64-v3`）を加えて警告を出します。`slurm.conf` に `Features=znver3` / `znver4` を設定しておけば, ジョブを投げずに `sinfo` の参照だけで済みます。

  *Intel のノード* は世代名（`sapphirerapids` など）ではなく, CPU フラグから決まる x86-64 のレベル（AVX-512 F/BW/CD/DQ/VL があれば `x86-64-v4`, AVX2/FMA/BMI2 なら `x86-64-v3`）で報告します。世代ごとの細かな拡張命令（AVX-VNNI, AMX など）の違いで, 古いノードに新しい世代向けのライブラリが渡るのを避けるためです。`x86-64-v3` は Zen でも動くので, AMD と Intel が混ざったクラスタはこれ 1 つで済みます。

  *実行時の選択.* `_loader.py` は CPU フラグを読み, 各ターゲットに必要なフラグが全て揃っている中で最上位の `libfmod_*.so` を読み込みます（無ければ `libfmod.so`）。`make ARCH=sapphirerapids` のように世代名でビルドしたライブラリも, その世代の拡張命令（AMX, AVX512-BF16 など）があるノードでだけ選ばれます。]

=== gfortran でビルドする場合

特別な設定は要りません。MKL と組み合わせる場合（`SL=MKL`）は GNU のスレッド層（`-lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core`）にリンクします。`_loader.py` は gfortran 版のライブラリを読み込む前に, 次の 2 つの環境変数を自動で設定します（ユーザーが値を設定していればそちらを優先）。

- `OMP_WAIT_POLICY=PASSIVE`: gfortran 版は OpenMP に libgomp を, MKL 版の numpy は libiomp5 を使います。同じプロセスで 2 つの OpenMP ランタイムが待機中のスレッドを回し続けてコアを奪い合い, 何十倍も遅くなるのを防ぎます。ifx/ifort 版は numpy と同じ libiomp5 なので影響を受けません。
- `GFORTRAN_UNBUFFERED_PRECONNECTED=y`: libgfortran は標準出力が端末でないと出力を終了までため込むため, FLEX の反復のログが途中で見られるようにします（AiiDA インターフェースはこのログで収束を判定します）。

=== CMake を使う場合

`libs/CMakeLists.txt` でも同じ判定・同じ成果物が得られます（IDE 連携などで CMake を使いたい場合）。

```bash
cmake -S libs -B build && cmake --build build -j     # = make auto
cmake -S libs -B build -DARCH=znver4                 # = make ARCH=znver4（libfmod.so 単体）
cmake -S libs -B build -DFC_COMPILER=gfortran -DSL=MKL -DARCH_LIST="znver3;znver4"
```

#detail[
  各変数は空なら自動判定です（`FC_COMPILER` は `PATH` から, `SL` は `MKLROOT`/`AOCL_ROOT` から, ターゲット一覧は `detect_arch.sh` から）。ライブラリはビルドツリーではなく `libs/` に出力します（`_loader.py` がそこを見るため）。
]


= 入力ファイルについて

== ファイル形式の選択（`ftype`）

ハミルトニアンの入力ファイル形式を`ftype`変数で指定します。`fname`にファイル名またはディレクトリ名（拡張子なし）を文字列で代入してください。

#table(
  columns: (auto, 1fr),
  [*ftype*], [*ファイル形式*],
  [0], [`{fname}/` ディレクトリ内に `ham_r.txt`, `irvec.txt`, `ndegen.txt` が存在する形式],
  [1], [`fname` で指定したファイル（独自形式。`.input`, `.hop` など拡張子も含めて指定）],
  [2], [`{fname}_hr.dat` ファイル（Wannier90のデフォルトホッピングファイル）。*第一原理計算との連携では通常これを使用*],
  [3], [ecalj の MLO 基底（非直交）のテキスト形式。$H(R)$ と重なり積分 $S(R)$ を含む],
  [4], [`{fname}/` ディレクトリ内の `HamRsMLO`, `HamiltonianPMTInfo`（ecalj `job_mlo` のバイナリ出力）。ftype=3 と同じく非直交 MLO 基底で, 重なり積分 $S(R)$ も同時に読み込みます。ftype=3,4 とも, 読み込み時に各軌道をノルム 1 に規格化します（バンドは変わりません）],
  [その他], [`Hopping.dat` ファイル（ecaljのホッピングファイル）],
)

`ftype=2` では Wannier90 の各行の軌道番号を $H_(m n)(R)$ の添字として読み,
縮退度で割ります。この模型の Fourier 規約は
$H(k) = sum_R H(R) exp(-2 pi i k dot R)$ です。

#detail[
  MLO の規格化: $D = "diag" S(R=0)$ として $H(R), S(R) -> D^(-1\/2) (dot) D^(-1\/2)$ とします。ecalj の MLO は規格化されておらず, Cu では $S_(a a)(0)$ が s で 0.016, d で 0.565 です。規格化しないと U は実効的に $U \/ S_(a a)(0)^2$ 程度で効いてしまうため, U, J はこの規格化した軌道の上の局所相互作用として扱います（FLEX, Eliashberg, AiiDA インターフェースでも同じ規約です）。
]

*Wannier90を使用する場合の例：*

```python
fname = 'inputs/Sr2RuO4'   # _hr.dat を除いたパスを指定
ftype = 2
```

== スピン軌道相互作用（`sw_soc`）

`sw_soc = True` にするとスピン軌道相互作用（SOC）を考慮します。ハミルトニアンはスピン自由度を含む（軌道数が 2 倍の）形式である必要があります。

- 対応: FLEX（option=14）と線形 Eliashberg（option=15）。IR 基底と MLO 基底にも対応しています。
- 非対応: 非線形 Eliashberg（option=23）と超伝導状態の $chi$（option=12, 13）。
- スピンの並びは既定で「全軌道の up の後に全軌道の down」です。別の並びの場合は `slist`（各軌道のスピン $plus.minus 1$）と `invs`（1 始まりの逆スピンの相手）を定義してください。

= ブラベー格子の設定（`brav`）

`brav`変数でブラベー格子の種類を指定します。これは逆格子ベクトルの計算やsymmetry lineの自動生成に使用されます。量子ESPRESSO（QE）やWannier90の出力と対応する設定が用意されています。

#table(
  columns: (auto, 1fr),
  [*brav*], [*格子の種類（対応するコード）*],
  [0], [単純格子（Simple cubic / tetragonal / orthorhombic）],
  [1], [面心立方格子（QEデフォルト, ibrav=2 相当）],
  [2], [体心立方格子（QEデフォルト, ibrav=3 相当）],
  [3], [六方晶格子（Hexagonal）],
  [4], [三方晶格子（Trigonal, QEの ibrav=5 相当）],
  [5], [底心格子（Base-centered）],
  [6], [面心立方格子（慣例的な向き）],
  [7], [体心立方格子（慣例的な向き）],
  [その他], [単斜晶格子（Monoclinic）],
)

= 計算モード（`option`）

計算内容は `option` 変数の整数値で指定します。`CalcMode` の名前で指定することもできます（例: `option = CalcMode.LIN_ELIASHBERG`）。

#table(
  columns: (auto, 1fr, auto),
  table.header([*option*], [*計算内容*], [*よく使う関連パラメータ*]),
  table.cell(colspan: 3)[_電子構造_],
  [0], [バンド図], [`k_sets`, `color_option`],
  [1], [状態密度（DOS）], [`Emin`, `Emax`, `delta`],
  [2, 3], [2 次元 / 3 次元フェルミ面], [`kz`, `kscale`, `color_option`],
  [4], [スペクトル関数 $A(bold(k), omega)$（FLEX の $Sigma$ 込みも可）], [`sw_self`],
  [17], [バンドフィリング（電子的 / ホール的の判定）], [`fill`],
  [18, 19, 20], [サイクロトロン質量, dHvA 振動数, 有効質量], [`RotMat`, `avec`],
  table.cell(colspan: 3)[_輸送・光学_],
  [5], [ボルツマン理論の伝導係数とホール係数], [`tau_mode`, `tau_const`],
  [6], [久保公式の光学伝導度, 反射率, 物質の色], [`delta`, `Emax`],
  [21, 22], [不純物スーパーセル / CPA のスペクトル], [`imp_*`],
  table.cell(colspan: 3)[_揺らぎ（RPA）_],
  [7, 8, 9], [スピン感受率（対称線 / 1 点 / $bold(q)$ マップ）], [`U`, `J`, `at_point`, `Ecut`],
  [10, 11], [ペアリング感受率（対称線 / $bold(q)$ マップ）], [`Ecut`],
  table.cell(colspan: 3)[_FLEX と超伝導_],
  [14], [FLEX 自己エネルギー], [`U`, `J`, `Nw` / `freq_grid`],
  [15], [線形 Eliashberg 方程式（$lambda$ とギャップ）], [`gap_sym`, `sw_self`],
  [16], [ギャップ関数の後処理], [],
  [23], [非線形 Eliashberg 方程式（$T_c$ 以下のギャップ）], [`sw_check_only`],
  [12, 13], [超伝導状態のスピン感受率], [`delta0`, `gap_sym`],
  [27], [超伝導状態の NMR（Knight shift, $1\/T_1 T$）], [`nmr_*`],
  table.cell(colspan: 3)[_準古典 Eilenberger_],
  [24, 25, 26], [一様系 / 表面 / 渦・渦格子], [`eil_*`],
)


== option=0: バンドプロット (`BAND`)

symmetry line に沿ったバンド分散 $E_n(bold(k))$ を計算し描画します。対称線は `k_sets` と `xlabel` で指定するか, 指定がなければ `brav` の設定から自動で選びます。

== option=1: 状態密度 (`DOS`)

$ "DOS"(omega) = - 1/pi sum_(bold(k), n) "Im" G^0_n (bold(k), omega + i delta) $

全軌道の合計と軌道ごとの部分 DOS をプロットします。

== option=2: 2 次元フェルミ面 (`FERMI_2D`)

指定した $k_z$ 断面（デフォルト: $k_z = 0$, 変数 `kz` で変更可能）における $k_x$-$k_y$ 面のフェルミ面を描画します。電子数 `fill` またはフェルミエネルギーから化学ポテンシャルを自己無撞着に求め, $E_n(bold(k)) = mu$ となる等エネルギー線を skimage の find\_contours で抽出します。`RotMat` で回転行列を指定することでフェルミ面を回転できます。

== option=3: 3 次元フェルミ面 (`FERMI_3D`)

Marching cubes アルゴリズムを用いて 3 次元フェルミ面をポリゴンで描画します。`kscale` で各軸方向の表示スケールを変えることができます。

既定（`kscale=1`）では, フェルミ面を第一ブリルアンゾーンに折り畳んで描きます。

#detail[
等エネルギー面を切り出す領域は分数座標で $k_("frac") in [-1\/2, 1\/2]^3$ の平行六面体, すなわち逆格子の*primitive cell*です。単純格子ではこれが第一ブリルアンゾーン（BZ 枠として描かれる Wigner–Seitz 胞）と一致しますが, 体心・面心・底心など複合格子や斜交格子では一致せず, そのままではフェルミ面が BZ からはみ出します。そのため `kscale=1`（既定値）のときは, ゾーンの外に出た部分を逆格子ベクトルで平行移動し, ゾーン面で切断して第一 BZ 内に折り畳みます（平行六面体の並進コピーは $bold(k)$ 空間を隙間なく埋めるので, 面積の重複も欠落も生じません）。`kscale` を 1 以外にした場合は意図的に別の $bold(k)$ 領域を見ていることになるので折り畳みは行いません。
]

`color_option=ColorMode.GAP`（3）にすると, Eilenberger 計算を駆動している*まさにそのペアリング形状因子* $"Re"[phi(bold(k))]$（`gap_sym`/`delta0`, または設定されていれば `eil_gap_orbital`/`eil_gap_file`）で表面を着色します。渦/表面ソルバーが使う固定 $k_z$ 断面ではなく, 実際の 3 次元 Wannier フェルミ面（全シート・全 $k_z$）上でギャップ（符号・ノード・異方性）を手早く確認できます。

== option=4: スペクトル関数 (`SPECTRUM`)

一体グリーン関数のトレースの虚部から電子のスペクトル関数

$ A(bold(k), omega) = -1/pi "Im" "Tr" [S(bold(k)) G(bold(k), omega + i delta)] $

を計算し, symmetry line に沿って描画します（直交基底では $S = 1$）。`sw_self=True` の場合は FLEX 計算が書き出した `self_en_wannier.npz` の自己エネルギーを $bold(k)$ 経路上へフーリエ補間し, Padé で実軸へ接続して, 相互作用の効果を含めたスペクトルを計算します。IR 基底で解いた $Sigma$ はサンプリング点の値（`sigma_ir`）から接続します（option=14 の「周波数格子」の節）。

== option=5: ボルツマン伝導係数 (`CONDUCTIVITY_BT`)

緩和時間近似のもとでボルツマン輸送方程式から輸送係数を計算します。緩和時間のモデルは `tau_mode` で選択します（既定は定数 $tau$）。

熱電関連の出力：
- 電気伝導率テンソル $sigma_(i j)$ （単位: S/m）
- 熱伝導率テンソル $kappa_(i j)$ （単位: W/m/K）
- ゼーベック係数テンソル $S_(i j)$ （単位: V/K）
- パワーファクター $sigma S^2$ （単位: $upright(W \/ m \/K^2)$）
- ペルチェ係数, ローレンツ数

#detail[
$S = -sigma^(-1) dots$ などに現れる逆行列は, 伝導のある方向だけで取ります。厳密な2次元モデルや $N_z = 1$（$k_z = 0$ だけを取るので $v_z = 0$）では $sigma_(z z) = 0$ になりますが, 面内の $S$, ペルチェ係数, $kappa$, ローレンツ数はそのまま求まり, 伝導のない方向の成分は 0 として出力されます（その旨が表示されます）。
]

=== 弱磁場ホール応答

同じモードで $B$ に線形なホール応答も出力します。ホールカーネルは *磁場3方向すべて* についてバンド分解して評価します。

$ S_(a b)^((g)) = 1/2 sum_(bold(k),n) w_bold(k) tau^2 (-partial f \/ partial epsilon) epsilon_(g u v) [v_a (M^(-1))_(b u) v_v - (a <-> b)] $

#detail[
ホール係数はテンソルのまま構成します。

$ rho^((1)) = -sigma^(-1) sigma^((1)) sigma^(-1), quad R_H^((g)) = rho^((1))_(b a) \/ B quad (a\,b\,g " は巡回") $

$sigma_(x y)\/(sigma_(x x) sigma_(y y))$ で代用すると $sigma$ が $x y z$ で対角であることを仮定することになり, 軸が直交しない格子では成立しません（単斜晶では $B=0$ でも $sigma_(x y) eq.not 0$）。$sigma$ が特異になる場合（厳密な2次元バンドでは $sigma_(z z)=0$）は, その軸に垂直な面内 $2 times 2$ ブロックで評価し, 残り2方向は未定義として出力します。
]

出力される物理量：
- $bold(B) parallel x, y, z$ それぞれのホール係数 $R_H$ （単位: $upright(m^3 \/ C)$）
- ホールキャリア密度 $n_H = 1\/(e |R_H|)$ （単位: $upright(c m^(-3))$）と, $"sign"(R_H)$ から判定した `electron` / `hole` の別
- ホール移動度 $mu_H = R_H sigma_(a a)$ （単位: $upright(c m^2 \/ V s)$）と 1 テスラあたりの $omega_c tau$
- $sigma^((1))_(x y)\/B$ （単位: $upright(S \/ m T)$）
- 本モードが前提とする弱磁場展開 $omega_c tau < 0.1$ が成立する磁場範囲

#block(fill: luma(240), inset: 8pt, radius: 4pt, width: 100%)[
  $n_H = 1\/(e|R_H|)$ が真のキャリア密度と一致するのは *単一の閉じたフェルミ面シートの場合のみ* です。多バンド金属では二流体平均 $R_H = (p mu_h^2 - n mu_e^2)\/(e(p mu_h + n mu_e)^2)$ となり, 補償近傍では発散・符号反転します。密度として読む前に, 併記される Luttinger 体積（後述）と必ず突き合わせてください。
]

=== $bold(k)$ メッシュ収束

$sigma^((1))$ は $M^(-1)$ を通じて $sigma$ より $bold(k)$ 微分が1階多く, さらに被積分関数がフェルミ面の周りで符号を変えるため, 収束が著しく遅くなります。収束していないメッシュでは $R_H$ が *符号ごと誤る* ことがあります。`hall_mesh_list` を設定すると一連のメッシュで計算して収束表を出力し, 最も細かい2点間で符号反転または5%超の変化があれば警告します。

#detail[
各メッシュについて2つの診断値を出力します。
- $N_"eff"$ — $-partial f\/partial epsilon$ 重みの参加率, すなわちフェルミ窓を担う実効的な状態数。全ての量が $N_"eff"$ 個の平均なので, これがノイズフロアを決めます（$10^3$ 未満で警告）。
- ホール被積分関数の相殺率 $|sum| \/ sum |dot|$。これが小さいときは大きな正負の寄与のわずかな残差になっており（補償近傍）, 符号がまだメッシュ依存の可能性があります（$10^(-2)$ 未満で警告）。
]

=== バンド分解

$sigma$ も $sigma^((1))$ も単純なバンド和なので, フェルミ面シートごとの内訳も出力します（各閉ポケットの Luttinger 体積から求めたキャリア密度, $sigma_(x x)$, そのシート単独の $R_H$ と $n_H$, 移動度）。各シート自身の $n_H$ を Luttinger 体積と比べれば, そのシートが単純な閉じたポケットとして振る舞っているかが分かり, 全体の $n_H$ と $n_e\/n_h$ の食い違いが実験で見えている補償そのものになります。

== option=6: 線形応答理論による伝導係数 (`CONDUCTIVITY_PT`)

久保公式（Kubo formula）に基づく線形応答理論を用いて電気伝導率などを計算します。`delta` パラメータが実効的な緩和時間 $tau approx hbar / delta$ に対応します。

さらに, 得られた光学伝導率 $sigma(omega)$ から以下の光学量を求めます（対角成分 $x x, y y, z z$）：

- 誘電関数 $epsilon(omega) = 1 + i sigma(omega) \/ (epsilon_0 omega)$
- プラズマ振動数（$f$-sum rule による非遮蔽値, および $"Re" epsilon_(x x) = 0$ の遮蔽値）
- 垂直入射の反射率 $R(omega) = |(N - 1) \/ (N + 1)|^2$, ここで $N = sqrt(epsilon)$ は複素屈折率
- 金属光沢の色: $R(omega)$ を CIE 1931 2° 標準観測者と標準光源 D65 で畳み込んで XYZ 三刺激値を求め, sRGB に変換した値（16 進表記・RGB 8bit・色度座標 $x y$・輝度 $Y$）。第一原理計算による金属色の標準的な手法に従っています（Prandini _et al._, npj Comput. Mater. *5*, 129 (2019)）。`poly` 行は 3 方向平均の $R$ による多結晶相当の色です。

非金属の場合は鏡面反射率が小さく（$n approx 1.5$ なら $R approx 0.04$）ほぼ無彩色なので, 色は物質を*通過した*光が担います。そこで吸収係数 $alpha = 2 omega kappa \/ (planck.reduce c)$ から次の 2 経路も出力します：

- *透過色*: 厚さ $d$ の板の透過率 $T = (1-R)^2 e^(-alpha d) \/ (1 - R^2 e^(-2 alpha d))$（両面フレネル + 非干渉多重反射）。$alpha = 0$ で $T = (1-R)\/(1+R)$ に帰着します。*この色は厚さ $d$ に依存する示量的な量*で, 金属光沢のような物質固有値ではありません。
- *体色（粉体・顔料）*: Kubelka--Munk 二流束模型の拡散反射率 $R_infinity = 1 + K\/S - sqrt((K\/S)^2 + 2 K\/S)$, $K\/S = alpha l$。実効散乱長 $l$ が唯一の自由パラメータです（$alpha$ と Kubelka--Munk の $K$ の規約差と粒径を吸収）。

さらに吸収率 $A = 1 - R - T$ からの色も出力しますが, $A$ は試料から出ていく光の補数なので*ネガ（`hex_neg`, sRGB で $1 - "rgb"$）を取って表示*します。これらは `d`, `l` を設定パラメータ `color_thick`, `color_scatt` で指定します。金属は不透明なので $T$, $R_infinity$ は黒く出ます（正しい挙動で, 図の色見本からも自動的に除かれます）。

出力は `optical.dat`（$sigma$, $epsilon$, $R$, $alpha$, $T$, $R_infinity$ のスペクトルとヘッダに色情報）, `dielectric.pdf/png`, `reflectivity.pdf/png`（色見本付き）です。可視域 380--780 nm (1.59--3.26 eV) を含むように `Emax` を 3.3 eV 以上に設定してください（不足する場合は警告が出ます）。

#detail[
色の精度について: 金属の色は滑らかな Drude 応答が支配するので誤差が色相に効きにくい一方, 非金属の色相は*吸収端の位置そのもの*です。模型のバンドギャップの誤差はそのまま色相のずれになり, 励起子効果とフォノン介在の間接遷移は本計算（Kubo バブル）に含まれません。`delta` のブロードニングも吸収端を鈍らせるため, 収束と模型の妥当性を確認して使用してください。
]

== option=7: スピン感受率スペクトル (`CHIS_SPECTRUM`)

RPA によるスピン感受率

$ chi_s (bold(q), omega) = frac(chi^0 (bold(q), omega), 1 - S chi^0 (bold(q), omega)) $

を symmetry line に沿って計算し, 虚部のスペクトルをプロットします（パラマグノンスペクトル）。相互作用行列 $S$ はオンサイト相互作用 `U`, `J` から生成されます。結果は `chis_spec.png` に出力されます。

== option=8: 特定 $bold(q)$ 点でのスピン感受率 (`CHIS_QPOINT`)

`at_point` で指定した波数 $bold(q)$ 点における $chi_s (bold(q), omega)$ の周波数依存性を計算・プロットします。

== option=9: $bold(q)$ 空間スピン感受率マップ (`CHIS_QMAP`)

`Ecut` で指定したエネルギー $omega_0$ における $chi_s (bold(q), omega_0)$ の $k_x$-$k_y$ 空間分布をプロットします。ネスティングベクトルの可視化などに使用します。出力は `chi0map.png` および `chismap.png`。

== option=10: ペアリング感受率スペクトル (`PHI_SPECTRUM`)

ペアリング感受率 $phi(bold(q), omega)$ を symmetry line に沿って計算しプロット（`phi_spec.png`）。

== option=11: $bold(q)$ 空間ペアリング感受率マップ (`PHI_QMAP`)

`Ecut` で指定したエネルギーにおけるペアリング感受率の $bold(q)$ 空間分布をプロットします（`phimap.png`）。`sw_omega` で実周波数 / 松原周波数を切り替えます。

== option=12: 超伝導状態のスピン感受率スペクトル (`CHIS_SPECTRUM_SC`)

非ゼロのギャップ関数 $Delta(bold(k))$ を仮定し, 異常グリーン関数 $F(bold(k))$ を含めた SC 状態の既約感受率

$ chi_0^"SC" = chi^"GG" plus.minus chi^"FF" $

を構築した上で RPA で動的スピン感受率 $chi_s^"SC" = chi_0^"SC" \/ (1 - S chi_0^"SC")$ を求め, symmetry line 沿いにスペクトルを描画します。`gap_sym` で初期対称性を選びます（負値で triplet）。出力は `chis_sc_spec.png`。

初期ギャップ振幅 `delta0` は, 全バンド共通の実数か, バンドごとの振幅・符号のリストで指定します（「エネルギー範囲・収束パラメータ」の `delta0`）。

MLO 入力の BdG 方程式では粒子・正孔の重なり $S(k)$ と正常ブロック
$H(k)-mu S(k)$ を用います。直交基底の場合と同様に, スピン非依存で時間反転対称な
正常状態を仮定します。モデルのバンドギャップは重なりを含めて軌道基底へ変換し,
NMR の Fermi 面上のギャップ尺度には対称直交化したギャップを使います。
バンド別振幅の変換には複素軌道の場合も一般化固有ベクトルを使います。

== option=13: 超伝導状態の特定 $bold(q)$ 点スピン感受率 (`CHIS_QPOINT_SC`)

option=12 と同じ SC 状態の枠組みで, `at_point` で指定した $bold(q)$ 点での $chi_s^"SC"(bold(q), omega)$ の周波数依存性を計算・プロットします。`delta0` の指定方法は option=12 と共通です。

実行時には以下も合わせて出力されます。

- $bold(k) = (0,0,0) arrow.r (0,0.5,0)$ 方向の BdG（Bogoliubov–de Gennes）バンド分散（`BdG_band.png`）
- $chi_s^"SC"(omega)$ の全軌道トレース（`chisq.png`, `chis_sc.dat`）
- 軌道分解した $chi_s^"SC"(omega)$（`chisq_orb.png`, `chis_scorb.dat`）

== option=14: FLEX 自己エネルギー計算 (`FLEX`)

FLEX 近似の自己エネルギー $Sigma(bold(k), i omega_n)$ を自己無撞着に求めます。得られた $Sigma$ は線形・非線形 Eliashberg（option=15, 23）, スペクトル（option=4）, 輸送（`tau_mode='flex'`）で使えます。

*主なパラメータ*

- `U`, `J`（または `orb_dep=True` の `Umat`, `Jmat`）: 相互作用
- `Nw`: 松原周波数の点数。または `freq_grid='ir'` で IR 基底を使います（下の「周波数格子」）
- `sw_out_self=True`: $Sigma$ をファイルに保存します
  - `self_en.npz`: $Sigma$ と周波数 `iw`, 温度。option=15, 23 の `sw_from_file` が読みます
  - `self_en_wannier.npz`: Wannier 実空間形式。option=4 と `tau_mode='flex'` が読みます
  - `sigma.bin`: 次の計算の初期値用（一様格子のみ）
- `sw_in_self=True`: `sigma.bin` を初期値として読み込みます（下の「継続計算」）
- `m_diis_num`: DIIS 加速の履歴長（既定 5）
- `sw_rescale_flex=True`: 反復の序盤に Stoner 因子が 1 に達しても発散しないようにします
- `sw_chi0_tail`（既定 True）: 一様格子での $chi_0$ の打ち切り誤差を減らします

#caution[
  標準出力に `[FLEX] WARNING: Stoner factor >= 1 at the final iteration` が出た場合, その温度・相互作用では正常状態が磁気的に不安定で, 得られた $Sigma$ は物理的な解ではありません。温度を上げるか `U` を下げてください。
]

*一様格子での `Nw` の選び方.* 一様な松原周波数の窓は $omega_c = (2 N_w - 1) pi T$ で, 温度に比例して狭くなります。最低温度でも $omega_c$ がバンド幅と `U` を十分上回るように `Nw` を選んでください。低温ほど大きな `Nw` が必要になるので, 低温や多軌道の計算では IR 基底（`freq_grid='ir'`）が便利です。

*非直交 MLO 基底（`ftype=3,4`）.* $G^(-1) = (i omega + mu) S - H - Sigma$ を使い, 粒子数は $"tr"(S G)$ で数えます（`U`, `J` の規格化は「入力ファイルについて」の節）。線形・非線形 Eliashberg とスペクトルも同じ規約です。自己エネルギーのファイルには重なり行列を保存しないので, 読み込むときは計算時と同じ模型を指定してください。`tau_mode='flex'` の輸送計算は直交基底専用です。

化学ポテンシャル探索にはバンド端の外側の熱的な裾も含めます。有限温度で完全に空または満杯の場合,
厳密には有限の化学ポテンシャルは存在せず, 戻り値は Fermi 関数が数値的に飽和する有限値です。

=== 継続計算（U・温度アニーリング）

収束しにくい条件では, 収束しやすい条件（小さい U や高い温度）で求めた `sigma.bin` を初期値にすると安定します。`sw_in_self=True` で前回の `sigma.bin` を読み込みます。

- *U アニーリング*: `sigma_in_scale` に $(U_"new"\/U_"old")^2$ 程度を指定します（弱結合で $Sigma prop U^2$）。
- *温度アニーリング*: 少しずつ冷やす（$T_"new"\/T_"old" gt.eq 0.8$ 程度）なら, そのまま読み込めば十分です。
- *`Nw` が前回と違う場合*: 読み込む前に `sigma.bin` が今回の `Nw` の格子へ自動で補間し直されます（その旨が表示されます）。

軌道数や $bold(k)$ 点数が前回と違う `sigma.bin` は使えず, エラーで止まります。

#detail[
  温度だけが違う `sigma.bin` をそのまま読むと, 周波数軸を $T_"old"\/T_"new"$ 倍に圧縮した初期値に相当します。この圧縮は冷却に伴う $Sigma$ の増大を近似的に模倣するので, 温度アニーリングの初期値として働きます。`Nw` が違う場合の自動補間も同じ圧縮を再現するように行うため, `Nw` を変えても変えなくても同じように振る舞います。
]

=== 精度と実装上の近似

#detail[
  - *`sw_rescale_flex`* は $chi_0(bold(q), i nu_n)$ を全運動量・全周波数で一律に $(1 - 10^(-4))\/S$ 倍します（$S = max_bold(q) "eig"[chi_0(bold(q), 0) hat(S)]$）。反復の序盤を乗り切るための処置で, 通常は $Sigma$ のフィードバックにより収束時には Stoner 因子が 1 未満になってリスケールは働かず, 収束解はリスケール無しの FLEX の固定点そのものです。収束時にも働いていれば上の警告が出ます。
  - *静的成分の減算*: 各反復で $Sigma(bold(k), i omega_0)$ の Hermite 部分を差し引き, バンドの静的なシフトを化学ポテンシャル側に吸収させてフェルミ面を安定させます。厳密な Hartree--Fock 項だけの減算ではないので, 得られた $Sigma$ から絶対的なバンドシフトは議論できません。
  - *松原方向の畳み込み* は長さ $2 N_w$ の FFT による循環畳み込みで, 打ち切り誤差は $O(1\/N_w)$ です。`Nw` を固定した温度走査ではこの誤差が温度とともに変わります。
  - *`sw_chi0_tail`* は $chi_0 = "conv"[G] - "conv"[G_0] + chi_0^"ref"$ と分け, ゆっくり減衰する参照部分（$G_0$ のバブル）を虚時間の解析式で厳密に扱い, FFT には速く減衰する残差（$G - G_0 tilde 1\/omega^2$）だけを通します。$chi_0$ の打ち切り誤差が $O(1\/N_w)$ から $O(1\/N_w^2)$ に改善します。$chi_0$ の段のコストは約 3 倍で, $N_w gt.eq 64$ 程度で有効です（極端に小さい $N_w$ では逆に悪くなり得ます）。SOC なしの経路（option=14, 15, 23）と MLO 基底に対応し, 改善は $chi_s$, $chi_c$, $V_sigma$, $V_Delta$ に伝わります。IR 基底では打ち切り誤差がないので無視されます。
  - $Sigma$ 自身の畳み込みの tail 誤差は主に周波数によらない定数（静的シフト）で, 静的成分の減算と化学ポテンシャルに吸収されます。Eliashberg カーネルの tail 誤差はギャップの一様（$s$ 波）成分にだけ入り, 符号反転するギャップでは対称性で消えます。このため, これらには個別の補正を入れていません。
]

=== 周波数格子: 一様な松原周波数と IR 基底（`freq_grid`）

`freq_grid='ir'` にすると, 周波数方向を一様な $N_w$ 点の松原周波数ではなく, IR（intermediate representation）基底の少数のサンプリング点で扱います。option=14, 15, 16, 23 の全て（SOC あり・なし, `sw_nosym` も）に対応し, sparse-ir が必要です。

```python
freq_grid = 'ir'
ir_wmax = None     # スペクトルの上限 [eV]。None で 2 x max|e_k - mu|
ir_eps = 1.0e-8    # 基底の精度。実軸へ Padé 接続するなら 1e-14 程度
```

- *点数*: 基底が決め, 典型的に 20〜30 点です（$Lambda = omega_"max"\/T$ と `ir_eps` に対数的にしか依存しません）。`Nw` と `sw_chi0_tail` は使われません。
- *メモリと時間*: 大きな配列（$G$, $Sigma$, $chi$ など）は周波数の点数に比例するので, 一様格子の `Nw` に比べて「`Nw` ÷ IR の点数」倍だけ減ります。例えば 2 軌道 FLEX（$64 times 64$, $T = 0.01$ eV）のピークメモリは, `Nw=2048` で 3.9 GB, IR（22 点）で 0.18 GB でした。計算時間も数十分の 1 程度になります。
- *精度*: 一様格子の $N_w -> infinity$ の極限とよく一致します（$mu$ と Stoner 因子で $10^(-5)$ 以下, 線形 Eliashberg の $lambda$ で $10^(-7)$ 程度）。
- *制限*: `sigma.bin`（`sw_in_self`, および `sw_out_self` のうち `sigma.bin` の部分）は一様格子専用です。
- *実軸への接続*: Padé を使う場合は `ir_eps=1e-14` 程度にしてください（点数は数点増えるだけです）。

#detail[
  - *`ir_wmax` の既定値*: RPA の $chi_0$ はバンド幅程度で足りますが, FLEX では $Sigma$ がスペクトルの重みをバンドの外へ広げるため 2 倍程度が必要です。
  - *ファイル*: `self_en.npz` にはサンプリング点上の $Sigma$ と周波数（`iw`）を保存し, 格子が異なる場合は読み込みを拒否します（`sw_from_file`, option=16）。`self_en_wannier.npz` は一様な $omega_n = (2n+1) pi T$ の先頭 64 点の形式で, IR 基底から評価して書き出します。ギャップ（`gap.npy`, `gap_wannier.npz`）はサンプリング点のまま保存し, 周波数を `gap_iw.npy` と `gap_wannier.npz` の `iw` に記録します。$omega_n -> 0$ への外挿は記録された周波数で行います。
  - *Padé の入力*: Padé は悪条件なので, IR の値を一様な $omega_n$ で評価し直すと, 基底の滑らかな打ち切り誤差が実軸の偽の構造（負のスペクトル, ピークの消失）として現れます。そこで Padé はサンプリング点の値そのもの（全点）から取ります（option=16 の `Gpade.dat`; `self_en_wannier.npz` にはサンプリング点の値 `sigma_ir`, `iw_ir` も保存し, option=4 の `sw_self` と `tau_mode='flex'` はこちらを使います）。精度は `ir_eps` で決まります。松原軸の上の量（$mu$, $lambda$, ギャップ, $chi$）にはこの問題はありません。
]


== option=15: 線形 Eliashberg 方程式 (`LIN_ELIASHBERG`)

RPA または FLEX で求めた有効ペアリング相互作用を使い, 線形化された Eliashberg 方程式を固有値問題として解きます。最大固有値 $lambda$ と対応するギャップ関数 $Delta(bold(k), i omega_n)$ が得られます（Arnoldi 法）。$lambda = 1$ で超伝導転移温度 $T_c$ となります。

- `sw_self=False`: 自己エネルギー補正なし（RPA）
- `sw_self=True`: FLEX の自己エネルギーを用いる（FLEX をその場で解きます）
- `sw_from_file=True`（`sw_self=True` と併用）: option=14 が書き出した `self_en.npz` から自己エネルギーを読み込み, FLEX 計算をスキップ（周波数格子と温度が一致しない場合は拒否されます）
- `gap_sym`: 初期ギャップ関数の対称性を指定（詳細は後述）
- *SOC あり*（`sw_soc=True`）では `gap_sym` の符号はギャップの空間反転パリティ（正: 偶, 負: 奇）を指定します（下の詳細）
- ギャップ関数は `gap.npy` と, 軌道ペアごとの `gap_11.dat`, `gap_12.dat`, …（軌道インデックスを 1 始まりにした `gap_{i j}.dat`）に出力されます（`sw_out_self=True` では `gap_iw.npy` と Wannier 実空間形式の `gap_wannier.npz` も）

#detail[
  SOC があると `gap_sym` はスピン構造を指定しません。SOC があるとパリティセクターの中でスピン singlet と triplet が混ざるためです。同じセクターには指定した状態の相手（異スピン・奇周波数）も同居するので, ソルバーはセクターの最大固有値ではなく `gap_sym` のスピン性格に一致するブランチを返します（SOC $arrow.r 0$ で $lambda$ が連続につながります）。SOC が無い場合は singlet/triplet が厳密なセクターなので, そのまま射影されます。
]

=== 対称性を仮定しない線形 Eliashberg（`sw_nosym`）

`sw_nosym=True` にすると, option=15 は時間反転・空間反転対称性を仮定せずに全ブリルアンゾーン上でギャップ方程式を解きます（裸の $G$, RPA）。Zeeman 場や有限運動量の対（FFLO）のように, それらの対称性が破れる状況のためのものです。スピンを持たない模型（`sw_soc=False`）は自動的にスピン 2 重化されます。singlet/triplet への射影は行わず, `gap_sym` は初期値にだけ使われます。FLEX の自己エネルギーには対応していません（`sw_self` を指定しても裸の $G$ を使います）。

- `zeeman_h`（3 成分, eV）: Zeeman 場 $H_Z = bold(h) dot bold(sigma)$（$g mu_B B\/2$ を含めた値）。化学ポテンシャルは同じ `fill` で解き直します。
- `fflo_q`（整数 3 成分のリスト）: 対の運動量 $bold(Q) = (i_(Q x)\/N_x, i_(Q y)\/N_y, i_(Q z)\/N_z)$ の一覧。例えば `[[q, 0, 0] for q in range(5)]` で $k_x$ 方向に走査します。$chi_0$ と $V_Delta$ は一度だけ計算して全ての $bold(Q)$ で使い回し, $lambda(bold(Q))$ を `lambda_Q.dat` に出力します。ギャップは `gap_nosym.npz`, `gap_nosym.dat` に出力されます。
- `nosym_arnoldi_m`（整数, 既定 60）: Arnoldi の Krylov 空間の次元。$lambda$ が小さい場合, $10^(-6)$ の精度には 60 程度以上が必要なことがあります。

== option=16: ギャップ関数の後処理 (`GAP_FUNCTION`)

option=15 で計算したギャップ関数 (`gap.npy`) を読み込み, 以下の量を計算・出力します。

- ギャップ関数の異常グリーン関数 $F(bold(k), i omega_n)$
- Padé 近似による実周波数への解析接続スペクトル（`freq_grid='ir'` ではサンプリング点から接続します。`ir_eps=1e-14` 程度を推奨。option=14 の「周波数格子」の節を参照）

== option=17: バンドフィリング (`BAND_FILLING`)

バンドの占有状況を手早く確認するための診断用モードです。各バンドについて, $bold(k)$ メッシュのうち $mu$ より上・下にある割合を出力します。非占有側がそのバンドのホール的な側, 占有側が電子的な側なので, この2つでシートがバンドの底寄りか天井寄りかが一目で分かります。占有割合の総和は単位胞あたりの電子数になり, `fill` との整合性チェックになります。 この2つの割合はどちらも Luttinger 体積であり, $upright(c m^3)$ あたりではなく単位胞・スピンあたりで規格化したものです。フェルミ面を持たないバンドは 1（単位胞に電子1個）または 0 になり, ポケットを持つバンドは対応する側にそのポケットの $V_"occ"\/V_"BZ"$ が出ます。option=5 で使う `fs_carrier_density` は同じ体積に $g_s\/V_"uc"$ を掛けたものです。`fill` と突き合わせるにはこちらの規格化が, 実験のキャリア濃度と比べるには $upright(c m^(-3))$ の方が対応します。

== option=18: サイクロトロン質量 (`CYCLOTRON_MASS`)

オンサーガー関係 $m^*_c = (hbar^2 \/ 2 pi)(partial S \/ partial E)$ に基づきサイクロトロン質量を計算します（$S$ はフェルミ面断面積）。磁場方向は回転行列 `RotMat` で指定します（実行時に対応する $hat(B)$ と $theta$, $phi$ が出力されます）。`theta`/`phi` を直接走査する option=19 とは指定方法が異なる点に注意してください。

- Phase 1: $k_z in [0, 0.5]$（簡約座標）を `meshkz=20` 点で走査し $S(k_z)$ をマップ
- Phase 2: 極値 $k_z$ で central finite difference により $partial S \/ partial E$ を計算し $m^*_c$ を $m_e$ 単位で出力

skimage の find\_contours と marching\_cubes を内部で使用します。

== option=19: dHvA 振動数 vs 角度 (`DHVA`)

閉じたフェルミ面断面の極値から, オンサーガー関係 $F = (hbar \/ (2 pi e)) A_"ext"$ により幾何学的な dHvA 振動数を計算します。$A_"ext"$ は磁場に垂直な面上で極値をとる断面積です。`theta` は `avec` を与えた座標系の Cartesian $z$ 軸から測る極角, `phi` は $x y$ 面内の方位角です。通常の設定のように $c parallel z$ であれば, $theta=0$ が $B parallel c$ に対応します。通常の `main.py` の実行経路では `phi=0` とし, $theta=0$–$90$ 度を等間隔の 40 点で走査します。

このルーチンはフェルミ面を Cartesian な逆空間で扱い, 指定した物理的な磁場方向に垂直な平面を掃引します。各平面上の閉曲線を検出して面の変位とともに追跡し, 連続する各軌道枝について最小・最大断面を抽出します。逆格子ベクトルだけ異なるゾーンコピーは統合します。このため非立方格子や傾いた磁場にも使えますが, `avec` には物理的な原始格子ベクトルを angstrom 単位で与える必要があります。断面積は $angstrom^(-2)$, 振動数は tesla 単位で出力されます。

結果は `dhva_band.png` に保存され, 横軸は `theta`, 縦軸は tesla 単位の振動数です。端末にはフェルミ面を持つバンド, 掃引窓, 軌道枝の数, 抽出した極値軌道数も出力されます。この option が求めるのは振動数のみです。振動振幅, Dingle 因子, Zeeman 分裂, サイクロトロン質量は計算しません（質量には option=18 を使用）。

開いた曲線および断面の掃引窓を広げたときに面積が変わる曲線は, 閉じた dHvA 軌道を定義しないため棄却します。棄却された曲線についての警告は直ちにエラーを意味しませんが, 振動数を使う前には面内メッシュ `Nx`, 変位平面数 `meshkz`, 補間格子 `grid_mesh`（既定値 120）を増やして収束を確認してください。特に $B parallel a b$ の近傍, Lifshitz 転移の近傍, 小さなポケットでは重要です。

計算した振動数は, 同じ種類のバンド計算の結果と比べてください。DFT のバンドは電子相関を含まないため, 実験の振動数とは一般に一致しません。実験を再現するには, 相関を取り込んだバンド構造（QSGW, LDA+DMFT）や経験的なバンドシフトに基づく Wannier 模型が必要です。

== option=20: 電子質量 (`ELECTRON_MASS`)

symmetry line 上で $m^* = hbar^2 (partial^2 E \/ partial bold(k)^2)^(-1)$ を計算します（$m_e$ 単位）。

== option=21: 不純物入りスペクトル (`SPECTRUM_IMPURITY`)

周期的なスーパーセル上にオンサイトの不純物ポテンシャルを置いたハミルトニアンを対角化し, 元の単位胞の $bold(k)$ へ unfold したスペクトル関数 $A(bold(k), omega)$ を計算します。

- `imp_cell`（3 成分）: 各軸方向のスーパーセルの大きさ。ホッピングの到達距離の 2 倍より大きくしてください（$N_"cell" times N_"orb"$ 次元の行列を対角化します）。
- `imp_sites`（`None` または整数のリスト）: スーパーセル内の不純物サイト（0 始まり, `gen_rlist` の順）。`None` では `imp_conc` の濃度でランダムに置きます（少なくとも 1 個）。
- `imp_conc`, `imp_seed`: ランダム配置の濃度と乱数の種（`None` で毎回異なる配置）。
- `imp_V`（実数 または $N_"orb" times N_"orb"$ 行列, eV）: オンサイトの不純物ポテンシャル（実数は単位行列の倍）。

MLO 入力では $H(R)$ と $S(R)$ の両方をスーパーセルへ折りたたみ, 一般化固有値問題を解きます。
unfolding にはスーパーセル重なりの対称 Löwdin 直交化を用い, 軌道分解重みは直交化後の軌道に対応します。
`imp_V` は $S$ ではなく $H$ に加える局所行列摂動です。

== option=22: CPA 伝導 / スペクトル (`SIGMA_CPA`)

Coherent Potential Approximation により混晶の自己エネルギーを反復解で求め, 実周波数のスペクトル関数（`cpa_spectrum.png`）と松原周波数の自己エネルギー（`sigma.bin`, `self_en.npz`）を出力します。濃度とオンサイト準位のずれは `main.py` の `calc_cpa` 呼び出しで固定値（`x_cpa=0.5`, `VB_shift=1.0` eV）を渡しています。変える場合はその行を編集してください。

`x_cpa` は準位がずれた成分の濃度です。MLO 入力の CPA 媒質は
$(z+mu) S(k)-H(k)-Sigma(z)$ を用い, スペクトルは $-"Im" "tr"(S G) \/ pi$ で評価します。
非直交 CPA の反復には NumPy 行列演算を使うため, 直交基底のネイティブ実装より遅くなる場合があります。

== option=23: 非線形 Eliashberg 方程式 (`NONLIN_ELIASHBERG`)

非線形（自己無撞着）SC FLEX-Eliashberg ループを解きます。線形 Eliashberg と異なり Δ を $T_c$ 以下で有限振幅まで成長させ, 異常グリーン関数 $F$ を含めた SC Dyson 方程式と FLEX 自己エネルギーを反復します。

初期ギャップは内部で自動生成されます（`gap.npy` を事前に用意する必要はありません）。

1. まず線形化 Eliashberg 方程式を内部で解き, Stoner 因子 $S$ と最大固有値 $lambda_"eliash"$ を計算します。
   - $S >= 1$（SDW/CDW 不安定）の場合は非線形ループに入らず終了します。
   - $lambda_"eliash" < 1$（$T >= T_c$, 超伝導不安定性なし）の場合も終了します。
   - `sw_check_only=True` の場合はここで終了し, $S$ と $lambda_"eliash"$ のみを報告します（非線形ループをスキップして温度走査で $T_c$ を素早く確認する用途）。
2. 上記の線形固有ベクトルから対称性が正しい初期 Δ の形を取り, 振幅を BCS 弱結合比 $Delta_0 = 1.764 k_upright(B) T_c$ にスケールして非線形ループに渡します。
#detail[
  各反復: $Delta_"new" = T \/ N_k sum V_Delta dot F$ → 振幅方向 Newton（secant）加速 + DIIS 形状混合（バイパス時は線形混合 pp=0.3）→ $Sigma$ 更新（`sw_self=True`; 静的部分を差し引き, `pp` で混合）→ Nambu の Dyson 方程式 $G = [G_0^(-1) + Delta G_0^dagger Delta^dagger]^(-1)$, $F = -G_0 Delta G^dagger$ で $G, F$ 更新 → $chi^"GG", chi^"FF"$ から $V_sigma, V_Delta$ 再構築（$V$ を固定する設定ではこの段を省略）。]

- `sw_self=True`: FLEX 自己エネルギー込みで Σ-dressed Green's function を使う
- `sw_from_file=True`: option=14 が書き出した `self_en.npz` から自己エネルギーを読み込む
- `sw_out_self=True`: ギャップを `gap.npy`, `gap_wannier.npz` に, `sw_self=True` では超伝導状態の $Sigma$ を `self_en_sc.npz`, `self_en_sc_wannier.npz` に出力します
- 反復が発散した場合（$Delta$ が NaN/Inf）はエラーで停止し, ギャップを書き出しません。`itemax` 回で収束しなかった場合は警告を表示します
- IR 基底（`freq_grid='ir'`）に対応しています

== 準古典 Eilenberger モード（option=24,25,26）

option=24〜26 は超伝導の*準古典 Eilenberger 方程式*（エネルギー積分された Gor'kov 方程式）を解きます。ギャップや不純物がフェルミ波長スケールでゆっくり変化する（$Delta, T_c, hbar\/tau << E_F$）状況で有効で, $T_c$・状態密度・表面 Andreev 束縛状態・渦／渦格子状態を扱う自然な枠組みです。結合定数 `eil_coupling` は無次元の $lambda$, ギャップの形は `gap_sym` で決まります。

#detail[
  未知数は準古典伝播関数 $g(bold(k)_F, bold(r), omega_n)$（normal）と $f$（anomalous）で, 単一の*Riccati 振幅* $a$（$f = 2a\/(1+a a^*)$, $g = (1-a a^*)\/(1+a a^*)$）で表し, フェルミ速度に沿った直線軌道上を数値的に安定な Fortran カーネル（`riccati_chords`; d-vector 用には $2 times 2$ スピン版 `matrix_riccati_chords`）で積分します。ペアリングは可分離形 $Delta(bold(k)_F, bold(r)) = Delta(bold(r)) phi(bold(k)_F)$ で, 形状因子 $phi$ は `gap_sym` で決まり $ |phi|^2 _"FS" = 1$ に規格化されます。
]

フェルミ面は 3 モード共通で `eil_fs_kind` により選びます: `None` は等方シリンダー（解析的な角度平均）, `'iso'`/`'ellipse'`/`'tb'` は `eil_fs_params` から作るモデル FS, `'wannier'` は読み込んだ Wannier バンドの実 FS とフェルミ速度（ギャップ対称性／マルチバンド構造は `gap_sym`, `delta0`, `eil_gap_orbital` から）。温度はグローバルの `tempK`/`temp`, ギャップ対称性はグローバルの `gap_sym` を用います（モデル FS では整数インデックスを連続調和関数に写像し, $2$（$s^plus.minus$）$-> s$）。

=== option=24: 一様 Eilenberger (`EILENBERGER`)

バルク（空間一様）ソルバです。`eil_*` サブモードフラグが全て off ならギャップ $Delta(T)$ を自己無撞着に解いて $T_c$ を報告し, `eil_find_tc=True` で二分法により $T_c$ を求めます。サブモードフラグ（排他的）は以下を選びます:

- `eil_imp_sweep=True`: 非磁性不純物率 $Gamma$（`eil_imp_list`）を走査し $T_c(Gamma)$ を `eilenberger_tc.dat` に出力 — Abrikosov–Gor'kov 対破壊曲線（等方 $s$ 波は Anderson の定理により抑制なし, 符号反転ギャップは強く抑制）。`eil_imp_c` は Born（$-> infinity$）〜ユニタリ（$-> 0$）散乱を内挿。
- `eil_pauli=True`: Zeeman（Maki）Pauli 限界走査 — シングレットギャップ $Delta(h)$, 一次転移のスピノーダル（Chandrasekhar–Clogston）, Zeeman 分裂した DOS。
- `eil_spin=True`: スピン分解（$2 times 2$）Zeeman 応答 — シングレット／平行 d-vector（$bold(d) parallel bold(h)$, Pauli 限界）vs 垂直 d-vector（$bold(d) perp bold(h)$, Zeeman 不感）。
- `eil_lambda=True`: 超流動密度 $rho_s(T)$ と侵入長 $lambda(T)$（フルギャップは指数的に平坦, ノードギャップは $T$ 線形）→ `penetration_depth.dat`。
- `eil_fs=True`: フェルミ速度付きモデル FS で同様, 異方的な $lambda_(x x), lambda_(y y)$ → `fs_penetration.dat`。
- `eil_free_energy=True`: 凝縮自由エネルギー $(Omega_s - Omega_n)\/N_0$ の $T$ 依存（結合定数積分）→ `free_energy.dat`。

=== option=25: 表面 Andreev 束縛状態 (`EILENBERGER_SURFACE`)

鏡面反射軌道に沿った Riccati 積分で鏡面表面近傍の自己無撞着ギャップ profile $Delta(x)$ を解き,（`eil_ldos=True` で）表面 LDOS も計算します。表面方位は `eil_surf_beta`（$d$ 波では $0 = [100]$ は束縛状態なし, $pi\/4 = [110]$ は反射時に感じる符号反転からゼロエネルギー Andreev 束縛状態 ZEBS が出現）。Zeeman 場 `eil_zeeman` は ZEBS を $plus.minus h$ に分裂させます。`eil_surf_dvector=True` では代わりに表面の自己無撞着トリプレット*d-vector テクスチャ*（スピン行列 Riccati による主成分＋副成分, 結合比 `eil_dvec_subratio`）を解きます。

表面・渦モードで `eil_fs_kind='wannier'` を使う場合, `mu0` を指定するとその値で
フェルミ面を構築します。未指定の場合は `fill` から化学ポテンシャルを計算します。

=== option=26: 渦・渦格子 (`EILENBERGER_VORTEX`)

磁束渦まわりの非一様ソルバです。`eil_field` $= B\/H_(c 2)$ で形状を選びます: `0` は大きな円形セル内の孤立渦（半径 `eil_vort_lxi`（$xi$ 単位）, グリッド `eil_vort_ngrid`）, `>0` は Doppler シフト付き円形セル渦格子。ギャップ profile $Delta(rho)$ と（`eil_ldos` で）ゼロエネルギーコア LDOS（Caroli–de Gennes–Matricon 束縛状態; 格子では Volovik の $sqrt(B)$ DOS）を出力します。Zeeman 場 `eil_zeeman` はコア状態をスピン分裂させます。サブモード:

- `eil_vort_current=True`: 環状超流動電流 $j_phi(rho)$ → `vortex_current.dat`。
- `eil_vort_field=True` / `eil_vort_maxwell=True`: 自己無撞着な有限 $kappa$ 磁場 $B(rho)$ / ベクトルポテンシャル $bold(A)(bold(r))$（Maxwell 反作用; `eil_kappa` を使用）。
- `eil_vort_lattice_sc=True` ＋ `eil_field_list`: *真の周期*磁気 Bloch 渦格子（定式化 A, 極端第 II 種）: 各コアに実ノードを持つ複素秩序変数 $Psi(bold(r))$ と完全な Abrikosov 超流動位相, $B\/H_(c 2)$ を走査して $ N(0)  (B)$ を得る（$d$ 波 $tilde sqrt(B)$ Volovik）。`eil_lattice` は `'square'`/`'triangular'`, `eil_nvortex` はセルあたり磁束量子数, 有限 `eil_kappa` で London 遮蔽を追加（`eil_vort_scA=True` で $bold(A)$ を準古典電流 $bold(j)_s =  bold(v)_F "Im" g $ から完全自己無撞着に, `je` の `A_renew` 法）。
- `eil_vort_dvector=True`: 自己無撞着トリプレット d-vector 渦／格子テクスチャ（主成分の winding ＋コア局在副成分; スピン行列 Riccati）。
- `eil_gap_orbital` / `eil_gap_file`: 軌道基底のギャップを FS バンドに射影して形状因子に使います（「Eilenberger パラメータ」の節）。


== option=27: 超伝導状態の NMR (`NMR_SC`)

option=12,13 と同じ BdG スピン帯磁率から, 超伝導状態の代表的な NMR 観測量である Knight shift とスピン格子緩和率を計算します。

$ K(T) prop chi'_s (bold(q)=0, omega -> 0), quad 1\/(T_1 T) prop angle.l |A(bold(q))|^2 chi''_s (bold(q), omega_0) angle.r_bold(q) \/ omega_0 . $

同じ温度・$bold(k)$ メッシュ・$bold(q)$ メッシュ・スメアリングで *正常状態* も評価され, 物理的に意味を持つのはその比です。$K_s\/K_n$ はシングレットギャップに対する Yosida 関数, $(1\/T_1T)_s\/(1\/T_1T)_n$ はフルギャップの Hebel--Slichter ピーク, あるいは節を持つギャップのべき則（ライン節なら $T^3$）を与えます。`nmr_sc.dat` の生の `K_sc`, `K_n`, `(1/T1T)` 列は診断用と考え, 比を引用してください。

#detail[
  比を取ると物質固有の係数（$gamma_n$, 超微細結合のスケール, $g$）が消え, メッシュ誤差の大部分も相殺するため, 比は各々の生の値よりはるかに速く収束します。正常状態の参照値は正常状態のバブルで計算します（SC バブルの厳密な $Delta -> 0$ 極限であり, 4 倍高速です）。
]

ギャップは `gap_sym`/`delta0` の形状因子（option=12,13 と同じもの）か, `nmr_gap_file` を指定した場合は option=15/23 が出力した自己無撞着な RPA/FLEX/Eliashberg ギャップです。いずれも $Delta(0)$ として読まれ, 各温度では BCS 内挿 $Delta(T) = Delta(0) tanh(1.74 sqrt(T_c\/T - 1))$ によって *振幅だけ* がスケールされ, 形状は固定されます。$T_c$ の既定値は *フェルミ面上* の $|Delta|$ の最大値から得る弱結合の見積り $max |Delta|_"FS" \/ 1.764$ です（ゾーン全体の最大値は, 状態を担うシートに形状因子の重みが乏しい場合に何倍も大きくなり, 掃引全体が誤った換算温度になってしまいます）。`nmr_tc` で上書きできます。

掃引は `nmr_trange`（$T_c$ 単位）を `nmr_nt` 点で行い, `nmr_sc.dat` と `nmr_sc.png` を出力します。

#detail[
  開始前に, 使用するスピンチャンネル, $Delta(0)$ と $T_c$, $bold(k)$ メッシュとスメアリングの階層のチェック（$delta gt.tilde v_F d k$ が必要）, コスト見積り（総コストは $N_T N_q N_k$ に比例）, および `nmr_qconv=True` のときは `nmr_qsize` とその半分のメッシュでの $bold(q)$ 和の変化が出力されます。
]

= カラープロット設定（`color_option`）

option=0,2,3では, バンドやフェルミ面の各点に色をつけることができます。

#table(
  columns: (auto, 1fr),
  [*color_option*], [*意味*],
  [0], [色なし（黒単色）],
  [1], [`olist`で指定した軌道成分をRGB（赤/緑/青）の混合比で表示],
  [2], [群速度 $|bold(v)(bold(k))|$ の大きさを色で表示],
  [3], [（option=3, `FERMI_3D` 限定）Eilenberger ペアリングギャップ $"Re"[phi(bold(k))]$ — `gap_sym`/`delta0`, または設定されていれば `eil_gap_orbital`/`eil_gap_file` の Nagai–Nakamura 射影],
)

`color_option=1`の場合, `olist`で軌道インデックスを `[R成分, G成分, B成分]` の形で指定します。複数の軌道を同じ色に割り当てる場合は, リストの要素をさらにリストにします。例:

```python
olist = [[0, 3], [1, 4], [2, 5]]
# 軌道0と3をRed, 軌道1と4をGreen, 軌道2と5をBlueに対応
```

= 各種パラメータの詳細説明

以下では, `main.py`の上部に記載されている全パラメータについて物理的な意味を含めて説明します。

== 基本設定

- `fname` （文字列）: 入力ファイルまたはディレクトリのパス（`ftype` に応じた形式で指定。ftype=1 は拡張子込みのファイル名, ftype=2 は `_hr.dat` を除いた名前, ftype=0, 4 はディレクトリ）

- `ftype` （整数）: 入力ハミルトニアンのファイル形式（「入力ファイルについて」の節を参照）

- `brav` （整数）: ブラベー格子の種類（「ブラベー格子の設定」の節を参照）

- `sw_soc` （bool）: スピン軌道相互作用を考慮するかどうかのスイッチ

- `inv_tol` （実数, 既定 $3 times 10^(-2)$）: 反転操作を $H(bold(R))$ の対称性として受け入れる相対残差のしきい値（option=15,16,23 は $Delta(bold(k))$ から $Delta(-bold(k))$ を再構成するのにこの操作を使います）。第一原理計算のフィットモデルには数 % の誤差があるため既定値は意図的に緩めてあります。フィットモデルが棄却される場合は大きくしてください。大きな値が必要だったときは, 出力される残差を確認してから結果を使ってください。

== メッシュ設定

- `Nx, Ny, Nz` （整数）: 第一ブリルアンゾーン内の $bold(k)$ 点メッシュ数。フェルミ面積分や FFT-convolution（FLEX, Eliashberg, chi）で使用します。2 次元系では `Nz=1` とします。FLEX/Eliashberg 計算では FFT 効率の都合上 2 のべき乗（32, 64 など）が望ましく, メモリ消費は $tilde.equiv N_x N_y N_z dot N_w dot N_"orb"^2$ で増えます。

- `Nw` （整数）: 松原周波数のメッシュ数（option=12,13,14,15,16,23 の FLEX/Eliashberg/SC chi 計算; `freq_grid='ir'` では使われません）または実周波数のメッシュ数（option=1,4 など）。温度 $T$ に対して $omega_n = (2n+1) pi T$ の形で $n=0,1,...,N_w-1$ まで取ります。低温ほど大きな $N_w$ が必要です（典型値: 256–1024。選び方は option=14 の節）。

- `kmesh` （整数）: バンド・スペクトルをsymmetry lineに沿って計算する際の$bold(k)$点数。大きいほど滑らかなプロットが得られます（例: 200～500程度）。

== 格子定数

- `abc` （リスト, 単位: Å）: 格子定数$a, b, c$。速度 $v = 1/hbar (partial E)/( partial bold(k))$ の計算やsymmetry lineの幅の計算（実空間の長さ換算）に使用します。Wannier90の場合はWIN/WOUTファイルの値に合わせてください。

- `alpha_beta_gamma` （リスト, 単位: degree）: 格子角度$alpha, beta, gamma$。直方晶・立方晶の場合は `[90., 90., 90.]` です。

== 温度・化学ポテンシャル

- `tempK` （実数, 単位: K）: 温度をケルビンで指定します。コード内で $T = k_B$`tempK`（単位: eV）に変換されます。

- `temp` （実数, 単位: eV）: $k_B T$を直接eV単位で指定する場合に使います。`tempK`と同時に定義した場合は`temp`が優先されます。

- `fill` （実数）: バンドフィリング。化学ポテンシャル $mu$ を $sum_(bold(k), n) f(epsilon_n (bold(k)) - mu) = N_k dot f"ill"$ で求めるための目標占有数で, 0 から軌道数 `no` までの値を取ります。`fill = no` で full-filled。物理的な 1 ユニットセルあたりの電子数との対応はハミルトニアンの形式に依存します:
  - SOC なし (`sw_soc=False`): スピン自由度は陽に含まれないので, `fill` はスピンあたりの電子数に対応。たとえば 3 軌道モデルで `fill=3` は full-filled（全スピン込みで 6 電子相当）, `fill=1.5` が half-filled。
  - SOC あり (`sw_soc=True`): ハミルトニアンが既にスピンを含む（`no = 2 N_"orb"`）ので, `fill` は 1 セルあたりの全電子数そのもの。3 軌道（`no=6`）で `fill=6` が full-filled, `fill=3` が half-filled。

- `mu0` （実数, 単位: eV）: 化学ポテンシャルを直接指定したい場合に定義します。定義されている場合, `fill`からの自己無撞着計算をスキップして`mu0`の値を使用します。

== エネルギー範囲・収束パラメータ

- `Emin, Emax` （実数, 単位: eV）: DOSやスペクトルを計算する際のエネルギー範囲の下限・上限。

- `delta` （実数, 単位: eV）: スペクトル計算のブロードニングパラメータ $delta$。グリーン関数の虚部に加える小量で, ディラックのデルタ関数を有限幅のローレンツ型に広げます。大きすぎると特徴が潰れ, 小さすぎると数値ノイズが出ます（0.01〜0.05 eV程度が目安）。

- `Ecut` （実数, 単位: eV）: option=9,11 で $bold(q)$ 空間の感受率マップを作成する際の固定エネルギー $omega_0$（フェルミ面近傍は 0 付近に設定）。

- `delta0` （実数 または 長さ`Norb`のリスト, 単位: eV）: option=12,13,27 の SC chi 計算で用いる初期ギャップ関数の振幅。物理的には超伝導ギャップサイズに対応します（典型値: $10^(-3) - 10^(-2)$ eV ≒ 1–10 meV）。
  - 実数を指定した場合: 全バンド共通の単一ギャップとして, 内部で生成した対称性形状を最大振幅 `delta0` にスケールします。0 にすると常状態と等価になります。
  - リストを指定した場合（例: `delta0=[0.,0.2,0.3,-0.1,0.]`）: バンド毎に振幅と符号を個別に指定する multi-gap モードになります。符号が異なる成分を混在させることで $s^plus.minus$ 型など符号反転ギャップを表現できます。

== 伝導計算パラメータ

- `tau_const` （実数, 単位: fs）: `tau_mode='const'` のときに使う定数緩和時間 $tau$（`'e2'`, `'dos'` ではフェルミ準位での $tau$）。`'epa'`, `'flex'` では使われないので正の値である必要はなく, ライブラリ関数に直接渡す場合は `None` でも構いません。定数近似では自由パラメータであり, 実験値との比較から決めます（例: 金属では 1〜100 fs 程度）。なお $R_H$ は $tau^2\/tau^2$ を持つため定数 $tau$ には *依存しません*。`tau_mode='const'` では, ホール係数はフェルミ窓の広がり以外に温度依存性を持たない純粋なバンド構造由来の量になります。

- `tau_mode` （文字列, 既定 `'const'`）: option=5 の緩和時間モデル。
  - `'const'`: 定数 `tau_const`（$bold(k)$・バンドに依らない）。
  - `'epa'`: 電子フォノン平均結合（Samsonidze--Kozinsky の EPA）による $tau(bold(k))$。`epa_file` が必要です。$R_H (T)$ やシートごとの移動度差が意味を持つのはこのモードです。
  - `'e2'`: フェルミ液体型 $tau$。$1\/tau prop (epsilon-mu)^2 + (pi k_B T)^2$（電子間散乱）で, フェルミ準位で `tau_const` になるよう規格化。
  - `'dos'`: $tau prop 1\/D(epsilon)$（散乱率が終状態の DOS に比例するフェルミの黄金律）。フェルミ準位との比として使い, $tau(epsilon) = tau_"const" D(mu)\/D(epsilon)$。
  - `'flex'`: FLEX 自己エネルギーによる準粒子 Boltzmann。同じ温度で `sw_out_self=True` として FLEX を解いた `self_en_wannier.npz` が必要です。化学ポテンシャルは FLEX のもの（`fill` は使わない）で, `tau_const` は無視されます。電流バーテックス補正は含みません（バブル近似）。
    #detail[軌道表示の行列 $Sigma(bold(k), i omega_n)$ を *FLEX を解いた k メッシュ上で* Padé により実軸へ接続し（補間した値ではなくソルバーのデータそのものなので, 補間ノイズによる偽の極が出ません）, 実エネルギー格子上の $Sigma^R (bold(k), epsilon)$ を輸送の k メッシュへ多重線形補間してバンドに射影します（凸結合なので因果性が保たれます）。そこから $Z = [1 - partial "Re" Sigma^R \/ partial epsilon]^(-1)$, $epsilon^* = Z(xi + "Re" Sigma^R (0))$, $Gamma = -"Im" Sigma^R (epsilon^*)$ から, エネルギー $epsilon^*$, 速度 $Z v$, 寿命 $tau^* = hbar \/ (2 Z Gamma)$ を使います。$Gamma$ を $epsilon^*$ で評価するので, ゼーベック係数に効く散乱の粒子--正孔非対称性が入ります。]
  実行時にフェルミ面平均 $angle.l tau angle.r$, フェルミ窓での範囲, および $angle.l tau^2 angle.r \/ angle.l tau angle.r^2$ が出力されます。最後の量は定数 $tau$ でちょうど 1 になり, $tau(bold(k))$ が $R_H$ を定数 $tau$ の値からどれだけずらすかを表します。

- `epa_file` （文字列 または `None`）: `tau_mode='epa'` のときに読み込む `epa.x`（`job='egrid'`）の出力ファイルパス。ハミルトニアンのエネルギーは QE 計算と同じ基準（フェルミ準位を 0 にずらしていないもの）である必要があります。
  #detail[散乱率は EPA の原著者による BoltzTraP 実装（`epa_calc`）と同じ形で評価します: $|g|^2$ はエネルギー区間ごとに一定, 終状態の状態密度 $rho(epsilon plus.minus omega)$ は輸送の k メッシュから作った細かい DOS を補間, ギャップ中央より下を価電子帯・上を伝導帯の格子に割り当て, 格子の外の状態は端の区間を使い, $tau lt.eq 1$ ns。区間の端 `edge` は QE 計算の絶対エネルギー（VBM/CBM）なので, Wannier ハミルトニアンのエネルギーも同じ基準である必要があります（フェルミ準位を 0 にずらしたモデルでは区間がずれます）。]
- `epa_pair_file` （文字列 または `None`, 既定 `None`）: 線形 Eliashberg 方程式（option 15, SOC なし）に, この `epa.x`（`job='egrid'`）出力から作ったフォノン媒介の引力を, スピン・電荷揺らぎの対形成相互作用に加えて取り込みます。効くのは s 波（等方的）の成分だけです。
  #detail[相互作用は遅延のある密度--密度型 $W(i nu) = -(1\/eta) sum_nu angle.l g_nu^2 angle.r_"FS" 2 omega_nu \/(nu^2 + omega_nu^2)$ で, 各サイト内の全軌道に等しく働く局所的なもの（Kanamori 表記で $U = U' = W$, $J = 0$）です。$eta$ はこの軌道表示の相互作用のバンド表示での振幅をフェルミ面平均したときの重なりの因子で, フェルミ面平均した等方的な結合が EPA の $lambda = 2 N_F sum_nu angle.l g_nu^2 angle.r \/ omega_nu$ に一致するよう規格化しています。$bold(q)$ に依存しないので効くのは s 波（等方的）の成分だけで, d 波などの符号反転するギャップの固有値は変わりません。フォノンの自己エネルギー（$1 + lambda$ の質量繰り込み）は含まないので, フォノンの寄与はやや過大評価になります。EPA の格子の端は絶対エネルギーなので, ハミルトニアンも同じ基準である必要があります。フェルミ面が EPA の二つの格子にまたがる場合は警告が出ます（金属では一つの格子がフェルミ準位を含むように取ってください）。]
- `epa_pair_scale` （実数, 既定 1）: そのフォノン結合に掛ける倍率。スピン揺らぎとフォノンの競合を調べるときに使います。

- `flex_gamma_imp` （実数, 単位: eV, 既定 0）: `tau_mode='flex'` で $-"Im" Sigma^R$ に加える不純物減衰 $Gamma_"imp"$（マティーセン則）。低温で FLEX の $Gamma prop T^2$ が小さくなりすぎるときの残留抵抗に相当します。

- `flex_npade` （整数, 既定 40）: `tau_mode='flex'` の Padé 接続に使う松原点の最大数（次数 `flex_npade/2`〜`flex_npade` の因果的な結果の中央値を使います）。どの次数も因果的にならない k 点は, 非因果なエネルギー区間を両端から補間して修復し, その数を警告として表示します。FLEX の k メッシュは `self_en_wannier.npz` に記録されます。`freq_grid='ir'` で解いた FLEX の $Sigma$ はサンプリング点の値（`sigma_ir`）から接続し, 点数が `flex_npade` より少なければ全点を使います。`ir_eps=1e-14` 程度を推奨します（option=14 の「周波数格子」の節を参照）。

- `hall_mesh_list` （リスト または `None`, 既定 `None`）: option=5 のホール係数に対する $bold(k)$ メッシュ収束スイープ。例: `[20,30,40]` や `[(20,20,10),(30,30,15)]`。各要素は整数（立方メッシュ）か `(Nx,Ny,Nz)` タプルです。`None` の場合は `(Nx,Ny,Nz)` 単一メッシュのみで, 収束情報は得られません。粗いメッシュでは $sigma^((1))$ が符号ごと誤り得るため, 新しい物質では最低一度はスイープすることを強く推奨します。

- `sw_tdf` （bool）: `True` にすると輸送分布関数（TDF）を先に計算し, そのエネルギー積分から輸送係数を求めます（エネルギー依存緩和時間に対応する場合に有効）。

伝導モード（option=5,6）のスピン縮退因子は `sw_soc` に従います。`sw_soc=False` ではバンドがスピン縮退しているとして $g_s=2$ を, `sw_soc=True` ではハミルトニアンに両スピンが既に含まれているとして $g_s=1$ を使います。

== 軌道・相互作用パラメータ

- `olist` （リスト）: カラープロット（color_option=1）で使用する軌道インデックスのリスト（「カラープロット設定」の節を参照）。

- `U` （実数, 単位: eV）: オンサイトクーロン斥力（ハバード$U$）。FLEX/RPA計算で使用。磁性・超伝導の計算には重要なパラメータです。

- `J` （実数, 単位: eV）: オンサイトフント結合定数。$U' = U - 2J$（金森遮蔽）が自動で計算されます。

- `chiolist`, `site` （整数配列, 既定は未定義）: $chi$ の軌道対の基底（1始まりの軌道番号の組 $(a, b)$ の並び）と, 各組のサイト番号。未定義なら `site_prof` から全軌道の組を作ります。d-p 模型や周期アンダーソン模型のように相互作用が一部の軌道（d, f）にしかない場合は, *相関軌道の組だけ* を並べると, $U$ をその軌道だけに入れた全軌道の計算と同じ RPA・FLEX・線形 Eliashberg の結果を, 組の数 $N_"corr"^2$ の小さな基底で得られます（$Sigma$, $Delta$ は相関軌道にだけ現れます）。SOC ありでは相関軌道の全スピンの組を含めてください（相関軌道 $n$ 個なら $4 n^2$ 組）。出力される $chi$ は相関軌道の部分だけで, 伝導電子を含む全体の感受率ではありません。また FLEX の $Sigma$ の静的部分は相関軌道にだけ入り, $epsilon_d - epsilon_p$ を実質的に動かすので, 静的部分の扱い（二重計上）に注意してください。
- `orb_dep` （bool）: `True`にすると軌道依存の相互作用行列 `Umat`, `Jmat` を使用します。`False`（デフォルト）では `U`, `J` の一定値を全軌道に適用します。`Umat` の対角は軌道内 $U_(i i)$, 非対角は軌道間 $U'_(i j)$ そのもの（$U' = U - 2J$ の金森関係は再適用されません）, `Jmat` は $J_(i j)$ です。いずれも実数行列で, 大きさは $N_"orb" times N_"orb"$, または `site_prof` が同じ軌道数のサイトを複数並べている場合は1サイト分の大きさでも構いません（その場合は各サイトに同じ行列が複製されます。例: 5軌道のセットを2-Fe・10軌道模型へ）。頂点関数が読むのは同一サイト上の軌道対のみなので, サイト間ブロックは参照されません。

- `UJ_material` （文字列, 任意）: `orb_dep=True` のとき, `Umat`, `Jmat` を直接書く代わりに保存済みのパラメータセット（鉄系超伝導体の cRPA $U$/$J$ など）を読み込みます。読み込み元は `UJ_dir`（デフォルト `UJ`）で, 物質ごとのファイルを置いたディレクトリ `{UJ_dir}/{UJ_material}.toml`, または複数物質をテーブル（`[FeSe]` など）で納めた単一の TOML ファイルのどちらでも構いません。各エントリの `U`, `J` はフラット（行優先で $N_"orb"^2$ 個）でも入れ子の行リストでも記述でき, `source`/`doi`（読み込み時に表示）と `orbitals`（ファイルの軌道順の名前。長さを検査）を持たせられます。同梱セット: `FeSe`, `FeTe`, `LiFeAs`, `BaFe2As2`, `LaFeAsO`, `LaFeAsOnakamura`（Fe-$d$ 5軌道の cRPA 値。出典は各ファイルの `source`/`doi` に記載しており, 使用する際はその論文を引用してください）。いずれも1サイト分のセットで, 多サイト模型では `orb_dep` の項で述べたように各サイトへ複製されます。ファイルの軌道順はハミルトニアンの軌道順と一致している必要があります。`None`（デフォルト）: ヘッダに書かれた `Umat`, `Jmat` を使用します。

- `m_diis_num` （整数, 任意）: FLEX (option=14) および非線形 Eliashberg (option=23) で使う DIIS（Pulay 加速混合）の履歴長。$2$ 以上で Pulay 外挿が有効, $1$ で線形混合へフォールバック。未定義の場合は 5 になります。値を大きくするほど収束は速まりますがメモリ消費が増えます（$N_x N_y N_z dot N_w dot N_"orb"^2$ 倍）。

== ギャップ関数の対称性（`gap_sym`）

Eliashberg方程式（option=15,23）を解く際, または SC chi 計算（option=12,13）の初期ギャップ形状を生成する際の対称性を指定します。

#table(
  columns: (auto, 1fr),
  [*gap_sym*], [*対称性*],
  [0], [$s$波（全ての$bold(k)$で正）],
  [1], [$d_(x^2-y^2)$波（$cos k_x - cos k_y$型）],
  [2], [$s^plus.minus$波（ネスティングベクトルを挟んで符号反転）],
  [3], [$d_(x y)$波（$sin k_x sin k_y$型）],
  [4], [$d_(x z)$波（$sin k_x sin k_z$型）],
  [5], [$d_(y z)$波（$sin k_y sin k_z$型）],
  [6], [$d + i d$（カイラル, 複素: index 1 $+ i dot$ index 3）],
  [7], [$d_(x z) + i d_(y z)$（カイラル, 複素: index 4 $+ i dot$ index 5）。$|phi|$ が $k_z=0$ 面全体で消えるため *水平* ライン節を持つ],
  [-1], [$p_x$波],
  [-2], [$p_y$波],
  [-3], [$p + i p$（カイラル, 複素: $p_x + i p_y$）],
)

カイラルな指定（$6, 7, -3$）は複素であり, カーネルもシードも実数である*線形 Eliashberg 方程式（option=15）の初期値には使えません*。実数の相方を別々に計算し（$d+i d$ なら 1 と 3, $d_(x z)+i d_(y z)$ なら 4 と 5, $p+i p$ なら $-1$ と $-2$）, 縮退と直交性を確認したうえで $Delta_1 plus.minus i Delta_2$ を後から組んでください（カイラリティは線形理論では決まらず, $T_c$ 以下の凝縮エネルギーが選びます）。これらは Eilenberger の形状因子（option=24,25,26）用で, そちらではそのまま使えます。

形状因子は模型の格子ベクトルから自動で作られるので, fcc・bcc・六方晶などでも表の対称性を正しく持ちます。鉄系の $s^plus.minus$（`gap_sym=2`）では, 符号を反転させるボンドの長さを `gap_rcut` [Å] で指定できます（例: 結晶軸 $a$）。2-Fe セルのように反転がサイトを入れ替える多サイト模型では, 軌道基底のギャップ（`eil_gap_orbital` / `eil_gap_file`）を使ってください。

#detail[
$cos(2 pi k_x)$ などの簡約座標で書いた形では, これらの調和関数がラベル通りの対称性を持つのは正方晶・直交系の 1 サイト格子だけです。fcc・bcc・六方晶のセルでは $cos(2 pi k_x) - cos(2 pi k_y)$ は周期的ではあっても $B_(1g)$ *ではなく*, 結晶のデカルト $C_(4z)$ で符号が反転しません。そこで形状因子は, 模型の R シェル上の*デカルト格子調和関数*として構成しています。

$ phi_Gamma (bold(k)) = sum_(bold(R) in "shell") K_Gamma (hat(bold(R))) e^(i bold(k) dot bold(R)), quad bold(R)_"cart" = bold(n) dot "avec" $

ここで $K_Gamma$ は各隣接方向について評価したデカルトの調和関数（$hat(n)_x^2 - hat(n)_y^2$, $hat(n)_x hat(n)_y$, $hat(n)_x$ など）です。*軸の分解は一切不要*である点に注意してください。位相 $bold(k) dot bold(R) = 2 pi bold(k)_"frac" dot bold(n)$ は基底に依らないため, デカルト幾何が必要なのはシェルの*係数*だけです。したがってギャップのために `sw_dec_axis` を使う必要はなく, BZ 和や FFT が前提とする $bold(k)$ メッシュの周期性もそのまま保たれます（$bold(R)$ は真の格子ベクトルなので $phi(bold(k)+bold(G)) = phi(bold(k))$ が構成上成立します）。

模型の $bold(R)$ 集合が指定した対称性を担えない場合は, 警告を出して簡約座標の形状因子にフォールバックします。

$s^plus.minus$（`gap_sym=2`）だけは例外で, シェルを対称性からは決められません。$A_(1g)$ はどのシェルも担うため $K_Gamma$ で選別できず, 既定では*簡約座標*で書いたシードパターン $(1, plus.minus 1, 0)$ の長さを代用しています。この代用は基底に依存し, 中心格子では $(1,1,0)$ がデカルトの $[110]$ ではないためカットオフが膨らんで遠方のシェルが選ばれます（$"KFe"_2"As"_2$: 2-Fe セルで $4a$, 1-Fe セルで $2 sqrt(2) a_"Fe"$。どちらも $Gamma$ とゾーンコーナーのポケットの間で符号が変わりません）。どのボンドに符号反転を置くかは対称性ではなく物理の選択なので, シェルの $|bold(R)|$ は `gap_rcut`（`main.py`）で [Å] 単位で明示できます。鉄系では結晶軸 $a$（$"KFe"_2"As"_2$ なら 3.825 Å）が $s^plus.minus$ のシェルで, この 1 つの値で 2-Fe セル（面内第 1 シェル）と 1-Fe セル（第 2 = 対角シェル）の両方が同じ $2 cos(k'_x a_"Fe") cos(k'_y a_"Fe")$ になります。実際に使われたシェルの $|bold(R)|$ はログに出力されます。

注意点が 2 つあります。カイラルな組は 2 つの相方で同じシェルを共有する必要があります（そうしないと, 両者を混ぜる回転に対して $|phi|$ が不変になりません。六方晶の $E_(2g)$ は相方ごとに別シェルを使うと 77% ずれます）。これはコード側で自動的に処理します。もう一つは*多サイト*セルで, ペアリングのボンドは $bold(R)$ ではなく $bold(R) + bold(tau)_j - bold(tau)_i$ なので, サブ格子の位相が依然として欠けます。鉄系の 2-Fe セルのように反転がサイトを入れ替える模型では, 軌道基底のルート（`eil_gap_orbital` / `eil_gap_file` と Nagai–Nakamura のバンド射影）を使ってください。
]

== $bold(k)$空間の設定

- `kz` （実数, 規格化単位）: option=2のフェルミ面プロットにおける$k_z$の値。$0 \sim 0.5$の範囲で指定します（$k_z = 0$は$Gamma$点面, $k_z = 0.5$はゾーン境界面）。

- `kscale` （リスト or 実数）: option=3の3次元フェルミ面表示における各軸の表示スケール。例: `kscale=[1.0, 1.0, 0.5]`とすると$k_z$方向を半分に縮めて表示。既定値の`[1.0, 1.0, 1.0]`のときのみフェルミ面は第一 BZ（Wigner–Seitz 胞）に折り畳まれます。

- `k_sets` （リストのリスト）: symmetry lineの各端点のk座標（規格化単位, $0 \sim 1$）。自動生成ではなく手動で対称線を設定したい場合に定義します。例: `k_sets=[[0,0,0],[0.5,0,0],[0.5,0.5,0]]`

- `xlabel` （文字列のリスト）: `k_sets`に対応するラベル。例: `xlabel=[r'$\Gamma$','X','M']`

- `at_point` （リスト）: option=8で感受率を計算する$bold(q)$点の座標（規格化単位）。

== スイッチ類

- `sw_unit` （bool）: `True`（デフォルト）で国際単位系（SI）の物理定数を使用し, 出力も実際の物理単位で得られます。`False`にすると$hbar = k_upright(B) = e = 1$の無次元系になります。

- `sw_omega` （bool）: option=11 でペアリング感受率を実周波数（`True`）または松原周波数（`False`）で計算するかのスイッチ。

- `sw_self` （bool）: option=4 (スペクトル) / option=15 (線形 Eliashberg) / option=23 (非線形 Eliashberg) で FLEX の自己エネルギーを使うかどうか。`True` にすると繰り込まれたグリーン関数を用います。

- `sw_out_self` （bool）: `True` にすると FLEX の自己エネルギーを `self_en.npz`, `self_en_wannier.npz`（一様格子では `sigma.bin` も）に, option=15, 23 ではギャップを `gap.npy`, `gap_iw.npy`, `gap_wannier.npz` に出力します。

- `sw_in_self` （bool）: `True` にすると FLEX 計算の初期自己エネルギーとして前回の `sigma.bin` を読み込みます（自己無撞着ループの収束を早める）。

- `sigma_in_scale` （float, 既定 1.0）: `sw_in_self=True` で読み込んだ $Sigma$ seed に掛ける係数。$U$ アニーリング用で, $(U_"new"\/U_"old")^2$ が目安（option=14 の説明参照）。

- `sw_from_file` （bool）: `True`（`sw_self=True` と併用）にすると自己エネルギーを `self_en.npz` から読み込み, FLEX 計算を実行せずに直接 Eliashberg 計算を行います。保存時と周波数格子・温度が異なる場合は読み込みを拒否します。

- `sw_check_only` （bool）: option=23（非線形 Eliashberg）専用。`True` で線形化方程式の段（Stoner 因子と $lambda_"eliash"$）で止めます（option=23 の節）。

- `sw_rescale_flex` （bool）: option=14 の FLEX 計算で, 反復中に Stoner 因子が 1 に達した場合に $chi_0$ 全体を一律に $(1-10^(-4))\/S$ 倍して発散を防ぎます。収束解にバイアスは残りません（option=14 の「精度と実装上の近似」）。

- `freq_grid` （`'matsubara'` / `'ir'`, 既定 `'matsubara'`）, `ir_wmax` （float または None）, `ir_eps` （float, 既定 $10^(-8)$）: option=14, 15, 16, 23 の周波数格子。`'ir'` で IR 基底のスパースサンプリングを使います（option=14 の「周波数格子」の節を参照。Padé で実軸へ接続するなら `ir_eps=1e-14` 程度）。

- `sw_chi0_tail` （bool）: option=14,15,23（SOC なし経路）で $chi_0$ を tail 補正付きで評価します（既定 True。IR 基底では無視されます）。打ち切り誤差が $O(1\/N_w)$ から $O(1\/N_w^2)$ に改善し, 同じ精度なら小さい `Nw` で済みます（option=14 の「精度と実装上の近似」）。

- `sw_nosym`, `zeeman_h`, `fflo_q`, `nosym_arnoldi_m`: 対称性を仮定しない線形 Eliashberg（option=15 の節を参照）。

- `imp_cell`, `imp_sites`, `imp_conc`, `imp_seed`, `imp_V`: 不純物スーパーセル（option=21 の節を参照）。

- `sw_dec_axis` （bool, *フェルミ面プロット専用*）: True にすると格子ベクトルを慣用的な直交軸で表し直し, $bold(R)$ も合わせて変換します。これは*描画*のための機能で, 中心格子（とりわけ体心の 122 セル）を見慣れた箱に展開（unfold）し, 10 軌道（2-Fe）のフェルミ面を 5 軌道（1-Fe）のものと直接比較できるようにするためのものです。*物理的な* $bold(k)$ における $E(bold(k))$ は厳密に保存されます。

  バンド図（option=0）とフェルミ面（option=2, 3）でのみ使えます。

  #detail[
  `brav=1,2,6,7` では変換後の $bold(R)$ が半整数成分を持ちます。これは unfold そのものですが, 同時に $H(bold(k))$ が分解後の簡約座標 $bold(k)$ について周期的でなくなることを意味します。したがって BZ 和・FFT 畳み込み・既約 wedge の折り畳み・セル体積（$det("avec")$ は化学式単位あたりの体積ではなくなります）はいずれもこのスイッチ下では無効です。そのため描画系の option=0（バンド図）, option=2, option=3（フェルミ面）でのみ受け付け, それ以外のモードはエラーで停止します。バンド図が許可されるのは経路上で $E(bold(k))$ を評価するだけだからで, DOS（option=1）は BZ 和なので対象外です。option=0 で対称線を自動生成している場合は, 分解後（単純直交）の軸に合わせた経路を生成し直します（中心格子の経路は分解後の軸では別の物理的 $bold(k)$ を指すため）。`k_sets` を自分で指定した場合はそのまま使われるので, 分解後の軸で与える必要があります。六方晶・三方晶（$sqrt(3)\/2$ 成分）は, そもそも整数の直交格子に載せられません。ギャップの対称性についても本スイッチは不要です（デカルトの $bold(R)$ シェルから構成されます。上の `gap_sym` 参照）。
]

=== Eilenberger パラメータ（option=24,25,26）

準古典 Eilenberger ソルバを駆動します。

*共通（3 モード）:*

- `eil_coupling` （float）: 無次元の可分離ペアリング結合 $lambda$（$ |phi|^2 _"FS" = 1$）。大きいほど $T_c$ 上昇。
- `eil_wc` （float, eV）: 固定 Matsubara カットオフエネルギー。ペアリングスケール／$T_c$ を決める。
- `eil_fs_kind` （`None`/`'iso'`/`'ellipse'`/`'tb'`/`'cyl'`/`'sphere'`/`'spheroid'`/`'wannier'`）: フェルミ面。`None`=等方シリンダー（一様侵入長計算は `'ellipse'` にフォールバック）; 面内モデル FS `'iso'`/`'ellipse'`/`'tb'` と 3D モデル FS `'cyl'`（波打ちシリンダー）/`'sphere'`/`'spheroid'` は `eil_fs_params` から生成（3D の場合は `eil_fs_nkz` $> 1$ が必要）; `'wannier'`=読み込んだバンドの実 FS ＋フェルミ速度（対称性／マルチバンドは `gap_sym`, `delta0`, `eil_gap_orbital` から）。
- `eil_fs_params` （tuple）: モデル FS パラメータ — 楕円質量 $(m_x, m_y)$, `tb` ホッピング, `'cyl'` の $(t, t_z)$, `'spheroid'` の $(m_x, m_y, m_z)$。
- `eil_fs_nkz` （整数, 既定 1）: フェルミ面に積む $k_z$ スライス数（`'wannier'` または 3D モデル FS）。$1$=$k_z=0$ の 1 断面（擬 2D）, $> 1$=真の 3D FS。$k_z$ 依存のギャップでは *必須*（$gt.eq 16$）: 水平ライン節は $k_z=0$ 面上にあるため, 1 断面の FS ではギャップが恒等的に 0 になります。
- `eil_fs_traj` （`None`/整数/tuple）: 渦・格子ソルバの軌道削減。コストは FS 点数に比例し, 3D FS では $10^3$–$10^4$ 点（モデルシリンダーは 24 方向）になります。`None`=全点使用, 整数=$beta$（方向）ビン数, `(n_beta, n_v, n_phi)`=方向／$|v_parallel|$ 分位／$phi$ ビンの完全指定。目安として $(48,4,8)$ は 1920 方向の FS を 176 点に削減し $times 17.7$ 高速化, バルクギャップのずれ 0.5%, 渦芯 LDOS ピークのずれ 0.2% です。答えが動かなくなるまで増やしてください。表面ソルバでは使いません（鏡面反射の相方を求めるのに $k_parallel$ とバンド指標が要るため）。
- `eil_pair_gauge` （`'trs'`/`'soc'`/`'diag'`）: バンドギャップ射影におけるペア相方のゲージ。`'trs'`（既定）はスピンレス／実ホッピング用で $phi = u^dagger Delta u$, ゲージ不変。`'soc'` はスピン込み基底で時間反転の相方 $(i sigma_y) u^*$ を使用。`'diag'` は旧実装（$u(-bold(k))$ を独立に対角化）で *ゲージ依存*, 比較用です。要求されるのは *正常状態* の $H(bold(k))$ の時間反転対称性であって $Delta$ のそれではありません — カイラル秩序変数は問題なく（巻きはそのまま $phi$ に渡ります）, 時間反転を破る正常状態（強磁性超伝導, $H$ 内の Zeeman, Haldane 型ホッピング）は不可で, $|angle.l u(-bold(k))| T |u(bold(k)) angle.r| < 1$ のとき警告が出ます。
- `eil_spin_order` （`'block'`/`'interleave'`）: `eil_pair_gauge='soc'` が仮定する基底の並び — `'block'`=[$"orb"_1 dots "orb"_N$ up, $"orb"_1 dots "orb"_N$ down], `'interleave'`=[$"orb"_1$ up, $"orb"_1$ down, $dots$]。
- `eil_imp_gamma` （float, eV）: 非磁性不純物散乱率 $Gamma$（$0$=clean）。
- `eil_imp_c` （float）: T 行列 $cot delta_0$ — 大=Born 極限, $0$=ユニタリ極限。
- `eil_fs_width` （float, eV）: フェルミ面のガウス幅。
- `eil_zeeman` （float, eV）: LDOS 用 Zeeman（Maki）場（表面: $d_[110]$ ZEBS を $plus.minus h$ に分裂; 渦: コア状態をスピン分裂; 自己無撞着周期格子 `eil_vort_lattice_sc`: ギャップ方程式に入る Pauli 対破壊として自己無撞着に扱う）。

*一様（option=24）:*

- `eil_method` （`'normalization'`/`'riccati'`）: $(g, f)$ ルート — `'normalization'` は高速, `'riccati'` は非一様カーネルと整合。
- `eil_find_tc`, `eil_imp_sweep`（＋ `eil_imp_list`）, `eil_pauli`, `eil_spin`, `eil_lambda`, `eil_fs`, `eil_free_energy` （bool）: option=24 の節で述べたサブモード。

*表面（option=25）:*

- `eil_surf_beta` （float, rad）: 表面方位 — $0 = [100]$, $pi\/4 approx 0.785 = [110]$（$d$ 波 ZEBS）。
- `eil_surf_dvector` （bool）: 自己無撞着トリプレット d-vector 表面テクスチャ。
- `eil_dvec_subratio` （float）: d-vector テクスチャの副成分／主成分結合比（$tilde 0.85$ がバルク閾値）。
- `eil_ldos` （bool）: 実周波数の表面／コア LDOS も計算。

*渦・渦格子（option=26）:*

- `eil_field` （float）: $B\/H_(c 2)$ — $0$=孤立渦, $>0$=Doppler 付き円形セル格子。
- `eil_field_list` （list）: *真の周期*格子で走査する $B\/H_(c 2)$ 値（例 `[0.04,0.08,0.16,0.32]`）; `None`=単一磁場。
- `eil_lattice` （`'square'`/`'triangular'`）: 周期格子の形状。
- `eil_kappa` （float）: GL パラメータ $kappa = lambda\/xi$ — 大（$gt.eq 10^3$）=極端第 II 種（遮蔽なし）, 有限=London 遮蔽／Maxwell 反作用。
- `eil_nvortex` （int）: 計算セルあたりの磁束量子数（スーパーセル）。
- `eil_vort_lxi` （float）, `eil_vort_ngrid` （int）: 孤立渦セルの半幅（$xi$ 単位）と 2D グリッドサイズ。
- `eil_vort_field`, `eil_vort_maxwell` （bool）: 自己無撞着な有限 $kappa$ 磁場 $B(rho)$ / ベクトルポテンシャル $bold(A)(bold(r))$。
- `eil_vort_current` （bool）: 環状超流動電流 $j_phi(rho)$。
- `eil_vort_lattice_sc` （bool）: je 流の完全自己無撞着な真の周期格子; `eil_vort_scA=True` で $bold(A)$ を準古典電流から自己無撞着化。
- `eil_vort_scA` （bool）: `eil_vort_lattice_sc` かつ `eil_kappa` 有限のとき, 解析的な London の $bold(A)$ ではなく準古典電流 $bold(j)_s = angle.l bold(v)_F "Im" g angle.r$ から $bold(A)$ を完全に自己無撞着化（je の `A_renew` 方式）。
- `eil_vort_dvector` （bool）: 自己無撞着トリプレット d-vector 渦／格子テクスチャ。
- `eil_vort_chiral` （bool）: 孤立 *カイラル* 渦の自己無撞着解（`eil_chiral_ell=1` で $p+i p$, $2$ で $d+i d$）を多成分・複素振幅のスピン行列ソルバで計算します。芯が *逆* カイラリティを誘起し振幅が本質的に複素になるため, スカラー渦ソルバでは扱えない場合です。両カイラリティは同一の結合定数を共有し（同じ既約表現の縮退した相方）, バルクは等振幅（ネマティック）の重ね合わせではなく単一チャンネルの等方ギャップになります。
- `eil_chiral_ell` （整数）: `eil_vort_chiral` のカイラリティ — $1 = p+i p$, $2 = d+i d$。
- `eil_chiral_m` （整数）: *主成分* の巻き数（$+1$ でカイラリティと同符号, $-1$ で逆）。誘起される逆成分は $m_- = m_+ + 2 ell$ に従います。
- `eil_chiral_dvec` （`'z'`/`'x'`）: カイラルトリプレットの equal-spin 方向（`'y'` は単位行列に比例し, ユニタリな行列 Riccati の初期値としては非対応）。
- `eil_field_dir` （`None` または 3 成分ベクトル）: 磁場／渦糸の方向（孤立渦と有限磁場の格子の両方）。例: $(0,0,1)$=$c$ 軸（既定）, $(1,0,0)$=面内。渦糸は $bold(B)$ に沿って走るため秩序変数はそれに垂直な面内で変化し, 問題は 2D のままです。その面の 2 軸は各々の rms フェルミ速度で規格化され, 楕円的な芯でも正方グリッドに収まります（$xi_1\/xi_2$ が出力されます）。$bold(B) parallel c$ 以外では 3D フェルミ面（`eil_fs_nkz` $> 1$ または 3D の `eil_fs_kind`）が必要です。有限の `eil_field` と異方的な面の組み合わせは未対応で（円形セルの Doppler シフトが規格化前の面を仮定するため）エラーになります。
- `eil_vort_tilt` （float, deg）: $c$ 軸からの磁場傾斜（擬 2D: 軌道 $B_z = B cos theta$, Zeeman $-> h\/cos theta$）。
- `eil_gap_orbital` （`None` / $N_"orb" times N_"orb"$ 行列または callable）: 軌道基底ペアポテンシャルで, その FS バンドへの低エネルギー射影がギャップを決める（Nagai–Nakamura, JPSJ *85*, 074707 (2016), Eq. 43; Wannier FS が必要）, `gap_sym`/`delta0` に優先。
- `eil_gap_file` （`None` / 文字列）: option=15/23（`LIN_ELIASHBERG`/`NONLIN_ELIASHBERG`）を `sw_out_self=True` で実行したとき `output_gap_wannier` が Wannier 実空間「ホッピング」形式で書き出した自己無撞着 RPA/FLEX ギャップの基底名（拡張子なし, 例 `'gap_wannier'`）。指定すると $Delta(bold(R), i omega_n)$ を読み込み, その逆フーリエ変換 $Delta_"orb"(bold(k)) = sum_bold(R) e^(-i 2 pi bold(k) dot bold(R)) Delta(bold(R))$ を FS バンドへ射影して `eil_gap_orbital` として用いる。*以前計算した RPA ギャップ*（例: $"KFe"_2"As"_2$, PRB *84*, 144514）を渦のペアリング形状因子に使うための経路。ギャップを書き出した RPA 計算と Eilenberger 計算は*同一*の Wannier Hamiltonian（同じ軌道基底・$bold(R)$/$bold(a)$ 規約, できれば $mu$/filling）を使う必要がある（バンド固有ベクトルと $Delta_"orb"$ が同じ基底であるため）。`eil_gap_orbital`/`gap_sym` に優先。
- `eil_gap_iw` （int）: `eil_gap_file` の開始 Matsubara インデックス（$0$=最低 $i omega_0$）。Eilenberger の形状因子は静的で, $i omega_0$ が対称性・符号・ノード・異方性を最も鮮鋭に持ち, FS 上で通常図示されるギャップに対応。
- `eil_gap_navg` （int）: `eil_gap_file` で平均する連続 Matsubara スライス数（$1$=単一スライス）。$> 1$ でノイズを平滑化するが異方性をわずかに希釈（$Delta(bold(k), i omega_n)$ は $n$ とともに等方化）。絶対値スケールは無関係（射影後の $phi$ は $ |phi|^2  = 1$ に規格化）。

== NMR パラメータ（option=27）

- `nmr_tc` （実数 または `None`）: $T_c$ [eV]。`None` はフェルミ面上のギャップからの弱結合見積り $max |Delta|_"FS" \/ 1.764$ を使用。
- `nmr_trange` （リスト）: 掃引する換算温度の範囲 $[T\/T_c "最小", T\/T_c "最大"]$。
- `nmr_nt` （整数）: 掃引の温度点数。
- `nmr_qsize` （リスト）: $1\/T_1$ の BZ 和に使う $bold(q)$ メッシュ。`Nx,Ny,Nz` から間引いて作ります（各成分は約数に切り下げ）。コストは $N_q$ に比例し（全体で $tilde N_T N_q N_k$）, `nmr_qsize=[Nx,Ny,Nz]` は $O(N_k^2)$ となりトイモデルでしか現実的ではありません。*分解能* を決めるのは $bold(k)$ 和（$delta gt.tilde v_F d k$ が必要）で, $bold(q)$ 和は滑らかな関数を *積分* するだけ（誤差 $tilde N_q^(-2)$）なので, 通常は各軸 8–16 点で十分です。ただし磁気不安定性の近傍では $chi_s(bold(q))$ が $bold(q)_"AF"$ で鋭くなるため注意が必要です。`nmr_qconv` で確認してください。
- `nmr_qfold` （bool）: $bold(q)$ と $-bold(q)$ を folding してコストを半減（空間反転対称かつ時間反転対称な系で厳密。BdG 基底は既にそれを仮定しています）。
- `nmr_qconv` （bool）: 掃引前に `nmr_qsize` とその半分のメッシュで $bold(q)$ 和を比較し, 相対変化を出力。
- `nmr_w0` （実数 または `None`, 単位: eV）: $chi''(bold(q), omega_0)\/omega_0$ に使う NMR プローブ振動数 $omega_0$。`None` は $0.5 delta$。$omega_0 lt.tilde delta << Delta$ を満たす必要があります。
- `nmr_hf_A`, `nmr_hf_B` （実数）: Mila--Rice の超微細形状因子 $A(bold(q)) = A + 2B(cos q_x + cos q_y)$。$B=0$ で一様な重み, $A = -4B$ で $(pi,pi)$ 応答が消えます。
- `nmr_gap_file` （文字列 または `None`）: option=15/23 が `sw_out_self=True` で出力したギャップ（`output_gap_wannier`）のベース名（拡張子なし）。`gap_sym`/`delta0` の形状因子の代わりに使います。$Delta(bold(R))$ はメッシュに依存しないため, Eliashberg 計算はこの掃引より粗い $bold(k)$ メッシュでも構いません。`gap_sym`/`delta0` は *有効な値* である必要があります（共通の SC chi ブロックが先にそれらからギャップを作り, その後で置き換えるため）が, 値自体は無意味になります。`gap_sym` がペアリングチャンネルを表さなくなるので `nmr_spsym` を明示してください。
- `nmr_gap_extrapolate` （bool）: $Delta(i omega_n) -> Delta(0)$ を外挿し, 最低松原スライスの $O((pi T)^2)$ のバイアスを除去します。`False` の場合は `nmr_gap_iw` のスライスを `nmr_gap_navg` 個平均して使用。
- `nmr_gap_iw`, `nmr_gap_navg` （整数）: `nmr_gap_extrapolate=False` のときに使う松原スライスと平均するスライス数。
- `nmr_gap_max` （実数 または `None`, 単位: eV）: 読み込んだギャップを規格化する全体振幅 $Delta(0)$。*フェルミ面上* の最大値として測ります。線形 Eliashberg のギャップ（固有ベクトルなのでスケールは任意）では必須です。`None` は非線形／自己無撞着ギャップの振幅をそのまま使用。
- `nmr_spsym` （bool または `None`）: `None` は `gap_sym` の符号からシングレット／トリプレットのチャンネルを決めます。`True`/`False` で上書き（`nmr_gap_file` 使用時に必要）。異常項の符号, すなわち Yosida 的に抑制されるチャンネル（シングレット, またはトリプレットで $bold(h) parallel bold(d)$）と $T=0$ で保たれるチャンネル（$bold(h) perp bold(d)$）のどちらかを選びます。

= 典型的な計算フローの例

== フェルミ面とバンドの確認

まずバンド構造とフェルミ面を確認します。

```python
fname, ftype, brav, sw_soc = 'inputs/SomeMaterial', 2, 0, False
option = 0          # バンドプロット
fill = 2.0          # 電子数
abc = [4.0, 4.0, 4.0]
alpha_beta_gamma = [90., 90., 90.]
tempK = 300
```

次にフェルミ面を確認します。

```python
option = 2          # 2次元フェルミ面
color_option = 1    # 軌道成分を色で表示
olist = [0, 1, 2]   # 各軌道をR/G/Bに対応
```

== 超伝導ギャップ関数の計算（RPA 法）

```python
option = 15         # 線形 Eliashberg 方程式
Nx, Ny, Nz, Nw = 32, 32, 1, 512
tempK = 50          # やや低温で計算
fill = 2.0
U, J = 1.0, 0.1
gap_sym = 1         # d_{x^2-y^2} 波で初期化
sw_self = False     # RPA で計算（FLEX 自己エネルギーなし）
sw_out_self = True  # ギャップ関数を出力
```

計算が終わると標準出力に固有値 $lambda$ が表示され, `gap_{ij}.dat` および `gap.npy` にギャップ関数が出力されます。

== FLEX + 線形 Eliashberg

自己エネルギーを繰り込んだより精密な計算の場合は, まず option=14 で FLEX 計算を行い, 次に option=15 を `sw_self=True`, `sw_from_file=True` で実行します。

```python
# ステップ 1: FLEX 自己エネルギー計算
option = 14
sw_out_self = True
m_diis_num = 5

# ステップ 2: 線形 Eliashberg 方程式（FLEX 自己エネルギー込み）
option = 15
sw_self = True
sw_from_file = True
```

== 非線形 Eliashberg 方程式（自己無撞着 SC ループ）

線形 Eliashberg で $lambda approx 1$ となる温度より低温に降ろし, $Delta$ を有限振幅まで成長させる場合は option=23 を使用します。初期 Δ（対称性形状と BCS スケールの振幅）は内部で自動生成されるため, 事前に option=15 を実行して `gap.npy` を用意する必要はありません。

```python
# まず sw_check_only=True で Tc 近傍を確認（非線形ループはスキップ）
option = 23
gap_sym = 1
sw_self = True            # FLEX 自己エネルギー込み
sw_from_file = True
sw_check_only = True
```

$lambda_"eliash" > 1$ となる温度まで降温したら, `sw_check_only=False` にして非線形ループを実行します。

```python
option = 23
sw_self = True            # FLEX 自己エネルギー込み
sw_from_file = True
sw_check_only = False
m_diis_num = 5            # DIIS Pulay 加速 + 振幅方向 Newton 加速
```

低温（$T tilde.equiv T_c \/ 5$）では DIIS 履歴を 5–10 程度にすると収束が速くなります。

== IR 基底での計算

低温や, バンドが高エネルギーまで延びる模型（MLO 基底など）では, 一様な松原周波数より IR 基底の方が少ない周波数点で正確です。上の FLEX + 線形 Eliashberg の流れで `freq_grid = 'ir'` を設定するだけです（設定は option=14 の「周波数格子」の節）。

== 超伝導状態の動的スピン感受率

ギャップが開いた状態で $chi_s^"SC"(bold(q), omega)$ を計算し, スピン励起の SC ギャップ依存性を見る場合：

```python
option = 12               # symmetry line 沿いの SC chi スペクトル
delta0 = 1.e-2            # 初期ギャップ振幅 (eV) ≒ 10 meV
gap_sym = 1               # d_{x^2-y^2} 波
U, J = 0.8, 0.1
tempK = 50
```

`option = 13` にすると `at_point` で指定した単一 $bold(q)$ 点での $chi_s^"SC"(omega)$ が得られます。

= AiiDA インターフェース

`examples/aiida/pyrpa_aiida.py` は, pyrpa の計算を AiiDA の `calcfunction` / `workfunction` として実行するサンプルです。各段階の入力・出力とそのつながり（来歴, provenance）が AiiDA のデータベースに記録されます。計算そのものは同じプロセスの中で `libs.plibs` / `libs.flibs` を直接呼んで行うため, computer・code・デーモンの設定は要りません。`main.py` をクラスタに投入する `CalcJob` は未実装です。

== 導入

```bash
pip install aiida-core            # 2.9.4 で動作確認
```

普段の計算環境を変えたくない場合は, MKL 版の numpy/scipy をそのまま使える仮想環境に入れてください。

```bash
python -m venv ~/venvs/pyrpa-aiida --system-site-packages
source ~/venvs/pyrpa-aiida/bin/activate
pip install aiida-core
```

テストと `run_example.py` は既定で一時プロファイル（`SqliteTempBackend`）を使うので, PostgreSQL・RabbitMQ・`verdi presto` は不要です（終了時に何も残りません）。結果を残す場合は `verdi presto` などでプロファイルを作り, `run_example.py --profile NAME` で指定します。

== 提供するステップ

#table(
  columns: (32%, 1fr, 26%),
  align: (left, left, left),
  table.header([*ステップ*], [*内容*], [*主な出力*]),
  [`load_hoppings`], [模型の読み込み（`ftype` は `main.py` と同じ。ディレクトリ形式は `FolderData`, ftype=2 は `_hr.dat` ファイルを `SinglefileData` で渡す）], [`rvec`, `ham_r` (, `S_r`, `avec`, `s_rcut`)],
  [`chemical_potential`], [裸のバンドの化学ポテンシャル（全 $bold(k)$ メッシュ）], [`mu`],
  [`flex_self_energy`], [FLEX 自己エネルギー（`mkself` / SOC では `mkself_soc`）], [`sigma`, `mu_self`, `flex_info`],
  [`rpa_susceptibility`], [$nu = 0$ の $chi_s$, $chi_c$（`sigma` を渡すと FLEX）], [`chi`, `stoner`],
  [`linearized_eliashberg`], [線形 Eliashberg 方程式（RPA / FLEX）], [`lambda_eliash`, `gap`, `gap_wannier`],
  [`flex_eliashberg`], [FLEX → $chi$ → $lambda$ を 1 ステップで実行し, $Sigma$ を保存しない], [上 3 つの結果],
  [`eilenberger_bulk`], [均一系の Eilenberger 方程式（`gap_wannier` を形状因子にできる）], [`result`, `f_avg`],
  [`eilenberger_surface`, \ `eilenberger_vortex`], [表面・渦糸（モデルフェルミ面）], [`profile`, `files`],
)

これらをつないだ workflow として `flex_workflow`（模型 → $mu$ → FLEX）, `eliashberg_workflow`（模型 → $mu$ → [FLEX] → $chi$, $lambda$）, `eliashberg_eilenberger_workflow`（その後 `eilenberger_bulk`）があります。RPA/FLEX のステップには全て `calcfunction` の出力として `runtime`（OpenMP のスレッド数と待ち方, `libfmod` のビルド, numpy の版, git のコミット）が付きます。

配列は pyrpa の配置のまま保存します（$Sigma$, $Delta$ は既約 $bold(k)$ メッシュ上の `[Norb, Norb, Nw, Nk]`）。MLO 基底（ftype=3,4）の扱いは `main.py` と同じです。

格子調和関数を使う場合, `load_hoppings` または workflow に `lattice=orm.Dict({"avec": [[...], [...], [...]], "s_rcut": 4.0})` を渡します。`avec` は原始並進ベクトルを行に並べたもの（angstrom）, `s_rcut` は省略可能な $s_+$-$s_-$ の殻の距離（angstrom）です。格子はホッピングと共に記録され, Eliashberg の初期ギャップと `eilenberger_bulk` の解析的形状因子に使われます。未指定の場合は警告の上で分数座標の式を使います。この式は非直交格子や複数サイトの格子では対称性が正しくない場合があります。他の計算が登録した格子は引き継ぎません。`gap_wannier` を指定した場合はそれが Eilenberger の形状因子になります。

== オプション（`options`）

RPA/FLEX のステップと workflow は省略可能な `options`（`orm.Dict`）を取ります。`main.py` の冒頭で設定する項目に対応し, 既定値は `OPTION_DEFAULTS` にあります。知らないキーはエラーになります。

各ステップはプロセス内ですでに有効な格子に依存せず指定された周波数格子を使い,
終了時や失敗時には以前の格子を復元します。格子はプロセス全体で共有されるため,
同じプロセスの複数スレッドで同時に切り替えてはいけません。

#table(
  columns: (24%, 14%, 1fr),
  align: (left, left, left),
  table.header([*キー*], [*既定値*], [*内容（`main.py` の対応項目）*]),
  [`soc`], [`False`], [スピン軌道相互作用（`sw_soc`）。模型はスピン込みで, U, J・`Umat`・`site_prof` はスピン軌道について解釈します],
  [`slist`, `invs`], [`None`], [スピンのラベル（$plus.minus 1$）と 1 始まりの逆スピン相手。`None` で `main.py` と同じ「全 up の後に全 down」],
  [`site_prof`], [`None`], [サイトごとの軌道数（例: `[5, 5]`）。`None` は全軌道を 1 サイト],
  [`chiolist`, `site`], [`None`], [$chi$ の軌道対（1 始まり）とそのサイト。例えば d 軌道だけに U を置くとき d--d の対だけにすると速くなります],
  [`Umat`, `Jmat`], [`None`], [軌道依存の U, J（`orb_dep`）。1 サイト分を与えると全サイトに複製します（`plibs.expand_UJ_sites`）],
  [`freq_grid`], [`'matsubara'`], [`'ir'` で IR 基底（`freq_grid`）。周波数点の数は基底が決めます],
  [`ir_wmax`, `ir_eps`], [`None`, `1e-8`], [IR 基底の打ち切り（`None` で $2 max|epsilon_bold(k) - mu|$）と精度],
  [`sw_tail`], [`False`], [テイル補正 $chi_0$（`sw_chi0_tail`）。SOC なし・一様格子のみ],
  [`eps`, `pp`, `m_diis`, \ `scf_loop`, `sub_sigma`, \ `sw_rescale`], [`1e-4`, `0.3`, `5`, \ `300`, `1`, \ `False`], [FLEX ループの設定（`flibs.mkself`）],
  [`require_converged`], [`True`], [FLEX が収束しない, または Stoner 因子 $>= 1$ で終わった場合にステップを失敗させる],
  [`store_gap`], [`True`], [`False` で全周波数の $Delta$ を保存せず `gap_wannier` のみ],
  [`store_sigma`], [`True`], [`False` で FLEX 付き `eliashberg_workflow` を `flex_eliashberg` の 1 ステップで実行し $Sigma$ を保存しない],
)

- *FLEX の収束（`flex_info`）.* Fortran の FLEX ループは収束したかどうかを返さないため, その出力ログ（`iter=`, `esterr=`）から反復回数・最後の esterr・収束の可否を記録します。Stoner 因子 $>= 1$ で終わった解は磁気的に不安定な非物理的な不動点で, 不安定性の近くでは反復の初期の振れでそこに落ちることがあるため, 収束とはみなしません。
- *IR 基底の引き継ぎ.* FLEX の $Sigma$ には解いた IR 基底（`ir_wmax`, `ir_eps`）が記録され, それを使う後段のステップは同じ基底を使います（$mu_"self"$ に合わせて打ち切りを選び直すことはしません）。
- *SOC.* $chi$ はスピンのチャンネル `chiszz`, `chispm`, `chisxx`, `chisyy` と `chic` を持ちます（`chis_q` は zz の先頭固有値）。Stoner 因子は出力しません（SOC の $chi_0$ には無く, `main.py` も出しません）。Eliashberg は対角パリティの模型のみで, `sw_tail` とは併用できません。SOC のギャップを `eilenberger_bulk` の形状因子にすることは未対応です。
- *MLO と周波数の打ち切り.* 一様格子では `sw_tail=True` か十分大きな `Nw`, あるいは `freq_grid='ir'` を使ってください（トラブルシューティングの MLO の項）。

== 実行例

`examples/aiida/run_example.py` で主なステップを試せます。

```bash
python examples/aiida/run_example.py --calc flex                      # square.hop: 模型 -> mu -> FLEX
python examples/aiida/run_example.py --calc eliashberg --flex         # FLEX + 線形 Eliashberg
python examples/aiida/run_example.py --calc eliashberg --flex --soc   # SOC（square_2orb_soc.hop, fill = 2）
python examples/aiida/run_example.py --calc eliashberg --flex --mlo-demo   # 2 軌道 MLO 模型, IR 基底
python examples/aiida/run_example.py --calc eliashberg --hop model_hr.dat --ftype 2 \
    --lattice '{"avec": [[4,0,0],[0,4,0],[0,0,6]], "s_rcut": 4.0}'
python examples/aiida/run_example.py --calc eliashberg --hop inputs/Cu --ftype 4 --fill 5.5 \
    --options '{"freq_grid": "ir", "require_converged": true}'   # ecalj job_mlo のディレクトリ
```

`--options` には上の表のキーを JSON で渡します。Python から直接使う場合は次のようになります。

```python
import aiida
from aiida import orm
from aiida.storage.sqlite_temp import SqliteTempBackend
aiida.load_profile(SqliteTempBackend.create_profile('pyrpa'), allow_switch=True)
import pyrpa_aiida as pa          # examples/aiida を sys.path に入れておく

res = pa.eliashberg_workflow(
    hop_file=orm.SinglefileData('inputs/square.hop'), ftype=orm.Int(1), fill=orm.Float(0.45),
    temp=orm.Float(0.02), U=orm.Float(0.3), J=orm.Float(0.0), mesh=orm.List([16, 16, 1]),
    Nw=orm.Int(128), gap_sym=orm.Int(1), use_flex=orm.Bool(True),
    options=orm.Dict({'freq_grid': 'ir', 'store_sigma': False}))
print(res['lambda_eliash'].value, res['flex_info'].get_dict())
```

= トラブルシューティング

- *化学ポテンシャルが収束しない*: `Nx, Ny, Nz` を増やすか, `tempK` を少し上げてください。

- *スペクトルにノイズが多い*: `delta` を少し大きくする, または `Nw` を増やしてください。

- *バンドがおかしい*: `brav` の設定が合っているか確認してください。特に QE との連携では `brav=1`（FCC）や `brav=2`（BCC）の convention に注意が必要です。

- *線形 Eliashberg が収束しない*: `U` を小さくする, `tempK` を高くする, メッシュ `Nx,Ny,Nz` を増やすなどを試してください。固有値 $lambda$ が 1 を超えている温度は超伝導相内に対応します。

- *非線形 Eliashberg が発散する*: `m_diis_num` を 5–10 に増やす, `tempK` を上げる, `Nw` を増やすなどを試してください。`sw_check_only=True` で先に Stoner 因子と $lambda_"eliash"$ を確認するのも有効です（option=23 の節）。

- *FLEX で max$|Sigma|$ が発散*: `sw_rescale_flex=True` を設定するか, `U` を下げてください。

- *`[FLEX] WARNING: Stoner factor >= 1 at the final iteration` が出る*: option=14 の注意枠を参照してください。不安定性の近くでは反復の初期の振れで非物理的な不動点に落ち, 混合率や周波数格子の違いで別の値に「収束」することもあります。

- *gfortran でビルドしたライブラリが非常に遅い*: OpenMP ランタイムの競合です（「Fortran カーネルのビルド」の節）。`_loader.py` が `OMP_WAIT_POLICY=PASSIVE` を自動で設定しますが, 環境変数で別の値を指定している場合はそれが優先されます。

- *計算ノードで `Illegal instruction`（SIGILL）で落ちる*: ログインノードの CPU 向けにビルドしたライブラリです。`make auto` でビルドし直してください（「クラスタで使う場合」）。

- *MLO 模型で Stoner 因子や $lambda$ が異常に大きい / 周波数の端でギャップや $Sigma$ が乱れる*: MLO のバンドは高エネルギーまで延びるため（Cu で約 $+30$ eV）, 一様格子で窓 $(2 N_w - 1) pi T$ が足りないと $chi_0$ の打ち切り誤差が大きくなります。`sw_chi0_tail=True` のまま `Nw` を増やすか, `freq_grid='ir'` を使ってください。U を d 軌道だけに置く dp 模型では, `orb_dep` の `Umat`/`Jmat` で s, p の行と列を 0 にし, `chiolist` を d--d の組に限ると速くなります。

- *伝導率の単位が合わない*: `sw_unit=True` になっているか確認してください。`False` にすると無次元系になります。

- *ホール係数が $bold(k)$ メッシュで符号が変わる*: これは $sigma^((1))$ の典型的な破綻の仕方であり, バグではありません。`hall_mesh_list` でスイープし, 2つの診断値を見てください。$N_"eff"$ が $10^3$ を下回っていればフェルミ窓を担う状態数が足りず, 相殺率が $10^(-2)$ を下回っていれば補償近傍で $sigma_(x x)$ よりはるかに細かいメッシュが必要です。`tempK` を上げてフェルミ窓を広げるのも有効です。

- *ホールキャリア密度が実験と合わない*: 同じ出力に併記される Luttinger 体積と突き合わせてください。$n_e \/ n_h approx 1$ なら補償金属であり, $n_H$ はそもそもキャリア密度ではありません（移動度の比較や二流体フィットを行ってください）。また `sw_soc` が正しく設定されているか（スピン縮退因子を決めます）, 実験磁場での $omega_c tau$ が $lt.tilde 0.1$ に収まっているかも確認してください。

- *$R_H$ に温度依存性が出ない*: `tau_mode='const'` では当然の挙動です（`tau_const` の項）。`tau_mode='epa'` などで $bold(k)$ 依存の $tau$ を入れてください。

= 付録: テストコード

`tests/` に数値カーネルと物理ベンチマークの回帰テストがあります。Fortran ライブラリをビルドしてから, リポジトリのルートで次を実行します（`pytest` が無くても動きます）。

```bash
python tests/run_all.py            # 全ファイル
python tests/run_all.py rpa        # ファイル名に "rpa" を含むものだけ
python tests/run_all.py --runslow  # slow マーカー付きのテストも含める
```

`pytest` があれば `pytest tests` でも実行できます。各ファイルの内容は `tests/README.md` と `Doc/internals_JP.typ` を参照してください。
