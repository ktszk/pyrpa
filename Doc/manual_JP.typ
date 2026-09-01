//文章フォーマット設定
#set text(font: "Times New Roman", size: 11pt)
#show regex("[\u3000-\u9fa5]"): set text(font: "Yu Gothic", lang: "ja")
#set page(paper: "a4", margin: (x: 2cm, y: 2cm))
#set math.equation(numbering: numbering.with("(1)"), supplement: "式")
#set heading(numbering: "1.")
//外部パッケージ
#import "@preview/physica:0.9.3": *

#align(center)[
  #text(size: 18pt, weight: "bold")[pyrpa マニュアル]
  #v(0.5em)
  #text(size: 11pt)[2025年]
]

#v(1em)

= pyrpaについて

pyrpaは, モデルハミルトニアン（タイトバインディング模型）を入力として, 様々な物性量を計算・可視化するPythonコードです。主にWannier関数に基づいた第一原理計算の結果をポスト処理する用途を想定しています。

現在サポートしている主な計算機能は以下の通りです。

- バンド構造・フェルミ面（2D / 3D）の計算と可視化
- 状態密度（DOS）の計算
- スペクトル関数 $A(bold(k), omega)$ の計算
- ボルツマン理論・線形応答理論による各種伝導係数（電気伝導率, 熱伝導率, ゼーベック係数など）の計算
- 弱磁場ホール係数 $R_H$, ホールキャリア密度, ホール移動度の計算（バンド分解による二流体分析付き）
- RPA（乱雑位相近似）に基づくスピン感受率 $chi_s$ およびペアリング感受率 $phi$ の計算
- 超伝導状態における動的スピン感受率 $chi_s^"SC"$ の計算
- FLEX（揺らぎ交換近似）による自己エネルギーの計算
- 線形 Eliashberg 方程式による超伝導固有値・ギャップ関数の計算
- 非線形 Eliashberg 方程式による自己無撞着 SC ギャップ関数の計算
- バンドフィリング・サイクロトロン質量・dHvA 振動数の計算
- 準古典 Eilenberger/Riccati ソルバ: $T_c$, 磁場侵入長, 表面 Andreev 束縛状態, 渦・渦格子
- 超伝導状態の NMR: BdG スピン帯磁率から Knight shift $K(T)$ と $1\/(T_1 T)$
- 不純物・CPA を含むスペクトル関数の計算

pyrpaはすべての設定を`main.py`の上部にある変数を書き換えることで行います。ライブラリのimport行より上にある変数のみを変更することを想定しています。

= 本改訂での変更点

以下の変更は過去のバージョンで得られた結果を変えます。特に最初の 2 点は過去の計算結果そのものを無効にし, ギャップ調和関数の変更は中心格子・六方晶セルでの結果を全て変えます。残りは出力内容, あるいは受け付ける設定が変わります。本改訂より前の出力がある場合は, 下記に照らして確認してください。

- *時間反転対称性下の既約 $bold(k)$ メッシュ.* `generate_irr_kpoint_inv` は $bold(k)$ 点の座標を許容誤差 $1\/max(N_x,N_y,N_z)$ — ちょうどメッシュ間隔 1 つ分 — で照合しており, $1\/N$ が2進小数で表せない場合に *隣の点* に誤マッチしていました。$N$ が 2 のべき乗（$4,8,16,32,64$）のメッシュは影響を受けませんが, それ以外は影響を受けます（例: $12^3$ では 1728 点中 919 点が誤対応）。現在は整数メッシュ添字から構成しており, 全てのメッシュで厳密です。2 のべき乗以外のメッシュで得た FLEX / Eliashberg / SC-$chi$ の結果は信頼できないものとして再計算してください。同ルーチンは $O(N_"kall"^2)$ でしたが $O(N_"kall")$ になりました（$64^3$ で 9.2 s $arrow$ 0.009 s）。

- *$N_x$ が奇数かつ $N_y$ が偶数の既約ウェッジ.* このパリティの組み合わせは `gen_irr_k` の未実装分岐であり, ウェッジが埋められないまま `klist` の範囲外に書き込むため $N_z >= 3$ でクラッシュしていました。現在は実装済みで, $8 times 8 times 6$ までの 384 通りのメッシュ／パリティ全てについて `TRS_irr.typ` の記述と一致することを確認しています。

- *伝導計算のスピン縮退.* option=5, 6 はスピン因子を `sw_soc` から取るようになりました。従来は無条件に $g_s=2$ としていたため, 両スピンを含むハミルトニアンでは結果が 2 倍ずれていました。

- *ホール係数.* $R_H$ は $sigma_(x y)\/(sigma_(x x) sigma_(y y))$ ではなく完全な伝導率テンソルから構成するようになり, ホールキャリア密度はキャリアの型を添えた正の値として出力されます。option=5 を参照してください。

- *ギャップ形状因子をデカルト格子調和関数で構成.* `gap_sym` は従来 *簡約座標* の関数（$cos(2 pi k_x) - cos(2 pi k_y)$ など）として評価しており, これがラベル通りの対称性を持つのは正方晶・直交系の 1 サイト格子だけでした。fcc・bcc・六方晶のセルでは, その関数は周期的ではあっても $B_(1g)$ ではなく, 結晶のデカルト $C_(4z)$ で符号が反転しません（実測: 破れは $O(1)$, 六方晶の $E_(2g)$ 二重項ノルムは 100% ずれます）。現在は模型の $bold(R)$ シェル上の $sum_bold(R) K_Gamma (hat(bold(R))) e^(i bold(k) dot bold(R))$ として構成し, $K_Gamma$ をデカルト空間で評価します。これによりラベル通りの既約表現と $bold(k)$ 空間の周期性が任意のブラベー格子で成立します（$10^(-15)$ で確認）。直交系の 1 サイト格子では全ての `gap_sym`・任意の軸比について従来の式をビット単位で再現するため, そこでの結果は変わりません。中心格子や六方晶のセルで自明でない `gap_sym` を使った計算 — SC-$chi$（option=12,13）, NMR（option=27）, Eilenberger の形状因子（option=24,25,26）, Eliashberg の初期値（option=15,23） — は値が変わり, 新しい方が正しい値です。構成方法と多サイトセルに残る制約は `gap_sym` の節を参照してください。

- *`sw_dec_axis` を描画モード限定に.* 軸の分解は中心格子を慣用的な直交軸に unfold するもので, 10 軌道（2-Fe）のフェルミ面やバンド構造を 5 軌道（1-Fe）のものと比較できるようにするための機能です。一方で $bold(R)$ に半整数成分が残るため, 分解後の簡約座標 $bold(k)$ について $H(bold(k))$ が周期的でなくなり, BZ 和・FFT 畳み込み・既約ウェッジの折り畳み・セル体積としての $det("avec")$ がいずれも無効になります。現在は option=0, 2, 3 でのみ受け付け, それ以外のモードは静かに誤った値を返す代わりにエラーで停止します。従来この機能の動機だったギャップの対称性については, 上項の通り不要になりました。

= 動作環境

以下のPythonパッケージが必要です。

- `numpy`
- `scipy`
- `matplotlib`
- `skimage`（フェルミ面プロット option=2,3, サイクロトロン質量 option=18, dHvA 軌道 option=19 で使用）

Fortran カーネル `libs/libfmod.so` をあらかじめコンパイルしておく必要があります。`libs/flibs`, `libs/plibs` はこれを呼び出す純 Python パッケージです。線形 / 非線形 Eliashberg・FLEX・SC chi 計算は OpenMP 並列の Fortran ルーチン経由で実行されます。

== Fortran カーネルのビルド

`libs` ディレクトリで以下のいずれかを実行します。

#table(
  columns: (auto, 1fr),
  [`make auto`], [*クラスタでの推奨.* ツールチェインはインストール／ロード状況から, CPU ターゲットは計算ノード自身から自動判定します。見つかった CPU の種類ごとに `libfmod_<target>.so` を 1 つずつビルドし, `libs/flibs/_loader.py` が実行ノードで使える最良のものを import 時に選びます。],
  [`make`], [単一ターゲットのライブラリを 1 つ。`ARCH`（既定 `znver3` = 可搬な AVX2 ベースライン）, `FC`, `SL` は自動判定またはコマンドラインの指定に従います。],
  [`make dispatch`], [`znver3` + `znver4` の固定ペア。すなわち検出を省いた `make auto` です。],
)

明示指定が常に優先されます: `make FC=ifx SL=MKL ARCH=znver4`, `make auto ARCH_LIST="znver3 znver4"` など。

*何を, どこから判定するか.* `FC` は `PATH` 上の `ifx`, `ifort`, `flang`, `amdflang`, `gfortran` の最初に見つかったもの, `SL` は実際にロードされているモジュール（`MKLROOT` $arrow$ MKL, `AOCL_ROOT` $arrow$ AOCL, どちらも無ければ OpenBLAS）に従います。root が設定されていないライブラリを選ぶと, BLAS シンボルが未解決のまま import 時に初めて落ちる `libfmod.so` ができてしまうためです。一方 *CPU ターゲット* はビルドホストから取ってはいけません。Zen 4 のログインノードで `-march=native` を使うと, Zen 3 の計算ノードでは SIGILL で落ちる AVX-512 コードができます。そこで `libs/detect_arch.sh` は, スケジューラの台帳（Slurm の `sinfo`, PBS の `pbsnodes`, いずれも無ければ自機のみ）にある各ノードに対して `gcc -march=native` がその CPU を何と呼ぶかを問い合わせ, 重複を除いたターゲットを出力します。

```
$ sh libs/detect_arch.sh --table
NODE         CPU                                      -march
salmon       AMD Ryzen 7 5700G with Radeon Graphics   znver3
unagi        AMD Ryzen 9 7900 12-Core Processor       znver4
yakisoba     AMD Ryzen 9 7900X 12-Core Processor      znver4
```

結果は `libs/.arch_cache` にキャッシュされます。`sh libs/detect_arch.sh --refresh` で再取得, `--local` で自機のみに限定できます。応答しないノードは黙って無視されるのではなく, 可搬な既定値 `znver3` を寄与したうえで警告が出ます。`slurm.conf` に `Features=znver3` / `znver4` を設定しておけば, 検出はジョブ投入なしの `sinfo` 参照だけで済みます（クラスタを管理できる場合はこちらを推奨）。

実行時に手で選ぶ必要はありません。`_loader.py` が CPU のフラグを読み, そのマシンで実際に実行できる最上位 ISA の `libfmod_*.so` を読み込みます（見つからなければ通常の `libfmod.so`）。AMD 上で MKL を使う場合は `MKL_ENABLE_INSTRUCTIONS=AVX512`（Zen 4）または `AVX2`（Zen 3）も export し, `OMP_NUM_THREADS` は決め打ちではなくスケジューラの値（`$SLURM_CPUS_PER_TASK`）から設定してください。

*CMake.* `libs/CMakeLists.txt` でも同じ判定・同じ成果物が得られます（CMake を好む場合や IDE 連携が必要な場合）。

```bash
cmake -S libs -B build && cmake --build build -j     # = make auto
cmake -S libs -B build -DARCH=znver4                 # = make ARCH=znver4（libfmod.so 単体）
cmake -S libs -B build -DFC_COMPILER=gfortran -DSL=MKL -DARCH_LIST="znver3;znver4"
```

各変数は空なら自動判定です（`FC_COMPILER` は `PATH` から, `SL` は `MKLROOT`/`AOCL_ROOT` から, ターゲット一覧は `detect_arch.sh` から）。ライブラリはビルドツリーではなく `libs/` に出力します（`_loader.py` がそこを見るため）。C++ コンパイラが必要なのは pocketfft 経路のみで, MKL 経路では CXX 言語自体を有効にしないため, `ifx` はあるが `icpx` が無い oneAPI 環境でも configure できます。なお参照実装は makefile 側です。`flang` と AOCC モジュールを併用する場合, makefile は子プロセスの `CPATH` と `LIBRARY_PATH` の掃除も行いますが, CMake 側は行いません。

= 入力ファイルについて

== ファイル形式の選択（`ftype`）

ハミルトニアンの入力ファイル形式を`ftype`変数で指定します。`fname`にファイル名またはディレクトリ名（拡張子なし）を文字列で代入してください。

#table(
  columns: (auto, 1fr),
  [*ftype*], [*ファイル形式*],
  [0], [`{fname}/` ディレクトリ内に `ham_r.txt`, `irvec.txt`, `ndegen.txt` が存在する形式],
  [1], [`{fname}.input` ファイル（独自形式）],
  [2], [`{fname}_hr.dat` ファイル（Wannier90のデフォルトホッピングファイル）。*第一原理計算との連携では通常これを使用*],
  [3], [MLO（最局在ワニエ）基底の非直交基底ホッピング形式],
  [4], [`{fname}/` ディレクトリ内の `HamRsMLO`, `HamiltonianPMTInfo`（ecalj `job_mlo` のバイナリ出力）。ftype=3 と同じく非直交 MLO 基底で, 重なり積分 $S(R)$ も同時に読み込みます],
  [その他], [`Hopping.dat` ファイル（ecaljのホッピングファイル）],
)

*Wannier90を使用する場合の例：*

```python
fname = 'inputs/Sr2RuO4'   # _hr.dat を除いたパスを指定
ftype = 2
```

== スピン軌道相互作用（`sw_soc`）

`sw_soc = True` にすることでスピン軌道相互作用（SOC）を考慮した計算を行います。SOC を含む系では, ハミルトニアンがスピン自由度を含む（軌道数が 2 倍になる）形式である必要があります。FLEX (option=14) および線形 Eliashberg (option=15) は SOC 対応版が用意されています（内部で `calc_flex_soc` / `calc_lin_eliash_soc` を呼びます）。非線形 Eliashberg (option=23) と SC chi (option=12,13) は現状 SOC 非対応です。

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

計算内容は `option` 変数の整数値で指定します。`CalcMode` IntEnum の値を直接指定することもできます（例: `option = CalcMode.LIN_ELIASHBERG`）。各モードの詳細を以下に示します。

#block(stroke: 1pt + red, inset: 8pt, radius: 4pt)[
  *注意: 未検証モードについて*

  以下のモードは `main.py` の `CalcMode` 定義で *(not implemented)* と注記されており, 実装の正確性が十分に検証されていません。結果の物理的妥当性・数値精度は保証されないため, 使用する場合は *自己責任* でご利用ください。出力を論文等で引用する前には独立な手段による検証を強く推奨します。

  - option=20 (`ELECTRON_MASS`): symmetry line 上の電子質量
  - option=21 (`SPECTRUM_IMPURITY`): 不純物入りスペクトル関数
  - option=23 (`NONLIN_ELIASHBERG`): 非線形 Eliashberg 方程式
]

== option=0: バンドプロット (`BAND`)

symmetry line に沿ったバンド分散 $E_n(bold(k))$ を計算し描画します。対称線は `k_sets` と `xlabel` で指定するか, 指定がなければ `brav` の設定から自動で選びます。

== option=1: 状態密度 (`DOS`)

$ "DOS"(omega) = - 1/pi sum_(bold(k), n) "Im" G^0_n (bold(k), omega + i delta) $

全軌道の合計と軌道ごとの部分 DOS をプロットします。

== option=2: 2 次元フェルミ面 (`FERMI_2D`)

指定した $k_z$ 断面（デフォルト: $k_z = 0$, 変数 `kz` で変更可能）における $k_x$-$k_y$ 面のフェルミ面を描画します。電子数 `fill` またはフェルミエネルギーから化学ポテンシャルを自己無撞着に求め, $E_n(bold(k)) = mu$ となる等エネルギー線を skimage の find\_contours で抽出します。`RotMat` で回転行列を指定することでフェルミ面を回転できます。

== option=3: 3 次元フェルミ面 (`FERMI_3D`)

Marching cubes アルゴリズムを用いて 3 次元フェルミ面をポリゴンで描画します。`kscale` で各軸方向の表示スケールを変えることができます。`color_option=ColorMode.GAP`（3）にすると, Eilenberger 計算を駆動している*まさにそのペアリング形状因子* $"Re"[phi(bold(k))]$（`gap_sym`/`delta0`, または設定されていれば `eil_gap_orbital`/`eil_gap_file`）で表面を着色します。渦/表面ソルバーが使う固定 $k_z$ 断面ではなく, 実際の 3 次元 Wannier フェルミ面（全シート・全 $k_z$）上でギャップ（符号・ノード・異方性）を手早く確認できます。

== option=4: スペクトル関数 (`SPECTRUM`)

一体グリーン関数のトレースの虚部から電子のスペクトル関数

$ A(bold(k), omega) = -1/pi "Im" "Tr" G(bold(k), omega + i delta) $

を計算し, symmetry line に沿って描画します。`sw_self=True` の場合は FLEX 計算で求めた自己エネルギーを取り込み, 相互作用の効果を含めたスペクトルを計算します。

== option=5: ボルツマン伝導係数 (`CONDUCTIVITY_BT`)

緩和時間近似のもとでボルツマン輸送方程式から輸送係数を計算します。緩和時間のモデルは `tau_mode` で選択します（既定は定数 $tau$）。

熱電関連の出力：
- 電気伝導率テンソル $sigma_(i j)$ （単位: S/m）
- 熱伝導率テンソル $kappa_(i j)$ （単位: W/m/K）
- ゼーベック係数テンソル $S_(i j)$ （単位: V/K）
- パワーファクター $sigma S^2$ （単位: $upright(W \/ m \/K^2)$）
- ペルチェ係数, ローレンツ数

=== 弱磁場ホール応答

同じモードで $B$ に線形なホール応答も出力します。ホールカーネルは *磁場3方向すべて* についてバンド分解して評価します。

$ S_(a b)^((g)) = 1/2 sum_(bold(k),n) w_bold(k) tau^2 (-partial f \/ partial epsilon) epsilon_(g u v) [v_a (M^(-1))_(b u) v_v - (a <-> b)] $

ホール係数はテンソルのまま構成します。

$ rho^((1)) = -sigma^(-1) sigma^((1)) sigma^(-1), quad R_H^((g)) = rho^((1))_(b a) \/ B quad (a\,b\,g " は巡回") $

$sigma_(x y)\/(sigma_(x x) sigma_(y y))$ で代用すると $sigma$ が $x y z$ で対角であることを仮定することになり, 軸が直交しない格子では成立しません（単斜晶では $B=0$ でも $sigma_(x y) eq.not 0$）。$sigma$ が特異になる場合（厳密な2次元バンドでは $sigma_(z z)=0$）は, その軸に垂直な面内 $2 times 2$ ブロックで評価し, 残り2方向は未定義として出力します。

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

各メッシュについて2つの診断値を出力します。
- $N_"eff"$ — $-partial f\/partial epsilon$ 重みの参加率, すなわちフェルミ窓を担う実効的な状態数。全ての量が $N_"eff"$ 個の平均なので, これがノイズフロアを決めます（$10^3$ 未満で警告）。
- ホール被積分関数の相殺率 $|sum| \/ sum |dot|$。これが小さいときは大きな正負の寄与のわずかな残差になっており（補償近傍）, 符号がまだメッシュ依存の可能性があります（$10^(-2)$ 未満で警告）。

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

- *透過色*: 厚さ $d$ の板の透過率 $T = (1-R)^2 e^(-alpha d) \/ (1 - R^2 e^(-2 alpha d))$（両面フレネル + 非干渉多重反射）。$alpha = 0$ で $T = (1-R)\/(1+R)$ に帰着します。**この色は厚さ $d$ に依存する示量的な量**で, 金属光沢のような物質固有値ではありません。
- *体色（粉体・顔料）*: Kubelka--Munk 二流束模型の拡散反射率 $R_infinity = 1 + K\/S - sqrt((K\/S)^2 + 2 K\/S)$, $K\/S = alpha l$。実効散乱長 $l$ が唯一の自由パラメータです（$alpha$ と Kubelka--Munk の $K$ の規約差と粒径を吸収）。

さらに吸収率 $A = 1 - R - T$ からの色も出力しますが, $A$ は試料から出ていく光の補数なので**ネガ（`hex_neg`, sRGB で $1 - "rgb"$）を取って表示**します。これらは `d`, `l` を設定パラメータ `color_thick`, `color_scatt` で指定します。金属は不透明なので $T$, $R_infinity$ は黒く出ます（正しい挙動で, 図の色見本からも自動的に除かれます）。

出力は `optical.dat`（$sigma$, $epsilon$, $R$, $alpha$, $T$, $R_infinity$ のスペクトルとヘッダに色情報）, `dielectric.pdf/png`, `reflectivity.pdf/png`（色見本付き）です。可視域 380--780 nm (1.59--3.26 eV) を含むように `Emax` を 3.3 eV 以上に設定してください（不足する場合は警告が出ます）。

色の精度について: 金属の色は滑らかな Drude 応答が支配するので誤差が色相に効きにくい一方, 非金属の色相は*吸収端の位置そのもの*です。模型のバンドギャップの誤差はそのまま色相のずれになり, 励起子効果とフォノン介在の間接遷移は本計算（Kubo バブル）に含まれません。`delta` のブロードニングも吸収端を鈍らせるため, 収束と模型の妥当性を確認して使用してください。

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

初期ギャップ振幅 `delta0` の指定方法は 2 通りあります（詳細は後述）：
- 単一の実数（例: `delta0=1.e-2`）: 全バンド共通の最大振幅としてスケールした単一バンドギャップを使用。
- バンドごとの値のリスト（例: `delta0=[0.,0.2,0.3,-0.1,0.]`, 長さ `Norb`）: バンド毎に振幅・符号（multi-gap, 符号反転を含む sign-changing multi-band gap）を指定可能。マルチギャップ系（例: 鉄系超伝導体）の符号反転 $s^plus.minus$ 状態などを扱う場合に使用します。

== option=13: 超伝導状態の特定 $bold(q)$ 点スピン感受率 (`CHIS_QPOINT_SC`)

option=12 と同じ SC 状態の枠組みで, `at_point` で指定した $bold(q)$ 点での $chi_s^"SC"(bold(q), omega)$ の周波数依存性を計算・プロットします。`delta0` の指定方法は option=12 と共通です。

実行時には以下も合わせて出力されます。

- $bold(k) = (0,0,0) arrow.r (0,0.5,0)$ 方向の BdG（Bogoliubov–de Gennes）バンド分散（`BdG_band.png`）
- $chi_s^"SC"(omega)$ の全軌道トレース（`chisq.png`, `chis_sc.dat`）
- 軌道分解した $chi_s^"SC"(omega)$（`chisq_orb.png`, `chis_scorb.dat`）

== option=14: FLEX 自己エネルギー計算 (`FLEX`)

FLEX アルゴリズムを自己無撞着に解き, 電子の自己エネルギー $Sigma(bold(k), i omega_n)$ を松原周波数空間で計算します。得られた自己エネルギーは `sigma.bin` および `self_en.npz` に保存されます。

- `sw_out_self=True`: 自己エネルギーをファイルへ出力
- `sw_in_self=True`: 前回の自己エネルギー (`sigma.bin`) を初期値として読み込む
- `sigma_in_scale`: `sw_in_self=True` で読み込んだ $Sigma$ seed に掛ける係数（既定 1.0）。$U$ アニーリング（小さい $U$ で収束させた `sigma.bin` を seed にして $U$ を段階的に上げる継続法）では, 弱結合の $Sigma prop U^2$ に基づく $(U_"new"\/U_"old")^2$ が目安です。温度アニーリングと異なり松原グリッドが変わらないため, seed の補間が不要という利点があります（`sigma.bin` にはヘッダが無いので `Nw` は前回と同一にしてください）
- 温度アニーリング: 同一 `Nw` での小刻みな冷却（$T_"new"\/T_"old" gt.eq 0.8$ 程度）なら `sigma.bin` の生読み込みで十分です。生読み込みは周波数軸を $T_"old"\/T_"new"$ 倍に圧縮した seed に相当し, この圧縮が冷却に伴う $Sigma$ の増大を近似的に模倣するため, 近臨界では忠実な補間よりむしろ収束が速いことを確認しています。`Nw` を変える場合（生読み込みでは読めない唯一のケース）は Python ヘルパー `plibs.regrid_sigma_bin(temp_new=..., Nw_new=..., w_scale=...)` で `sigma.bin` を新グリッドへ再補間できます。元の温度・`Nw`・軌道数は `sigma.bin` 末尾のメタデータレコード（最近の `io_sigma` が自動追記; 旧コードの読み込みには影響しない末尾レコード）から自動取得され, 明示指定した値がフッタと矛盾する場合はエラーになります（フッタの無い旧ファイルでは `temp_old`/`Norb`/`Nw_old` の明示指定が必要）。`sw_in_self` の読み込み側でもフッタとの照合を行い, メッシュ（`Nw`/軌道数/k点数）が一致しない `sigma.bin` はエラー停止します（Fortran の read はレコードが大きい分には黙って壊れた並びを読んでしまうため, その事故を防ぐ検査です）。温度だけが異なる場合は T アニーリングとみなして情報表示のみで続行します。なお SOC 版（`mkself_soc`）の `sigma.bin` も同一形式に統一済みで, 本ヘルパーがそのまま使えます（旧・要素別レコード形式のファイルは読めません）（共役対称性 $Sigma_(l m)(-i omega) = Sigma_(m l)(i omega)^*$ による負周波数拡張つき 3 次スプライン; 循環畳み込みの折返しアーチファクトを持つ最終周波数点は除去し, カットオフ超過分は内側窓でフィットした $c_0 + c_1\/(i omega)$ tail で補完）。近臨界の冷却では `w_scale=temp_old/temp_new`（生読み込みと同じ圧縮）を推奨, 忠実な補間は `w_scale=1`（既定）です
- `m_diis_num`: DIIS 履歴長（指定なしのデフォルトは 5）
- `sw_rescale_flex=True`: 反復中に Stoner 因子 $S = max_bold(q) "eig"[chi_0(bold(q), 0) hat(S)]$ が 1 に達した場合, $chi_0(bold(q), i nu_n)$ を一律に $(1 - 10^(-4))\/S$ 倍して磁気不安定性による発散を回避
- `sw_chi0_tail=True`: $chi_0$ を tail 補正付きで評価（既定 False; option=14,15,23 の SOC なし経路で有効）。$chi_0 = "conv"[G] - "conv"[G_0] + chi_0^"ref"$ と双線形分解し, ゆっくり減衰する $1\/(i omega)$ 参照部分（$G_0$ バブル）は虚時間の解析式で厳密に, FFT 畳み込みには速く減衰する残差（$G - G_0 tilde 1\/omega^2$）のみを通します。松原打ち切り誤差が $O(1\/N_w) -> O(1\/N_w^2)$ に改善します（検証済み: 厳密 Lindhard 関数に対し収束次数 1.0 → 2.0）。$chi_0$ 段のコストは約 3 倍。$N_w gt.eq 64$ 程度で有効（極端に小さい $N_w$ では参照バブルの $tau$ 離散化誤差の係数が大きく逆転し得ます）。改善された $chi_0$ は $chi_s$/$chi_c$, FLEX 頂点 $V_sigma$, ペアリング頂点 $V_Delta$ に自動的に伝播します。

*実装上の近似に関する注意*（いずれも計算コスト・数値安定性を考慮した意図的な仕様です）:

- `sw_rescale_flex` のリスケールは発散する特定モードだけでなく全運動量・全周波数の $chi_0$ を一様に縮めます。ただしこれは序盤の反復を生き延びるための一時的な処置で, 自己エネルギーのフィードバックにより通常は収束時点で Stoner 因子が 1 未満に下がり, リスケールは非アクティブになります（この場合の収束解はリスケール無しの正真の FLEX 固定点であり, バイアスはありません）。収束した最終反復でもリスケールが作動したままの場合は, その温度・相互作用で正常状態が磁気的に不安定（FLEX の範囲で $T < T_N$）であることを意味し, 標準出力に `[FLEX] WARNING: Stoner factor >= 1 at the final iteration` と表示されます。この警告が出た結果は人工的に弱めた相互作用に対するものなので, 定量的には使用しないでください（温度を上げるか, 高温で収束させた $Sigma$ を `sw_in_self` で読み込む継続法を推奨）。なおクランプ先が 1 の直下（$1 - 10^(-4)$）なのは意図的です: クランプ直後の $Sigma$ の過大応答が, リスケールされた偽の自己無撞着状態から抜け出す機構として働きます。クランプを 0.98 程度に緩めると, $S < 1$ の正しい固定点が存在する近臨界の系でも偽のリスケール固定点に速く「収束」してしまうことを確認しています（振動による反復数の増加は, この安全性の対価です）。
- 自己エネルギーからは静的成分が各反復で減算されます（ラッパー引数 `sub_sigma`; `1`: $Sigma(bold(k), i omega_0)$ の Hermite 部分（デフォルト）, `2`: 周波数平均, `0`: 減算なし）。これはバンド分散の静的シフトを化学ポテンシャル調整側に吸収させてフェルミ面を安定させるための処方であり, 厳密な Hartree--Fock 項のみの減算ではありません（$i omega_0$ での値には低エネルギーの動的成分も含まれます）。絶対的なバンドシフトを議論したい場合は `sub_sigma=0` を使用してください。
- 松原周波数方向の畳み込み（$chi_0$, $Sigma$, Eliashberg カーネル）は長さ $2 N_w$ の FFT による循環畳み込みで, 既定では高周波 tail 補正を行いません（$V$ の最遠点のみ bare 頂点で近似）。実効カットオフは $omega_c = (2 N_w - 1) pi T$ と温度に比例して縮むため, $N_w$ 固定の温度走査には $O(1\/N_w)$ の系統誤差が乗ります。最低温度でも $omega_c$ がバンド幅と $U$ を十分上回るように $N_w$ を選んでください。`sw_chi0_tail=True` で $chi_0$ の打ち切り誤差を $O(1\/N_w^2)$ に改善できます。なお $Sigma$ 自身の畳み込みの tail 誤差は $omega$ に依らない定数（HF 的静的シフト）が主で `sub_sigma`/化学ポテンシャル調整に吸収され, Eliashberg カーネルの tail 誤差はギャップの一様（$s$ 波）成分にのみ入り符号反転ギャップでは対称性で消えるため, これらには個別の補正を実装していません（$V_Delta$, $V_sigma$ は補正済み $chi_0$ を通じて改善されます）。

== option=15: 線形 Eliashberg 方程式 (`LIN_ELIASHBERG`)

RPA または FLEX で求めた有効ペアリング相互作用を使い, 線形化された Eliashberg 方程式を固有値問題として解きます。最大固有値 $lambda$ と対応するギャップ関数 $Delta(bold(k), i omega_n)$ が得られます。$lambda = 1$ で超伝導転移温度 $T_c$ となります。デフォルトは shift + deflation 付き power method, `arnoldi_m > 0` で Arnoldi 法も選択可能です。

- `sw_self=False`: 自己エネルギー補正なし（RPA）
- `sw_self=True`: FLEX の自己エネルギーを用いる（`sigma.bin` が必要）
- `sw_from_file=True`: `sigma.bin` からファイルを読み込み, FLEX 計算をスキップ
- `gap_sym`: 初期ギャップ関数の対称性を指定（詳細は後述）
- ギャップ関数は `gap.npy` と, 軌道ペアごとの `gap_11.dat`, `gap_12.dat`, …（軌道インデックスを 1 始まりにした `gap_{i j}.dat`）に出力される

== option=16: ギャップ関数の後処理 (`GAP_FUNCTION`)

option=15 で計算したギャップ関数 (`gap.npy`) を読み込み, 以下の量を計算・出力します。

- ギャップ関数の異常グリーン関数 $F(bold(k), i omega_n)$
- Padé 近似による実周波数への解析接続スペクトル

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

幾何と単位変換は, 正方晶格子での $F(theta)=F(0)/cos(theta)$ を含む解析的な円筒フェルミ面テストで検証しています。物質モデルでの確認として, 実験格子定数を用いた 10 軌道 FeS Wannier 模型では $B parallel c$ で $0.489$–$3.123$ kT となり, T. Terashima _et al._, Phys. Rev. B *99*, 134501 (2019) の Fig. 3(b) に示された同じ磁場方向の GGA 振動数（およそ $0.5$–$3$ kT）と整合します。同論文の実験値は $0.15$–$1.87$ kT とこれより小さいですが, この差は計算の不備ではなく想定されるものであり, 鉄系に限った話でもありません。強相関金属では DFT（LDA/GGA）のバンド構造が電子相関によるバンドの narrowing と軌道依存のバンドシフトを含まないため, 一般にフェルミ面のポケットを系統的に大きく, 質量繰り込み（したがって $E_F$ での状態密度）を小さく見積もります。この偏りは本コードの対相関の計算にもそのまま伝播します。繰り込みを受けていない DFT バンドの上で評価したクーパー不安定性は弱く出るため, RPA/Eliashberg の固有値や $T_c$ は過小評価されます。鉄系以外ではずれがさらに大きくなる例もあります。重い電子系反強磁性体 $"CeRhIn"_5$ では, $4f$ 電子を遍歴として扱うバンド計算は large Fermi surface を与えますが, dHvA で観測されるフェルミ面は $4f$ を持たない参照物質 $"LaRhIn"_5$ とほぼ同じで, サイクロトロン質量だけが約 1 桁重くなります（H. Shishido _et al._, J. Phys. Soc. Jpn. *71*, 162 (2002)）。$"Na"_x "CoO"_2$ では LDA が予言する 6 個の小さな $e'_g$ ホールポケット（D. J. Singh, Phys. Rev. B *61*, 13397 (2000)）が光電子分光では観測されず, $a_(1g)$ 由来のシートのみが見えます（H.-B. Yang _et al._, Phys. Rev. Lett. *92*, 246403 (2004)）。これを局所相関だけで説明できるかについては, 現在も議論が続いています。したがってこの option のベンチマークは同種のバンド構造から計算された振動数であり, 実験値そのものではありません。実験を再現するには, 相関を取り込んだバンド構造（QSGW, LDA+DMFT）や経験的なバンドシフトに基づく Wannier 模型が必要です。

== option=20: 電子質量 (`ELECTRON_MASS`)

symmetry line 上で $m^* = hbar^2 (partial^2 E \/ partial bold(k)^2)^(-1)$ を計算します（$m_e$ 単位）。

== option=21: 不純物入りスペクトル (`SPECTRUM_IMPURITY`)

実空間スーパーセル上で不純物を加えたハミルトニアンを構築し, 不純物グリーン関数を介してスペクトル関数 $A(bold(k), omega)$ を計算します。`imp_list` で不純物サイトを指定。

== option=22: CPA 伝導 / スペクトル (`SIGMA_CPA`)

Coherent Potential Approximation により混晶の自己エネルギーを反復解で求め, 実周波数のスペクトル関数（`cpa_spectrum.png`）と松原周波数の自己エネルギー（`sigma.bin`, `self_en.npz`）を出力します。`x_cpa`（不純物濃度）, `VA`, `VB`（オンサイト摂動）はサブルーチン内で固定値が設定されています。

== option=23: 非線形 Eliashberg 方程式 (`NONLIN_ELIASHBERG`)

非線形（自己無撞着）SC FLEX-Eliashberg ループを解きます。線形 Eliashberg と異なり Δ を $T_c$ 以下で有限振幅まで成長させ, 異常グリーン関数 $F$ を含めた SC Dyson 方程式と FLEX 自己エネルギーを反復します。

初期ギャップは内部で自動生成されます（`gap.npy` を事前に用意する必要はありません）。

1. まず線形化 Eliashberg 方程式を内部で解き, Stoner 因子 $S$ と最大固有値 $lambda_"eliash"$ を計算します。
   - $S >= 1$（SDW/CDW 不安定）の場合は非線形ループに入らず終了します。
   - $lambda_"eliash" < 1$（$T >= T_c$, 超伝導不安定性なし）の場合も終了します。
   - `sw_check_only=True` の場合はここで終了し, $S$ と $lambda_"eliash"$ のみを報告します（非線形ループをスキップして温度走査で $T_c$ を素早く確認する用途）。
2. 上記の線形固有ベクトルから対称性が正しい初期 Δ の形を取り, 振幅を BCS 弱結合比 $Delta_0 = 1.764 k_upright(B) T_c$ にスケールして非線形ループに渡します。
3. 各反復: $Delta_"new" = T \/ N_k sum V_Delta dot F$ → 振幅方向 Newton（secant）加速 + DIIS 形状混合（バイパス時は線形混合 pp=0.3）→ $Sigma$ 更新 → SC Dyson で $G, F$ 更新 → $chi^"GG", chi^"FF"$ から $V_sigma, V_Delta$ 再構築

- `m_diis_num`: DIIS 履歴長（$>= 2$ で Pulay 外挿, それ以下は線形混合）
- `sw_self=True`: FLEX 自己エネルギー込みで Σ-dressed Green's function を使う
- `sw_from_file=True`: `sigma.bin` から自己エネルギーを読み込む
- `sw_check_only`: 後述（5節「スイッチ類」）
- 振幅方向の Newton 加速によりギャップ振幅の収束が高速化されています（詳細は `libs/src/ffeliash.f90` の `sw_amp_newton` 関連コメント参照）
- 詳細は `libs/src/ffeliash.f90` 参照

== 準古典 Eilenberger モード（option=24,25,26）

option=24〜26 は超伝導の*準古典 Eilenberger 方程式*（エネルギー積分された Gor'kov 方程式）を解きます。ギャップや不純物がフェルミ波長スケールでゆっくり変化する（$Delta, T_c, hbar\/tau << E_F$）状況で有効で, $T_c$・状態密度・表面 Andreev 束縛状態・渦／渦格子状態を扱う自然な枠組みです。未知数は準古典伝播関数 $g(bold(k)_F, bold(r), omega_n)$（normal）と $f$（anomalous）で, 単一の*Riccati 振幅* $a$（$f = 2a\/(1+a a^*)$, $g = (1-a a^*)\/(1+a a^*)$）で表し, フェルミ速度に沿った直線軌道上を数値的に安定な Fortran カーネル（`riccati_chords`; d-vector 用には $2 times 2$ スピン版 `matrix_riccati_chords`）で積分します。ペアリングは可分離形 $Delta(bold(k)_F, bold(r)) = Delta(bold(r)) phi(bold(k)_F)$ で, 形状因子 $phi$ は `gap_sym` で決まり $ |phi|^2 _"FS" = 1$ に規格化されるため, 結合定数 `eil_coupling` は無次元の $lambda$ です。

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

=== option=26: 渦・渦格子 (`EILENBERGER_VORTEX`)

磁束渦まわりの非一様ソルバです。`eil_field` $= B\/H_(c 2)$ で形状を選びます: `0` は大きな円形セル内の孤立渦（半径 `eil_vort_lxi`（$xi$ 単位）, グリッド `eil_vort_ngrid`）, `>0` は Doppler シフト付き円形セル渦格子。ギャップ profile $Delta(rho)$ と（`eil_ldos` で）ゼロエネルギーコア LDOS（Caroli–de Gennes–Matricon 束縛状態; 格子では Volovik の $sqrt(B)$ DOS）を出力します。Zeeman 場 `eil_zeeman` はコア状態をスピン分裂させます。サブモード:

- `eil_vort_current=True`: 環状超流動電流 $j_phi(rho)$ → `vortex_current.dat`。
- `eil_vort_field=True` / `eil_vort_maxwell=True`: 自己無撞着な有限 $kappa$ 磁場 $B(rho)$ / ベクトルポテンシャル $bold(A)(bold(r))$（Maxwell 反作用; `eil_kappa` を使用）。
- `eil_vort_lattice_sc=True` ＋ `eil_field_list`: *真の周期*磁気 Bloch 渦格子（定式化 A, 極端第 II 種）: 各コアに実ノードを持つ複素秩序変数 $Psi(bold(r))$ と完全な Abrikosov 超流動位相, $B\/H_(c 2)$ を走査して $ N(0)  (B)$ を得る（$d$ 波 $tilde sqrt(B)$ Volovik）。`eil_lattice` は `'square'`/`'triangular'`, `eil_nvortex` はセルあたり磁束量子数, 有限 `eil_kappa` で London 遮蔽を追加（`eil_vort_scA=True` で $bold(A)$ を準古典電流 $bold(j)_s =  bold(v)_F "Im" g $ から完全自己無撞着に, `je` の `A_renew` 法）。
- `eil_vort_dvector=True`: 自己無撞着トリプレット d-vector 渦／格子テクスチャ（主成分の winding ＋コア局在副成分; スピン行列 Riccati）。
- `eil_gap_orbital`: 軌道基底のペアポテンシャルで, その FS バンドへの低エネルギー射影がギャップを決める（Nagai–Nakamura マルチバンド Eilenberger, JPSJ *85*, 074707 (2016), Eq. 43; Wannier FS が必要）, `gap_sym`/`delta0` に優先します。

付随ドライバ `calc_vortex_lattice_symmetry`（ライブラリから呼び出し）は Ichioka–Machida 格子自由エネルギーをセルの頂角とギャップ対格子方位 $theta_0$ について最小化し,*安定な渦格子対称性*とその磁場発展（例: $H_(c 2)$ 近傍での $d$ 波 三角→正方 転移）を決定します。Wannier FS を与えると $theta_0$ は結晶全体（FS ＋ギャップ）を剛体回転させるため, フェルミ速度異方性も対称性選択に寄与します。

== option=27: 超伝導状態の NMR (`NMR_SC`)

option=12,13 と同じ BdG スピン帯磁率から, 超伝導状態の代表的な NMR 観測量である Knight shift とスピン格子緩和率を計算します。

$ K(T) prop chi'_s (bold(q)=0, omega -> 0), quad 1\/(T_1 T) prop angle.l |A(bold(q))|^2 chi''_s (bold(q), omega_0) angle.r_bold(q) \/ omega_0 . $

同じ温度・$bold(k)$ メッシュ・$bold(q)$ メッシュ・スメアリングで *正常状態* も評価され, 物理的に意味を持つのはその比です。$K_s\/K_n$ はシングレットギャップに対する Yosida 関数, $(1\/T_1T)_s\/(1\/T_1T)_n$ はフルギャップの Hebel--Slichter ピーク, あるいは節を持つギャップのべき則（ライン節なら $T^3$）を与えます。比を取ると物質固有の係数（$gamma_n$, 超微細結合のスケール, $g$）が消え, メッシュ誤差の大部分も相殺するため, 比は各々の生の値よりはるかに速く収束します。`nmr_sc.dat` の生の `K_sc`, `K_n`, `(1/T1T)` 列は診断用と考え, 比を引用してください。正常状態の参照値は正常状態のバブルで計算します（SC バブルの厳密な $Delta -> 0$ 極限であり, 4 倍高速です）。

ギャップは `gap_sym`/`delta0` の形状因子（option=12,13 と同じもの）か, `nmr_gap_file` を指定した場合は option=15/23 が出力した自己無撞着な RPA/FLEX/Eliashberg ギャップです。いずれも $Delta(0)$ として読まれ, 各温度では BCS 内挿 $Delta(T) = Delta(0) tanh(1.74 sqrt(T_c\/T - 1))$ によって *振幅だけ* がスケールされ, 形状は固定されます。$T_c$ の既定値は *フェルミ面上* の $|Delta|$ の最大値から得る弱結合の見積り $max |Delta|_"FS" \/ 1.764$ です（ゾーン全体の最大値は, 状態を担うシートに形状因子の重みが乏しい場合に何倍も大きくなり, 掃引全体が誤った換算温度になってしまいます）。`nmr_tc` で上書きできます。

掃引は `nmr_trange`（$T_c$ 単位）を `nmr_nt` 点で行い, `nmr_sc.dat` と `nmr_sc.png` を出力します。開始前に, 使用するスピンチャンネル, $Delta(0)$ と $T_c$, $bold(k)$ メッシュとスメアリングの階層のチェック（$delta gt.tilde v_F d k$ が必要）, コスト見積り（総コストは $N_T N_q N_k$ に比例）, および `nmr_qconv=True` のときは `nmr_qsize` とその半分のメッシュでの $bold(q)$ 和の変化が出力されます。

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

- `fname` （文字列）: 入力ファイルのパス（`ftype`に応じた形式で指定, 拡張子は通常不要）

- `ftype` （整数）: 入力ハミルトニアンのファイル形式（4節参照）

- `brav` （整数）: ブラベー格子の種類（5節参照）

- `sw_soc` （bool）: スピン軌道相互作用を考慮するかどうかのスイッチ

- `inv_tol` （実数, 既定 $3 times 10^(-2)$）: 反転操作を $H(bold(R))$ の対称性として受け入れる相対残差のしきい値（option=14,15,16,23 は $Delta(bold(k))$ から $Delta(-bold(k))$ を再構成するのにこの操作を使います）。第一原理計算のフィットモデルには数 % の誤差があるため既定値は意図的に緩めてあります。フィットモデルが棄却される場合は大きくしてください。大きな値が必要だったときは, 出力される残差を確認してから結果を使ってください。

== メッシュ設定

- `Nx, Ny, Nz` （整数）: 第一ブリルアンゾーン内の $bold(k)$ 点メッシュ数。フェルミ面積分や FFT-convolution（FLEX, Eliashberg, chi）で使用します。2 次元系では `Nz=1` とします。FLEX/Eliashberg 計算では FFT 効率の都合上 2 のべき乗（32, 64 など）が望ましく, メモリ消費は $tilde.equiv N_x N_y N_z dot N_w dot N_"orb"^2$ で増えます。

- `Nw` （整数）: 松原周波数のメッシュ数（option=12,13,14,15,16,23 の FLEX/Eliashberg/SC chi 計算）または実周波数のメッシュ数（option=1,4 など）。温度 $T$ に対して $omega_n = (2n+1) pi T$ の形で $n=0,1,...,N_w-1$ まで取ります。低温ほど大きな $N_w$ が必要です（典型値: 256–1024）。実効的な周波数カットオフは $omega_c = (2 N_w - 1) pi T$ で温度に比例して縮み, 松原畳み込みは鋭い打ち切り（tail 補正なし）です。したがって `Nw` 固定のまま降温すると系統誤差が増えるため, 温度走査では最低温度でも $omega_c$ がバンド幅・`U` を十分上回るよう `Nw` を選んでください（option=14 の注意書きも参照）。

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

- `delta0` （実数 または 長さ`Norb`のリスト, 単位: eV）: option=12,13 の SC chi 計算で用いる初期ギャップ関数の振幅。物理的には超伝導ギャップサイズに対応します（典型値: $10^(-3) - 10^(-2)$ eV ≒ 1–10 meV）。
  - 実数を指定した場合: 全バンド共通の単一ギャップとして, 内部で生成した対称性形状を最大振幅 `delta0` にスケールします。0 にすると常状態と等価になります。
  - リストを指定した場合（例: `delta0=[0.,0.2,0.3,-0.1,0.]`）: バンド毎に振幅と符号を個別に指定する multi-gap モードになります。符号が異なる成分を混在させることで $s^plus.minus$ 型など符号反転ギャップを表現できます。

== 伝導計算パラメータ

- `tau_const` （実数, 単位: fs）: `tau_mode='const'` のときに使う定数緩和時間 $tau$。定数近似では自由パラメータであり, 実験値との比較から決めます（例: 金属では 1〜100 fs 程度）。なお $R_H$ は $tau^2\/tau^2$ を持つため定数 $tau$ には *依存しません*。`tau_mode='const'` では, ホール係数はフェルミ窓の広がり以外に温度依存性を持たない純粋なバンド構造由来の量になります。

- `tau_mode` （文字列, 既定 `'const'`）: option=5 の緩和時間モデル。
  - `'const'`: 定数 `tau_const`（$bold(k)$・バンドに依らない）。
  - `'epa'`: 電子フォノン平均結合（Samsonidze--Kozinsky の EPA）による $tau(bold(k))$。`epa_file` が必要です。$R_H (T)$ やシートごとの移動度差が意味を持つのはこのモードです。
  - `'dos1'` / `'dos2'`: DOS ベースの $tau$（`tau_const` でスケール）。
  実行時にフェルミ面平均 $angle.l tau angle.r$, フェルミ窓での範囲, および $angle.l tau^2 angle.r \/ angle.l tau angle.r^2$ が出力されます。最後の量は定数 $tau$ でちょうど 1 になり, $tau(bold(k))$ が $R_H$ を定数 $tau$ の値からどれだけずらすかを表します。

- `epa_file` （文字列 または `None`）: `tau_mode='epa'` のときに読み込む `epa.x`（`job='egrid'`）の出力ファイルパス。

- `hall_mesh_list` （リスト または `None`, 既定 `None`）: option=5 のホール係数に対する $bold(k)$ メッシュ収束スイープ。例: `[20,30,40]` や `[(20,20,10),(30,30,15)]`。各要素は整数（立方メッシュ）か `(Nx,Ny,Nz)` タプルです。`None` の場合は `(Nx,Ny,Nz)` 単一メッシュのみで, 収束情報は得られません。粗いメッシュでは $sigma^((1))$ が符号ごと誤り得るため, 新しい物質では最低一度はスイープすることを強く推奨します。

- `sw_tdf` （bool）: `True` にすると輸送分布関数（TDF）を先に計算し, そのエネルギー積分から輸送係数を求めます（エネルギー依存緩和時間に対応する場合に有効）。

伝導モード（option=5,6）のスピン縮退因子は `sw_soc` に従います。`sw_soc=False` ではバンドがスピン縮退しているとして $g_s=2$ を, `sw_soc=True` ではハミルトニアンに両スピンが既に含まれているとして $g_s=1$ を使います。

== 軌道・相互作用パラメータ

- `olist` （リスト）: カラープロット（color_option=1）で使用する軌道インデックスのリスト（7節参照）。

- `U` （実数, 単位: eV）: オンサイトクーロン斥力（ハバード$U$）。FLEX/RPA計算で使用。磁性・超伝導の計算には重要なパラメータです。

- `J` （実数, 単位: eV）: オンサイトフント結合定数。$U' = U - 2J$（金森遮蔽）が自動で計算されます。

- `orb_dep` （bool）: `True`にすると軌道依存の相互作用行列 `Umat`, `Jmat` を使用します。`False`（デフォルト）では `U`, `J` の一定値を全軌道に適用します。

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

$cos(2 pi k_x)$ などの簡約座標で書いた形では, これらの調和関数がラベル通りの対称性を持つのは正方晶・直交系の 1 サイト格子だけです。fcc・bcc・六方晶のセルでは $cos(2 pi k_x) - cos(2 pi k_y)$ は周期的ではあっても $B_(1g)$ *ではなく*, 結晶のデカルト $C_(4z)$ で符号が反転しません。そこで形状因子は, 模型の R シェル上の*デカルト格子調和関数*として構成しています。

$ phi_Gamma (bold(k)) = sum_(bold(R) in "shell") K_Gamma (hat(bold(R))) e^(i bold(k) dot bold(R)), quad bold(R)_"cart" = bold(n) dot "avec" $

ここで $K_Gamma$ は各隣接方向について評価したデカルトの調和関数（$hat(n)_x^2 - hat(n)_y^2$, $hat(n)_x hat(n)_y$, $hat(n)_x$ など）です。*軸の分解は一切不要*である点に注意してください。位相 $bold(k) dot bold(R) = 2 pi bold(k)_"frac" dot bold(n)$ は基底に依らないため, デカルト幾何が必要なのはシェルの*係数*だけです。したがってギャップのために `sw_dec_axis` を使う必要はなく, BZ 和や FFT が前提とする $bold(k)$ メッシュの周期性もそのまま保たれます（$bold(R)$ は真の格子ベクトルなので $phi(bold(k)+bold(G)) = phi(bold(k))$ が構成上成立します）。直交系の 1 サイト格子では, 上の表の全ての指定・任意の軸比について従来の簡約座標の式をビット単位で再現するので, 既存の入力に影響はありません。

`main.py` はこのための格子を自動的に登録します（`plibs.set_harmonic_lattice(avec, rvec)`）。`plibs.lattice_harmonic()` で単独の調和関数を直接構成でき, `gap_symms(..., avec=, rvec=)` で格子を明示的に渡すこともできます。模型の $bold(R)$ 集合が指定した対称性を担えない場合は, 警告を出して従来の簡約座標の形状因子にフォールバックします。

注意点が 2 つあります。カイラルな組は 2 つの相方で同じシェルを共有する必要があります（そうしないと, 両者を混ぜる回転に対して $|phi|$ が不変になりません。六方晶の $E_(2g)$ は相方ごとに別シェルを使うと 77% ずれます）。これはコード側で自動的に処理します。もう一つは*多サイト*セルで, ペアリングのボンドは $bold(R)$ ではなく $bold(R) + bold(tau)_j - bold(tau)_i$ なので, サブ格子の位相が依然として欠けます。鉄系の 2-Fe セルのように反転がサイトを入れ替える模型では, 軌道基底のルート（`eil_gap_orbital` / `eil_gap_file` と Nagai–Nakamura のバンド射影）を使ってください。

== $bold(k)$空間の設定

- `kz` （実数, 規格化単位）: option=2のフェルミ面プロットにおける$k_z$の値。$0 \sim 0.5$の範囲で指定します（$k_z = 0$は$Gamma$点面, $k_z = 0.5$はゾーン境界面）。

- `kscale` （リスト or 実数）: option=3の3次元フェルミ面表示における各軸の表示スケール。例: `kscale=[1.0, 1.0, 0.5]`とすると$k_z$方向を半分に縮めて表示。

- `k_sets` （リストのリスト）: symmetry lineの各端点のk座標（規格化単位, $0 \sim 1$）。自動生成ではなく手動で対称線を設定したい場合に定義します。例: `k_sets=[[0,0,0],[0.5,0,0],[0.5,0.5,0]]`

- `xlabel` （文字列のリスト）: `k_sets`に対応するラベル。例: `xlabel=[r'$\Gamma$','X','M']`

- `at_point` （リスト）: option=8で感受率を計算する$bold(q)$点の座標（規格化単位）。

== スイッチ類

- `sw_unit` （bool）: `True`（デフォルト）で国際単位系（SI）の物理定数を使用し, 出力も実際の物理単位で得られます。`False`にすると$hbar = k_upright(B) = e = 1$の無次元系になります。

- `sw_omega` （bool）: option=11 でペアリング感受率を実周波数（`True`）または松原周波数（`False`）で計算するかのスイッチ。

- `sw_self` （bool）: option=4 (スペクトル) / option=15 (線形 Eliashberg) / option=23 (非線形 Eliashberg) で FLEX の自己エネルギーを使うかどうか。`True` にすると繰り込まれたグリーン関数を用います。

- `sw_out_self` （bool）: `True` にすると FLEX 計算の自己エネルギーを `sigma.bin` および `self_en.npz` に出力します。

- `sw_in_self` （bool）: `True` にすると FLEX 計算の初期自己エネルギーとして前回の `sigma.bin` を読み込みます（自己無撞着ループの収束を早める）。

- `sigma_in_scale` （float, 既定 1.0）: `sw_in_self=True` で読み込んだ $Sigma$ seed に掛ける係数。$U$ アニーリング用で, $(U_"new"\/U_"old")^2$ が目安（option=14 の説明参照）。

- `sw_from_file` （bool）: `True` にすると自己エネルギーを `sigma.bin` から読み込み, FLEX 計算を実行せずに直接 Eliashberg 計算を行います。

- `sw_check_only` （bool）: option=23（非線形 Eliashberg）専用。`True` にすると内部の線形化 Eliashberg 解（Stoner 因子 $S$ と固有値 $lambda_"eliash"$）を計算した時点で停止し, 非線形ループを実行しません。温度走査で $T_c$（$lambda_"eliash"=1$ となる温度）の見積もりだけを高速に行いたい場合に使用します。$S >= 1$ または $lambda_"eliash" < 1$ の場合は, このフラグの値に関わらず非線形ループ手前で自動的に停止します。

- `sw_rescale_flex` （bool）: option=14 の FLEX 計算で, 反復中に Stoner 因子が 1 に達した場合に $chi_0$ 全体を一律に $(1-10^(-4))\/S$ 倍して発散を防ぎます。反復の序盤を安定化するための処置で, 通常は $Sigma$ のフィードバックにより収束時にはリスケールが非アクティブになり, 収束解にバイアスは残りません。最終反復でもリスケールが作動していた場合は正常状態が磁気不安定である旨の警告（`[FLEX] WARNING`）が表示されます（option=14 の注意書き参照）。

- `sw_chi0_tail` （bool）: option=14,15,23（SOC なし経路）で $chi_0$ を tail 補正付きで評価します（既定 False）。松原打ち切り誤差が $O(1\/N_w)$ から $O(1\/N_w^2)$ に改善し, 同じ精度なら小さい `Nw` で済みます（$chi_0$ 段のコスト約 3 倍, $N_w gt.eq 64$ 目安）。詳細は option=14 の注意書きを参照。

- `sw_dec_axis` （bool, *フェルミ面プロット専用*）: True にすると格子ベクトルを慣用的な直交軸で表し直し, $bold(R)$ も合わせて変換します。これは*描画*のための機能で, 中心格子（とりわけ体心の 122 セル）を見慣れた箱に展開（unfold）し, 10 軌道（2-Fe）のフェルミ面を 5 軌道（1-Fe）のものと直接比較できるようにするためのものです。*物理的な* $bold(k)$ における $E(bold(k))$ は厳密に保存されます（全ての `brav` について $10^(-14)$ で確認）。

  `brav=1,2,6,7` では変換後の $bold(R)$ が半整数成分を持ちます。これは unfold そのものですが, 同時に $H(bold(k))$ が分解後の簡約座標 $bold(k)$ について周期的でなくなることを意味します。したがって BZ 和・FFT 畳み込み・既約 wedge の折り畳み・セル体積（$det("avec")$ は化学式単位あたりの体積ではなくなります）はいずれもこのスイッチ下では無効です。そのため描画系の option=0（バンド図）, option=2, option=3（フェルミ面）でのみ受け付け, それ以外のモードはエラーで停止します。バンド図が許可されるのは経路上で $E(bold(k))$ を評価するだけだからで, DOS（option=1）は BZ 和なので対象外です。option=0 で対称線を自動生成している場合は, 分解後（単純直交）の軸に合わせた経路を生成し直します（中心格子の経路は分解後の軸では別の物理的 $bold(k)$ を指すため）。`k_sets` を自分で指定した場合はそのまま使われるので, 分解後の軸で与える必要があります。六方晶・三方晶（$sqrt(3)\/2$ 成分）は, そもそも整数の直交格子に載せられません。ギャップの対称性についても本スイッチは不要です（デカルトの $bold(R)$ シェルから構成されます。上の `gap_sym` 参照）。

=== Eilenberger パラメータ（option=24,25,26）

準古典 Eilenberger ソルバを駆動します。温度はグローバルの `tempK`/`temp`, ギャップ対称性はグローバルの `gap_sym` を用います。

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
- `eil_zeeman` （float, eV）: LDOS 用 Zeeman（Maki）場（表面: $d_[110]$ ZEBS を $plus.minus h$ に分裂; 渦: コア状態をスピン分裂）。

*一様（option=24）:*

- `eil_method` （`'normalization'`/`'riccati'`）: $(g, f)$ ルート — `'normalization'` は高速, `'riccati'` は非一様カーネルと整合。
- `eil_find_tc` （bool）: 現在の不純物設定で $T_c$ を二分探索。
- `eil_imp_sweep` （bool）, `eil_imp_list` （array）: `eil_imp_list` で $Gamma$ を走査し $T_c(Gamma)$ を `eilenberger_tc.dat` に出力。
- `eil_pauli`, `eil_spin`, `eil_lambda`, `eil_fs`, `eil_free_energy` （bool）: 上記 option=24 で述べた排他的サブモード。

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
- `eil_gap_file` （`None` / 文字列）: option=15/23（`LIN_ELIASHBERG`/`NONLIN_ELIASHBERG`）を `sw_out_self=True` で実行したとき `output_gap_wannier` が Wannier 実空間「ホッピング」形式で書き出した自己無撞着 RPA/FLEX ギャップの基底名（拡張子なし, 例 `'gap_wannier'`）。指定すると $Delta(bold(R), i omega_n)$ を読み込み, その逆フーリエ変換 $Delta_"orb"(bold(k)) = sum_bold(R) e^(i 2 pi bold(k) dot bold(R)) Delta(bold(R))$ を FS バンドへ射影して `eil_gap_orbital` として用いる。*以前計算した RPA ギャップ*（例: $"KFe"_2"As"_2$, PRB *84*, 144514）を渦のペアリング形状因子に使うための経路。ギャップを書き出した RPA 計算と Eilenberger 計算は*同一*の Wannier Hamiltonian（同じ軌道基底・$bold(R)$/$bold(a)$ 規約, できれば $mu$/filling）を使う必要がある（バンド固有ベクトルと $Delta_"orb"$ が同じ基底であるため）。`eil_gap_orbital`/`gap_sym` に優先。
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

標準出力に Stoner 因子 $S$ と固有値 $lambda_"eliash"$ が表示されます。$lambda_"eliash" > 1$ となる温度まで降温したら, `sw_check_only=False` にして非線形ループを実行します。

```python
option = 23
sw_self = True            # FLEX 自己エネルギー込み
sw_from_file = True
sw_check_only = False
m_diis_num = 5            # DIIS Pulay 加速 + 振幅方向 Newton 加速
```

低温（$T tilde.equiv T_c \/ 5$）では DIIS 履歴を 5–10 程度にすると収束が速くなります。

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

= テストコード

`tests/` ディレクトリには, 数値カーネルと物理ベンチマークの回帰テストが 7 ファイル・173 項目 含まれています。`pytest` がインストールされていない環境でも, 各テストファイルを Python スクリプトとして直接実行できます。

実行前に Fortran 共有ライブラリ `libs/libfmod.so` がコンパイル済みである必要があります。未コンパイルの場合は `libs` ディレクトリで `make FC=<compiler> SL=<library>` を実行してください。

== 実行方法

個別に実行する場合:

```bash
python tests/test_transport.py
python tests/test_rpa_flex.py
```

サマリ表付きで全て実行する場合（`pytest` 不要）:

```bash
python tests/run_all.py          # 全ファイル
python tests/run_all.py rpa      # ファイル名に "rpa" を含むものだけ
```

`pytest` が利用できる環境では `pytest tests` でも実行できます。一部のテストは `tmp_path` フィクスチャを取るため `pytest` でのみ実行され, 単体ランナーではスキップされます。

== ファイル別のカバー範囲

#table(
  columns: (auto, auto, auto, 1fr),
  align: (left, right, right, left),
  table.header([*ファイル*], [*項目数*], [*時間*], [*確認している内容*]),
  [`test_eilenberger.py`], [57], [437 s],
    [準古典 Eilenberger / Riccati ソルバー: 均一系の極限, Fortran カーネルと Python 参照実装の一致, 表面 Andreev 束縛状態, 渦糸および自己無撞着渦糸格子, 凝縮自由エネルギー, model / Wannier フェルミ面, Pauli 制限, triplet $d$-vector texture],
  [`test_rpa_flex.py`], [34], [0.7 s],
    [RPA / FLEX / Eliashberg の構成要素: 相互作用頂点, RPA 行列代数, 厳密 Lindhard 和・対バブルとの照合による $chi_0$ と $phi_0$, テイル補正 $chi_0$, 自己エネルギーの再グリッド, および時間反転対称性下の既約 $bold(k)$ メッシュ],
  [`test_effmass.py`], [29], [3.8 s],
    [直交基底・MLO 基底での逆有効質量（バンド間項, overlap 項, 縮退の扱い）と, dHvA 極値軌道の抽出: 断面の幾何, 開軌道の棄却, ゾーン還元],
  [`test_nmr.py`], [22], [16 s],
    [超伝導状態の NMR: BCS ギャップ内挿, $bold(q)$ メッシュの折り畳み, 自己無撞着 Eliashberg ギャップの読み込み, Knight shift（Yosida と triplet）, Hebel--Slichter ピーク, ライン節の $T^3$ 則],
  [`test_transport.py`], [18], [0.8 s],
    [直流極限での Kubo$arrow.l.r$ボルツマン一致, Wiedemann--Franz 則, 対称化されたバンド間熱流頂点, バンド縮退, および弱磁場ホール応答: スカラー参照との $sigma^((1))$ 一致, $R_H arrow -1\/(n e)$ と格子角非依存性, Luttinger 体積によるキャリア密度],
  [`test_velocity_mlo.py`], [7], [0.4 s],
    [非直交（MLO）基底でのバンド速度: $-epsilon_n C^dagger (partial S\/partial k) C$ の overlap 項, エルミート性, 直交経路への完全一致でのディスパッチ],
  [`test_tools.py`], [6], [0.5 s],
    [他ファイルが照合先とする参照実装: Fermi / Bose 恒等式, 閉形式のモデル束, Fortran 畳み込みの検証に使う厳密 $chi_0$],
)

スイート全体で約 7.5 分かかり, そのほぼ全て（96%）が `test_eilenberger.py` です（準古典ソルバは 2 次元グリッド上の自己無撞着場の問題のため）。突出した項目はなく, 重い方から 4 つ（`test_vortex_in_plane_field`, `test_chiral_vortex`, `test_surface_gap_heals_to_bulk`, `test_free_energy_selects_the_chiral_partner`）が各 57–77 秒です。さらに軽くしたい場合はこれらの `Ng` / `nbeta` / `ngrid` / `eil_fs_nkz` を下げてください。他のファイルはいずれも 20 秒以内で終わります。

= トラブルシューティング

- *化学ポテンシャルが収束しない*: `Nx, Ny, Nz` を増やすか, `tempK` を少し上げてください。

- *スペクトルにノイズが多い*: `delta` を少し大きくする, または `Nw` を増やしてください。

- *バンドがおかしい*: `brav` の設定が合っているか確認してください。特に QE との連携では `brav=1`（FCC）や `brav=2`（BCC）の convention に注意が必要です。

- *線形 Eliashberg が収束しない*: `U` を小さくする, `tempK` を高くする, メッシュ `Nx,Ny,Nz` を増やすなどを試してください。固有値 $lambda$ が 1 を超えている温度は超伝導相内に対応します。

- *非線形 Eliashberg が発散する*: `m_diis_num` を 5–10 に増やす, `tempK` を上げる, `Nw` を増やすなどを試してください。`sw_check_only=True` で先に Stoner 因子 $S$ と固有値 $lambda_"eliash"$ を確認するのも有効です。$S$ が 1 を超えている場合は磁気不安定なので `U` を下げる必要があり, $lambda_"eliash" < 1$ の場合はまだ $T_c$ 以上なので降温してください（いずれの場合も非線形ループは自動的に実行されません）。

- *FLEX で max$|Sigma|$ が発散*: `sw_rescale_flex=True` を設定するか, `U` を下げてください。

- *伝導率の単位が合わない*: `sw_unit=True` になっているか確認してください。`False` にすると無次元系になります。

- *ホール係数が $bold(k)$ メッシュで符号が変わる*: これは $sigma^((1))$ の典型的な破綻の仕方であり, バグではありません。`hall_mesh_list` でスイープし, 2つの診断値を見てください。$N_"eff"$ が $10^3$ を下回っていればフェルミ窓を担う状態数が足りず, 相殺率が $10^(-2)$ を下回っていれば補償近傍で $sigma_(x x)$ よりはるかに細かいメッシュが必要です。`tempK` を上げてフェルミ窓を広げるのも有効です。

- *ホールキャリア密度が実験と合わない*: 同じ出力に併記される Luttinger 体積と突き合わせてください。$n_e \/ n_h approx 1$ なら補償金属であり, $n_H$ はそもそもキャリア密度ではありません（移動度の比較や二流体フィットを行ってください）。また `sw_soc` が正しく設定されているか（スピン縮退因子を決めます）, 実験磁場での $omega_c tau$ が $lt.tilde 0.1$ に収まっているかも確認してください。

- *$R_H$ に温度依存性が出ない*: `tau_mode='const'` では当然の挙動です。$R_H$ は $tau^2\/tau^2$ を持つため定数 $tau$ は厳密に相殺します。`epa_file` を用意して `tau_mode='epa'` にすると $bold(k)$ 依存の $tau$ が入ります。
