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
- RPA（乱雑位相近似）に基づくスピン感受率 $chi_s$ およびペアリング感受率 $phi$ の計算
- 超伝導状態における動的スピン感受率 $chi_s^"SC"$ の計算
- FLEX（揺らぎ交換近似）による自己エネルギーの計算
- 線形 Eliashberg 方程式による超伝導固有値・ギャップ関数の計算
- 非線形 Eliashberg 方程式による自己無撞着 SC ギャップ関数の計算
- キャリア数・サイクロトロン質量・dHvA 振動数の計算
- 不純物・CPA を含むスペクトル関数の計算

pyrpaはすべての設定を`main.py`の上部にある変数を書き換えることで行います。ライブラリのimport行より上にある変数のみを変更することを想定しています。

= 動作環境

以下のPythonパッケージが必要です。

- `numpy`
- `scipy`
- `matplotlib`
- `skimage`（フェルミ面プロット option=2,3 およびサイクロトロン質量 option=18 で使用）

Fortran で書かれた内部ライブラリ（`libs/flibs`, `libs/plibs`）がコンパイルされている必要があります。線形 / 非線形 Eliashberg・FLEX・SC chi 計算は OpenMP 並列の Fortran ルーチン経由で実行されます。

= 入力ファイルについて

=== ファイル形式の選択（`ftype`）

ハミルトニアンの入力ファイル形式を`ftype`変数で指定します。`fname`にファイル名またはディレクトリ名（拡張子なし）を文字列で代入してください。

#table(
  columns: (auto, 1fr),
  [*ftype*], [*ファイル形式*],
  [0], [`{fname}/` ディレクトリ内に `ham_r.txt`, `irvec.txt`, `ndegen.txt` が存在する形式],
  [1], [`{fname}.input` ファイル（独自形式）],
  [2], [`{fname}_hr.dat` ファイル（Wannier90のデフォルトホッピングファイル）。*第一原理計算との連携では通常これを使用*],
  [3], [MLO（最局在ワニエ）基底の非直交基底ホッピング形式],
  [その他], [`Hopping.dat` ファイル（ecaljのホッピングファイル）],
)

*Wannier90を使用する場合の例：*

```python
fname = 'inputs/Sr2RuO4'   # _hr.dat を除いたパスを指定
ftype = 2
```

=== スピン軌道相互作用（`sw_soc`）

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
  [4], [三方晶格子（Trigonal, QEの sbrav=5 相当）],
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

  - option=19 (`DHVA`): dHvA 振動数 vs 角度
  - option=20 (`ELECTRON_MASS`): symmetry line 上の電子質量
  - option=21 (`SPECTRUM_IMPURITY`): 不純物入りスペクトル関数
  - option=23 (`NONLIN_ELIASHBERG`): 非線形 Eliashberg 方程式
]

=== option=0: バンドプロット (`BAND`)

symmetry line に沿ったバンド分散 $E_n(bold(k))$ を計算し描画します。対称線は `k_sets` と `xlabel` で指定するか, 指定がなければ `brav` の設定から自動で選びます。

=== option=1: 状態密度 (`DOS`)

$ "DOS"(omega) = - 1/pi sum_(bold(k), n) "Im" G^0_n (bold(k), omega + i delta) $

全軌道の合計と軌道ごとの部分 DOS をプロットします。

=== option=2: 2 次元フェルミ面 (`FERMI_2D`)

指定した $k_z$ 断面（デフォルト: $k_z = 0$, 変数 `kz` で変更可能）における $k_x$-$k_y$ 面のフェルミ面を描画します。電子数 `fill` またはフェルミエネルギーから化学ポテンシャルを自己無撞着に求め, $E_n(bold(k)) = mu$ となる等エネルギー線を skimage の find\_contours で抽出します。`RotMat` で回転行列を指定することでフェルミ面を回転できます。

=== option=3: 3 次元フェルミ面 (`FERMI_3D`)

Marching cubes アルゴリズムを用いて 3 次元フェルミ面をポリゴンで描画します。`kscale` で各軸方向の表示スケールを変えることができます。

=== option=4: スペクトル関数 (`SPECTRUM`)

一体グリーン関数のトレースの虚部から電子のスペクトル関数

$ A(bold(k), omega) = -1/pi "Im" "Tr" G(bold(k), omega + i delta) $

を計算し, symmetry line に沿って描画します。`sw_self=True` の場合は FLEX 計算で求めた自己エネルギーを取り込み, 相互作用の効果を含めたスペクトルを計算します。

=== option=5: ボルツマン伝導係数 (`CONDUCTIVITY_BT`)

ボルツマン輸送方程式を用いて以下の各種輸送係数を計算します。緩和時間近似（$tau = $ `tau_const` fs）を使用します。

出力される物理量：
- 電気伝導率テンソル $sigma_{i j}$ （単位: S/m）
- 熱伝導率テンソル $kappa_{i j}$ （単位: W/m/K）
- ゼーベック係数テンソル $S_{i j}$ （単位: V/K）
- パワーファクター $sigma S^2$ （単位: $upright(W \/ m \/K^2)$）
- ペルチェ係数, ローレンツ数

=== option=6: 線形応答理論による伝導係数 (`CONDUCTIVITY_PT`)

久保公式（Kubo formula）に基づく線形応答理論を用いて電気伝導率などを計算します。`delta` パラメータが実効的な緩和時間 $tau approx hbar / delta$ に対応します。

=== option=7: スピン感受率スペクトル (`CHIS_SPECTRUM`)

RPA によるスピン感受率

$ chi_s (bold(q), omega) = frac(chi^0 (bold(q), omega), 1 - S chi^0 (bold(q), omega)) $

を symmetry line に沿って計算し, 虚部のスペクトルをプロットします（パラマグノンスペクトル）。相互作用行列 $S$ はオンサイト相互作用 `U`, `J` から生成されます。結果は `chis_spec.png` に出力されます。

=== option=8: 特定 $bold(q)$ 点でのスピン感受率 (`CHIS_QPOINT`)

`at_point` で指定した波数 $bold(q)$ 点における $chi_s (bold(q), omega)$ の周波数依存性を計算・プロットします。

=== option=9: $bold(q)$ 空間スピン感受率マップ (`CHIS_QMAP`)

`Ecut` で指定したエネルギー $omega_0$ における $chi_s (bold(q), omega_0)$ の $k_x$-$k_y$ 空間分布をプロットします。ネスティングベクトルの可視化などに使用します。出力は `chi0map.png` および `chismap.png`。

=== option=10: ペアリング感受率スペクトル (`PHI_SPECTRUM`)

ペアリング感受率 $phi(bold(q), omega)$ を symmetry line に沿って計算しプロット（`phi_spec.png`）。

=== option=11: $bold(q)$ 空間ペアリング感受率マップ (`PHI_QMAP`)

`Ecut` で指定したエネルギーにおけるペアリング感受率の $bold(q)$ 空間分布をプロットします（`phimap.png`）。`sw_omega` で実周波数 / 松原周波数を切り替えます。

=== option=12: 超伝導状態のスピン感受率スペクトル (`CHIS_SPECTRUM_SC`)

非ゼロのギャップ関数 $Delta(bold(k))$ を仮定し, 異常グリーン関数 $F(bold(k))$ を含めた SC 状態の既約感受率

$ chi_0^"SC" = chi^"GG" plus.minus chi^"FF" $

を構築した上で RPA で動的スピン感受率 $chi_s^"SC" = chi_0^"SC" \/ (1 - S chi_0^"SC")$ を求め, symmetry line 沿いにスペクトルを描画します。`gap_sym` で初期対称性を選びます（負値で triplet）。出力は `chis_sc_spec.png`。

初期ギャップ振幅 `delta0` の指定方法は 2 通りあります（詳細は後述）：
- 単一の実数（例: `delta0=1.e-2`）: 全バンド共通の最大振幅としてスケールした単一バンドギャップを使用。
- バンドごとの値のリスト（例: `delta0=[0.,0.2,0.3,-0.1,0.]`, 長さ `Norb`）: バンド毎に振幅・符号（multi-gap, 符号反転を含む sign-changing multi-band gap）を指定可能。マルチギャップ系（例: 鉄系超伝導体）の符号反転 $s^plus.minus$ 状態などを扱う場合に使用します。

=== option=13: 超伝導状態の特定 $bold(q)$ 点スピン感受率 (`CHIS_QPOINT_SC`)

option=12 と同じ SC 状態の枠組みで, `at_point` で指定した $bold(q)$ 点での $chi_s^"SC"(bold(q), omega)$ の周波数依存性を計算・プロットします。`delta0` の指定方法は option=12 と共通です。

実行時には以下も合わせて出力されます。

- $bold(k) = (0,0,0) arrow.r (0,0.5,0)$ 方向の BdG（Bogoliubov–de Gennes）バンド分散（`BdG_band.png`）
- $chi_s^"SC"(omega)$ の全軌道トレース（`chisq.png`, `chis_sc.dat`）
- 軌道分解した $chi_s^"SC"(omega)$（`chisq_orb.png`, `chis_scorb.dat`）

=== option=14: FLEX 自己エネルギー計算 (`FLEX`)

FLEX アルゴリズムを自己無撞着に解き, 電子の自己エネルギー $Sigma(bold(k), i omega_n)$ を松原周波数空間で計算します。得られた自己エネルギーは `sigma.bin` および `self_en.npz` に保存されます。

- `sw_out_self=True`: 自己エネルギーをファイルへ出力
- `sw_in_self=True`: 前回の自己エネルギー (`sigma.bin`) を初期値として読み込む
- `m_diis_num`: DIIS 履歴長（指定なしのデフォルトは 5）
- `sw_rescale_flex=True`: 自己エネルギーが発散しそうな場合に max$|Sigma| tilde.eq U$ となるよう動的にスケーリング

=== option=15: 線形 Eliashberg 方程式 (`LIN_ELIASHBERG`)

RPA または FLEX で求めた有効ペアリング相互作用を使い, 線形化された Eliashberg 方程式を固有値問題として解きます。最大固有値 $lambda$ と対応するギャップ関数 $Delta(bold(k), i omega_n)$ が得られます。$lambda = 1$ で超伝導転移温度 $T_c$ となります。デフォルトは shift + deflation 付き power method, `arnoldi_m > 0` で Arnoldi 法も選択可能です。

- `sw_self=False`: 自己エネルギー補正なし（RPA）
- `sw_self=True`: FLEX の自己エネルギーを用いる（`sigma.bin` が必要）
- `sw_from_file=True`: `sigma.bin` からファイルを読み込み, FLEX 計算をスキップ
- `gap_sym`: 初期ギャップ関数の対称性を指定（詳細は後述）
- ギャップ関数は `gap_{ij}.dat` および `gap.npy` に出力される

=== option=16: ギャップ関数の後処理 (`GAP_FUNCTION`)

option=15 で計算したギャップ関数 (`gap.npy`) を読み込み, 以下の量を計算・出力します。

- ギャップ関数の異常グリーン関数 $F(bold(k), i omega_n)$
- Padé 近似による実周波数への解析接続スペクトル

=== option=17: キャリア数 (`CARRIER_NUM`)

各バンドのフェルミ面の体積から電子数（粒子数）・ホール数を計算します。`fill` との整合性確認に使います。

=== option=18: サイクロトロン質量 (`CYCLOTRON_MASS`)

オンサーガー関係 $m^*_c = (planck^2 \/ 2 pi)(partial S \/ partial E)$ に基づきサイクロトロン質量を計算します（$S$ はフェルミ面断面積）。

- Phase 1: $k_z in [0, pi/2]$ を `meshkz=20` 点で走査し $S(k_z)$ をマップ
- Phase 2: 極値 $k_z$ で central finite difference により $partial S \/ partial E$ を計算し $m^*_c$ を $m_e$ 単位で出力

skimage の find\_contours と marching\_cubes を内部で使用します。

=== option=19: dHvA 振動数 vs 角度 (`DHVA`)

オンサーガー関係 $F = planck \/ (2 pi e) dot A_"ext"$ を用いて, 磁場の極角 $theta$ に対する dHvA 振動数 $F(theta)$ を計算します。`theta_list` は 0–90° を 40 点で走査します。

=== option=20: 電子質量 (`ELECTRON_MASS`)

symmetry line 上で $m^* = planck^2 (partial^2 E \/ partial bold(k)^2)^{-1}$ を計算します（$m_e$ 単位）。

=== option=21: 不純物入りスペクトル (`SPECTRUM_IMPURITY`)

実空間スーパーセル上で不純物を加えたハミルトニアンを構築し, 不純物グリーン関数を介してスペクトル関数 $A(bold(k), omega)$ を計算します。`imp_list` で不純物サイトを指定。

=== option=22: CPA 伝導 / スペクトル (`SIGMA_CPA`)

Coherent Potential Approximation により混晶の自己エネルギーを反復解で求め, 実周波数のスペクトル関数（`cpa_spectrum.png`）と松原周波数の自己エネルギー（`sigma.bin`, `self_en.npz`）を出力します。`x_cpa`（不純物濃度）, `VA`, `VB`（オンサイト摂動）はサブルーチン内で固定値が設定されています。

=== option=23: 非線形 Eliashberg 方程式 (`NONLIN_ELIASHBERG`)

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

= カラープロット設定（`color_option`）

option=0,2,3では, バンドやフェルミ面の各点に色をつけることができます。

#table(
  columns: (auto, 1fr),
  [*color_option*], [*意味*],
  [0], [色なし（黒単色）],
  [1], [`olist`で指定した軌道成分をRGB（赤/緑/青）の混合比で表示],
  [2], [群速度 $|bold(v)(bold(k))|$ の大きさを色で表示],
)

`color_option=1`の場合, `olist`で軌道インデックスを `[R成分, G成分, B成分]` の形で指定します。複数の軌道を同じ色に割り当てる場合は, リストの要素をさらにリストにします。例:

```python
olist = [[0, 3], [1, 4], [2, 5]]
# 軌道0と3をRed, 軌道1と4をGreen, 軌道2と5をBlueに対応
```

= 各種パラメータの詳細説明

以下では, `main.py`の上部に記載されている全パラメータについて物理的な意味を含めて説明します。

=== 基本設定

- `fname` （文字列）: 入力ファイルのパス（`ftype`に応じた形式で指定, 拡張子は通常不要）

- `ftype` （整数）: 入力ハミルトニアンのファイル形式（3節参照）

- `brav` （整数）: ブラベー格子の種類（4節参照）

- `sw_soc` （bool）: スピン軌道相互作用を考慮するかどうかのスイッチ

=== メッシュ設定

- `Nx, Ny, Nz` （整数）: 第一ブリルアンゾーン内の $bold(k)$ 点メッシュ数。フェルミ面積分や FFT-convolution（FLEX, Eliashberg, chi）で使用します。2 次元系では `Nz=1` とします。FLEX/Eliashberg 計算では FFT 効率の都合上 2 のべき乗（32, 64 など）が望ましく, メモリ消費は $tilde.equiv N_x N_y N_z dot N_w dot N_"orb"^2$ で増えます。

- `Nw` （整数）: 松原周波数のメッシュ数（option=12,13,14,15,16,23 の FLEX/Eliashberg/SC chi 計算）または実周波数のメッシュ数（option=1,4 など）。温度 $T$ に対して $omega_n = (2n+1) pi T$ の形で $n=0,1,...,N_w-1$ まで取ります。低温ほど大きな $N_w$ が必要です（典型値: 256–1024）。

- `kmesh` （整数）: バンド・スペクトルをsymmetry lineに沿って計算する際の$bold(k)$点数。大きいほど滑らかなプロットが得られます（例: 200～500程度）。

=== 格子定数

- `abc` （リスト, 単位: Å）: 格子定数$a, b, c$。速度 $v = 1/hbar (partial E)/( partial bold(k))$ の計算やsymmetry lineの幅の計算（実空間の長さ換算）に使用します。Wannier90の場合はWIN/WOUTファイルの値に合わせてください。

- `alpha_beta_gamma` （リスト, 単位: degree）: 格子角度$alpha, beta, gamma$。直方晶・立方晶の場合は `[90., 90., 90.]` です。

=== 温度・化学ポテンシャル

- `tempK` （実数, 単位: K）: 温度をケルビンで指定します。コード内で $T = k_B$`tempK`（単位: eV）に変換されます。

- `temp` （実数, 単位: eV）: $k_B T$を直接eV単位で指定する場合に使います。`tempK`と同時に定義した場合は`temp`が優先されます。

- `fill` （実数）: バンドフィリング。化学ポテンシャル $mu$ を $sum_(bold(k), n) f(epsilon_n (bold(k)) - mu) = N_k dot f"ill"$ で求めるための目標占有数で, 0 から軌道数 `no` までの値を取ります。`fill = no` で full-filled。物理的な 1 ユニットセルあたりの電子数との対応はハミルトニアンの形式に依存します:
  - SOC なし (`sw_soc=False`): スピン自由度は陽に含まれないので, `fill` はスピンあたりの電子数に対応。たとえば 3 軌道モデルで `fill=3` は full-filled（全スピン込みで 6 電子相当）, `fill=1.5` が half-filled。
  - SOC あり (`sw_soc=True`): ハミルトニアンが既にスピンを含む（`no = 2 N_"orb"`）ので, `fill` は 1 セルあたりの全電子数そのもの。3 軌道（`no=6`）で `fill=6` が full-filled, `fill=3` が half-filled。

- `mu0` （実数, 単位: eV）: 化学ポテンシャルを直接指定したい場合に定義します。定義されている場合, `fill`からの自己無撞着計算をスキップして`mu0`の値を使用します。

=== エネルギー範囲・収束パラメータ

- `Emin, Emax` （実数, 単位: eV）: DOSやスペクトルを計算する際のエネルギー範囲の下限・上限。

- `delta` （実数, 単位: eV）: スペクトル計算のブロードニングパラメータ $delta$。グリーン関数の虚部に加える小量で, ディラックのデルタ関数を有限幅のローレンツ型に広げます。大きすぎると特徴が潰れ, 小さすぎると数値ノイズが出ます（0.01〜0.05 eV程度が目安）。

- `Ecut` （実数, 単位: eV）: option=9,11 で $bold(q)$ 空間の感受率マップを作成する際の固定エネルギー $omega_0$（フェルミ面近傍は 0 付近に設定）。

- `delta0` （実数 または 長さ`Norb`のリスト, 単位: eV）: option=12,13 の SC chi 計算で用いる初期ギャップ関数の振幅。物理的には超伝導ギャップサイズに対応します（典型値: $10^{-3} - 10^{-2}$ eV ≒ 1–10 meV）。
  - 実数を指定した場合: 全バンド共通の単一ギャップとして, 内部で生成した対称性形状を最大振幅 `delta0` にスケールします。0 にすると常状態と等価になります。
  - リストを指定した場合（例: `delta0=[0.,0.2,0.3,-0.1,0.]`）: バンド毎に振幅と符号を個別に指定する multi-gap モードになります。符号が異なる成分を混在させることで $s^plus.minus$ 型など符号反転ギャップを表現できます。

=== 伝導計算パラメータ

- `tau_const` （実数, 単位: fs）: ボルツマン理論（option=5）での定数緩和時間$tau$。実際の系では散乱機構に依存しますが, 定数近似では実験値との比較から決めます（例: 金属では1〜100 fs程度）。

- `sw_tdf` （bool）: `True`にすると輸送分布関数（TDF）を先に計算し, そのエネルギー積分から輸送係数を求めます（エネルギー依存緩和時間に対応する場合に有効）。

=== 軌道・相互作用パラメータ

- `olist` （リスト）: カラープロット（color_option=1）で使用する軌道インデックスのリスト（6節参照）。

- `U` （実数, 単位: eV）: オンサイトクーロン斥力（ハバード$U$）。FLEX/RPA計算で使用。磁性・超伝導の計算には重要なパラメータです。

- `J` （実数, 単位: eV）: オンサイトフント結合定数。$U' = U - 2J$（金森遮蔽）が自動で計算されます。

- `orb_dep` （bool）: `True`にすると軌道依存の相互作用行列 `Umat`, `Jmat` を使用します。`False`（デフォルト）では `U`, `J` の一定値を全軌道に適用します。

- `m_diis_num` （整数, 任意）: FLEX (option=14) および非線形 Eliashberg (option=23) で使う DIIS（Pulay 加速混合）の履歴長。$2$ 以上で Pulay 外挿が有効, $1$ で線形混合へフォールバック。未定義の場合は 5 になります。値を大きくするほど収束は速まりますがメモリ消費が増えます（$N_x N_y N_z dot N_w dot N_"orb"^2$ 倍）。

=== ギャップ関数の対称性（`gap_sym`）

Eliashberg方程式（option=15,23）を解く際, または SC chi 計算（option=12,13）の初期ギャップ形状を生成する際の対称性を指定します。

#table(
  columns: (auto, 1fr),
  [*gap_sym*], [*対称性*],
  [0], [$s$波（全ての$bold(k)$で正）],
  [1], [$d_{x^2-y^2}$波（$cos k_x - cos k_y$型）],
  [2], [$s^plus.minus$波（ネスティングベクトルを挟んで符号反転）],
  [3], [$d_{x y}$波（$sin k_x sin k_y$型）],
  [-1], [$p_x$波],
  [-2], [$p_y$波],
)

=== $bold(k)$空間の設定

- `kz` （実数, 規格化単位）: option=2のフェルミ面プロットにおける$k_z$の値。$0 \sim 0.5$の範囲で指定します（$k_z = 0$は$Gamma$点面, $k_z = 0.5$はゾーン境界面）。

- `kscale` （リスト or 実数）: option=3の3次元フェルミ面表示における各軸の表示スケール。例: `kscale=[1.0, 1.0, 0.5]`とすると$k_z$方向を半分に縮めて表示。

- `k_sets` （リストのリスト）: symmetry lineの各端点のk座標（規格化単位, $0 \sim 1$）。自動生成ではなく手動で対称線を設定したい場合に定義します。例: `k_sets=[[0,0,0],[0.5,0,0],[0.5,0.5,0]]`

- `xlabel` （文字列のリスト）: `k_sets`に対応するラベル。例: `xlabel=[r'$\Gamma$','X','M']`

- `at_point` （リスト）: option=8で感受率を計算する$bold(q)$点の座標（規格化単位）。

=== スイッチ類

- `sw_unit` （bool）: `True`（デフォルト）で国際単位系（SI）の物理定数を使用し, 出力も実際の物理単位で得られます。`False`にすると$hbar = k_upright(B) = e = 1$の無次元系になります。

- `sw_omega` （bool）: option=11 でペアリング感受率を実周波数（`True`）または松原周波数（`False`）で計算するかのスイッチ。

- `sw_self` （bool）: option=4 (スペクトル) / option=15 (線形 Eliashberg) / option=23 (非線形 Eliashberg) で FLEX の自己エネルギーを使うかどうか。`True` にすると繰り込まれたグリーン関数を用います。

- `sw_out_self` （bool）: `True` にすると FLEX 計算の自己エネルギーを `sigma.bin` および `self_en.npz` に出力します。

- `sw_in_self` （bool）: `True` にすると FLEX 計算の初期自己エネルギーとして前回の `sigma.bin` を読み込みます（自己無撞着ループの収束を早める）。

- `sw_from_file` （bool）: `True` にすると自己エネルギーを `sigma.bin` から読み込み, FLEX 計算を実行せずに直接 Eliashberg 計算を行います。

- `sw_check_only` （bool）: option=23（非線形 Eliashberg）専用。`True` にすると内部の線形化 Eliashberg 解（Stoner 因子 $S$ と固有値 $lambda_"eliash"$）を計算した時点で停止し, 非線形ループを実行しません。温度走査で $T_c$（$lambda_"eliash"=1$ となる温度）の見積もりだけを高速に行いたい場合に使用します。$S >= 1$ または $lambda_"eliash" < 1$ の場合は, このフラグの値に関わらず非線形ループ手前で自動的に停止します。

- `sw_rescale_flex` （bool）: option=14 の FLEX 計算で, 反復中に max$|Sigma|$ が `U` を超えそうな場合に動的にスケーリングを掛けて発散を防ぎます。Stoner 因子が 1 に近い系で有効。

- `sw_dec_axis` （bool）: True にすると格子ベクトルを適切に分解して逆格子ベクトルを設定します。

= 典型的な計算フローの例

=== フェルミ面とバンドの確認

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

=== 超伝導ギャップ関数の計算（RPA 法）

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

=== FLEX + 線形 Eliashberg

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

=== 非線形 Eliashberg 方程式（自己無撞着 SC ループ）

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

=== 超伝導状態の動的スピン感受率

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

`tests/` ディレクトリには, 主要な数値カーネルと物理ベンチマークを確認するための回帰テストが含まれています。`pytest` がインストールされていない環境でも, 各テストファイルを Python スクリプトとして直接実行できます。

実行前に Fortran 共有ライブラリ `libs/libfmod.so` がコンパイル済みである必要があります。未コンパイルの場合は `libs` ディレクトリで `make FC=<compiler> SL=<library>` を実行してください。

=== 実行方法

個別に実行する場合:

```bash
python tests/test_eilenberger.py
python tests/test_rpa_flex.py
```

`pytest` が利用できる環境では以下でも実行できます。

```bash
pytest tests
```

=== `tests/test_eilenberger.py`

準古典 Eilenberger / Riccati ソルバーのテストです。均一系, 表面, 渦糸, 渦糸格子, model Fermi surface, Pauli 制限, triplet $d$-vector texture など, Eilenberger 関連の物理ベンチマークを確認します。

主な確認項目:

- 松原周波数 cutoff の温度スケーリング
- Anderson theorem と Abrikosov--Gor'kov pair breaking の再現
- 弱結合 BCS 比 $2 Delta_0 / k_B T_c approx 3.53$
- Fortran Riccati カーネルと Python 参照実装の一致
  - scalar `riccati_chords`
  - spin-matrix `matrix_riccati_batch`
  - batched chord `matrix_riccati_chords`
- $d$波表面でのギャップ抑制と zero-energy bound state
- 渦糸芯での CdGM zero-energy peak
- 渦糸格子での Volovik 型 field dependence
- `build_model_fs` の規格化と等方 Fermi surface 極限
- singlet Pauli 抑制と triplet equal-spin pairing の耐性
- 表面・渦糸芯における triplet $d$-vector texture

このテストは比較的物理寄りのベンチマークであり, 実行時間は環境により数十秒程度かかる場合があります。

=== `tests/test_rpa_flex.py`

RPA / FLEX / Eliashberg の基本部品を対象にした軽量な回帰テストです。Eilenberger 以外の Fortran wrapper, RPA 行列演算, FLEX 用 bubble / vertex, Eliashberg solver の smoke test を確認します。

主な確認項目:

- `get_chi_orb_list` の multi-site orbital-pair / site index の整合性
- `gen_SCmatrix` の 2 軌道 Kanamori 型 vertex の基準値
- `gen_SCmatrix_orb` の軌道依存 `Umat`, `Jmat` vertex の基準値
- 1 軌道 RPA 公式
  $ chi_s = chi^0 / (1 - U chi^0) $
  の確認
- `S=C=0` のとき `get_chis_chic` が bare $chi^0$ に戻ること
- 1 軌道 tight-binding model に対する `gen_Green0` の解析式
  $ G^0(k,i omega_n) = 1 / (i omega_n + mu - epsilon_k) $
  の確認
- `get_chi0` と `get_chi0_conv` の一致
- `get_Vsigma_nosoc_flex` が相互作用ゼロでゼロ vertex を返すこと
- `linearized_eliashberg` が相互作用ゼロで $lambda = 0$ を返し, 有限な配列を返すこと
- `nonlinear_eliashberg` がゼロ seed / ゼロ相互作用で自明解 $Delta=0$ を保つこと
- `_load_sigma_from_file` が `self_en.npz` 不在時にクラッシュせず `None` を返すこと
- `output_gap_function` の 1 軌道出力 smoke test

このテストは小さな 1 軌道模型と小さい $k$ メッシュを使うため, RPA/FLEX 周辺の API 破壊を素早く検出する用途に向いています。

= トラブルシューティング

- *化学ポテンシャルが収束しない*: `Nx, Ny, Nz` を増やすか, `tempK` を少し上げてください。

- *スペクトルにノイズが多い*: `delta` を少し大きくする, または `Nw` を増やしてください。

- *バンドがおかしい*: `brav` の設定が合っているか確認してください。特に QE との連携では `brav=1`（FCC）や `brav=2`（BCC）の convention に注意が必要です。

- *線形 Eliashberg が収束しない*: `U` を小さくする, `tempK` を高くする, メッシュ `Nx,Ny,Nz` を増やすなどを試してください。固有値 $lambda$ が 1 を超えている温度は超伝導相内に対応します。

- *非線形 Eliashberg が発散する*: `m_diis_num` を 5–10 に増やす, `tempK` を上げる, `Nw` を増やすなどを試してください。`sw_check_only=True` で先に Stoner 因子 $S$ と固有値 $lambda_"eliash"$ を確認するのも有効です。$S$ が 1 を超えている場合は磁気不安定なので `U` を下げる必要があり, $lambda_"eliash" < 1$ の場合はまだ $T_c$ 以上なので降温してください（いずれの場合も非線形ループは自動的に実行されません）。

- *FLEX で max$|Sigma|$ が発散*: `sw_rescale_flex=True` を設定するか, `U` を下げてください。

- *伝導率の単位が合わない*: `sw_unit=True` になっているか確認してください。`False` にすると無次元系になります。
