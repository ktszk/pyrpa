# pyrpa test-suite

前提: Fortranライブラリのビルド (`cd libs && make FC=ifx SL=MKL`)。

## 実行方法

```sh
# 推奨: pytest (リポジトリルートから; pytest.ini が tests/ を拾う)
pytest                        # 全テスト
pytest -k chi0                # 名前でフィルタ
pytest --runslow              # @pytest.mark.slow 付きも実行

# pytest が無い環境 (クラスタ等): 各ファイルは単体実行できる
python tests/test_rpa_flex.py           # 1ファイル
python tests/test_rpa_flex.py chi0      # 名前に 'chi0' を含むテストのみ
python tests/run_all.py                 # 全ファイル一括 + サマリ表
python tests/run_all.py rpa trans -q    # ファイル名フィルタ + 静音
```

## ファイル構成

| ファイル | 対象 |
|---|---|
| `test_rpa_flex.py` | RPA / FLEX / Eliashberg (χ0, 頂点, tail補正, mkself) |
| `test_eilenberger.py` | 準古典 Eilenberger (一様系・表面・渦・渦格子) |
| `test_transport.py` | 輸送係数 (Boltzmann / Kubo, Wiedemann-Franz) |
| `test_tools.py` | 共有ツール `_tools.py` 自体の自己検証 (使用例を兼ねる) |
| `_tools.py` | 共有ユーティリティ (下記) |
| `conftest.py` | pytest 設定 (sys.path, `slow` マーカー, `silence` fixture) |
| `run_all.py` | pytest 不要の一括ランナー |

## 共有ツール `_tools.py`

テストファイル(tests/ 直下)から `import _tools as T` で利用する。

- **モデルビルダー** — Fortran 入口に渡す一式 (klist/kmap/invk/hamk/eig/uni/Gk/olist/Smat/Cmat...) を dict で返す
  - `one_orbital_square(Nx, Ny, Nz, Nw, temp, mu, tx, ty, U, J)` — 1軌道正方格子。`eps_full` (フルグリッド分散) 付きで厳密参照と直結
  - `two_orbital_square(..., t1, t2, de, V)` — 2軌道 + 混成。バンドは閉形式 `(e1+e2)/2 ± sqrt(((e1-e2)/2)² + V²)` で検算可能
  - `band_diagonal_model(Nk, Norb, seed)` — 輸送テスト用ランダム帯対角モデル
  - `cylindrical_fs(Nb, gap_sym)` — Eilenberger 用円筒 FS (`<φ²>=1` 規格化)
- **厳密参照解** — Matsubara カットオフ誤差を含まない解析解
  - `exact_chi0_1orb(eps_full, qlist, mlist, temp, mu, a, delta)` — Lindhard 泡 (コードの符号規約 χ>0)。`a≠0` で厳密可解な2極自己エネルギー Σ = a²/(iω+μ−δ) 入りにも対応 → χ0 の Nw 収束次数テストの基準
  - `two_pole_green(eig, temp, mu, Nw, a, delta)` — 対応する G (コード配置 [1,1,Nw,Nk])
- **数値ヘルパー**
  - `fit_power_order(ns, errs)` — err ∼ C/n^p の p を対数フィット (sharp cutoff ≈1, tail補正 ≈2 の固定に使用)
  - `rel_err(a, b)`
- **物理アサーション** — 「物理として成立すべき性質」を一行で検査
  - `assert_hermitian(m)` / `assert_positive_semidefinite(m)` — 静的 χs/χc 行列など
  - `assert_green_causal(Gk)` — Im G_aa(iω_n>0) ≤ 0 (因果律)
  - `assert_raises(exc, fn, ...)` — pytest 無しでも動く pytest.raises 代替
- **その他**
  - `fermi / bose / fermionic_matsubara / bosonic_matsubara` — オーバーフロー安全な分布関数
  - `silence()` — Fortran ラッパの進捗出力を抑制するコンテキストマネージャ
  - `run_standalone(globals())` — 各ファイル末尾の `__main__` ランナー共通化 (tmp_path 対応, 引数で名前フィルタ)

## 新しいテストを書くときの型

```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _tools as T
import libs.flibs as F

def test_something():
    md = T.one_orbital_square(Nx=8, Ny=8, Nw=64, temp=0.05, mu=0.3)
    with T.silence():
        chi, stoner = F.get_chi0(md['Smat'], md['Cmat'], md['Gk'], md['olist'],
                                 md['kmap'], md['invk'], md['temp'],
                                 md['Nx'], md['Ny'], md['Nz'])
    ref = T.exact_chi0_1orb(md['eps_full'], md['klist'], [0], md['temp'], md['mu'])
    assert np.abs(chi[0, 0, 0, :] - ref[:, 0]).max() < ...

if __name__ == '__main__':
    sys.exit(T.run_standalone(globals()))
```

時間のかかるテストは `@pytest.mark.slow` を付ける (pytest では `--runslow` 指定時のみ実行; 単体実行では常に走る)。
