//文章フォーマット設定
#set text(lang: "ja")
#set text( font: ("Times New Roman","Yu Mincho"), size: 11pt) //全体のフォントサイズと英語フォントの設定
#let heading_font(body)={
  set text(font: ("Arial", "Yu Gothic"))
  body
}
#set page(paper: "a4",margin: (x:2cm, y:2cm),)
#set math.equation(numbering: numbering.with("(1)"), supplement: "式")
//外部パッケージ
#import "@preview/physica:0.9.3": *

= ハミルトニアンの対称性とグリーン関数への影響
しょっちゅう使う割には毎度導出しているので, すぐ確認するように纏めた。基本的には周期性のある場合を考えています。
また軌道は実関数軌道を想定しています。複素軌道を用いる場合はスピン同様時間反転による符号の変化を考慮してください。
== エルミート性のみを考慮した時に現れる対称性
=== ハミルトニアン
ハミルトニアンは(最近は非エルミート系の議論もあるが, )通常エルミートであるので次の対称性を持つ
$ hat(H)^(dagger)=hat(H) $
波数で対角的なとき, これは各波数におけるハミルトニアン$H(k)$でも成り立つ為, 
各波数におけるWannier関数や軌道を基底とするハミルトニアンは次の対称性を持つ
$ hat(H)(k)=hat(H)^(dagger)(k) $
$ H^(l m)(k)=[H^(m l)(k)]^* $
強束縛模型で考えると
$ hat(H)(k)=sum_r hat(t)(r)exp(i k dot r) $
なのでホッピングは自動的に
$ [hat(t)(r)]^dagger=hat(t)(-r) $
を満たすこととなる。
=== 松原グリーン関数
相互作用の無い時の松原グリーン関数の定義は
$ hat(G)_0(k,i epsilon_m)&=(i epsilon_m hat(I) - hat(H)(k))^(-1)\
&=hat(U)(k)[i epsilon_m hat(I)-hat(E)(k)]^(-1)hat(U)^(dagger)(k) $
$ G_0^(l l')(k,i epsilon_m)=sum_n (U^(l n)(k)(U^(dagger))^(n l')(k))/(i epsilon_m -xi_n (k)) $
なので, 上記のハミルトニアンのエルミート性から以下の対称性が導かれる。
$ hat(G)_0^(dagger)(k,i epsilon_m)&=[(i epsilon_m hat(I) - hat(H)(k))^(-1)]^(dagger)\
&=(-i epsilon_m hat(I) - hat(H)^(dagger)(k))^(-1)\
&=hat(G)_0(k,-i epsilon_m) $
$ [G_0^(l l')(k,i epsilon_m)]^(*)= G_0^(l' l)(k,-i epsilon_m) $
通常自己エネルギーも同様の対称性を持つので, 相互作用による自己エネルギーを考慮した$G(k,i epsilon_n)$も同様の対称性を持つ
=== 既約感受率
松原グリーン関数の性質を用いると, 既約感受率はその定義
$ chi_(0)^(l_1,l_2,l_3,l_4)(q,i omega_n)=sum_(k,i epsilon_m) G^(l_1,l_3)(k+q,i(epsilon_m+omega_n))G^(l_4,l_2)(k,i epsilon_m) $
から, 次の対称性を持つことが分かる。
$ [chi_(0)^(l_1,l_2,l_3,l_4)(q,i omega_n)]^* &=[sum_(k,i epsilon_m) G^(l_1,l_3)(k+q,i(epsilon_m+omega_n))G^(l_4,l_2)(k,i epsilon_m)]^*\
&=sum_(k,i epsilon_m) [G^(l_1,l_3)(k+q,i(epsilon_m+omega_n))]^* [G^(l_4,l_2)(k,i epsilon_m)]^*\
&=sum_(k,-i epsilon_m) G^(l_3,l_1)(k+q,-i(epsilon_m+omega_n))G^(l_2,l_4)(k,-i epsilon_m)\
&=chi_0^(l_3,l_4,l_1,l_2)(q,-i omega_n) $
$l_1,l_2$及び$l_3,l_4$の軌道の足をペアに取り, 4階のテンソルから2階のテンソルへと射影してやると, 最終的に
$ [hat(chi)_0(q,i omega)]^dagger=hat(chi)_0(q,-i omega) $
のように, グリーン関数と同様の対称性が松原周波数方向で成り立つことが分かる。

同様の基底で相互作用バーテックス$hat(S),hat(C)$を表すと,
スピン/電荷感受率
$ hat(chi)_s&=(hat(I)-hat(chi)_0hat(S))^(-1)hat(chi)_0=hat(chi)_0(hat(I)-hat(S)hat(chi)_0)^(-1)\
  hat(chi)_c&=(hat(I)+hat(chi)_0hat(C))^(-1)hat(chi)_0=hat(chi)_0(hat(I)+hat(C)hat(chi)_0)^(-1) $
は既約感受率と同様の対称性を持つことになるので
$ [hat(chi)_s (q,i omega_n)]^dagger&=hat(chi)_s (q,-i omega_n)\
  [hat(chi)_c (q,i omega_n)]^dagger&=hat(chi)_c (q,-i omega_n) $
を満たす。

また, 既約関数は$k'=k+q,epsilon'_m=epsilon_m+omega_n$と置き,式変形を行うと 
$ chi_(0)^(l_1,l_2,l_3,l_4)(q,i omega_n)&=sum_(k,i epsilon_m) G^(l_1,l_3)(k+q,i(epsilon_m+omega_n))G^(l_4,l_2)(k,i epsilon_m)\ 
&=sum_(k,i epsilon_m) G^(l_1,l_3)(k',i epsilon'_m)G^(l_4,l_2)(k'-q,i(epsilon'_m-omega_n))\
&=sum_(k',i epsilon'_m)G^(l_4,l_2)(k'-q,i(epsilon'_m-omega_n))G^(l_1,l_3)(k',i epsilon'_m)\
&=chi_(0)^(l_4,l_3,l_2,l_1)(-q,-i omega_n) $
も成り立つ。
以上より既約感受率が満たす対称性は
$ chi_(0)^(l_1,l_2,l_3,l_4)(q,i omega_n)=chi_(0)^(l_4,l_3,l_2,l_1)(-q,-i omega_n)
=(chi^*_(0))^(l_3,l_4,l_1,l_2)(q,-i omega_n)=(chi^*_(0))^(l_2,l_1,l_4,l_3)(-q,i omega_n) $
これはスピン/電荷感受率も満たす。この対称性は二体の電子間相互作用が満たすべき対称性なので, 感受率は電子間の相互作用として用いる事が出来る。

この得られたグリーン関数と二体の相互作用の対称性を用いると, 自己エネルギーもグリーン関数と同様の対称性を持つことが分かる。

$ [hat(Sigma)(k,i epsilon_m)]^dagger=hat(Sigma)(k,-i epsilon_m) $
=== ギャップ関数・異常グリーン関数
松原方向の議論を行う前にBCSにおけるギャップ関数の対称性を考える。
一般的に超伝導ギャップは$bold(d)$ベクトルを用いて
$ Delta^(l m)(k)=(d^(l m)_0 (k) sigma_0+bold(d)^(l m)(k)dot bold(sigma))i sigma_y $
と表せる。これをスピンに関する行列で書き直すと,
$ Delta^(l m)(k)=mat(-d^(l m)_x (k)+i d^(l m)_y (k),d^(l m)_0 (k)+d^(l m)_z (k);-d^(l m)_0 (k)+d^(l m)_z (k),d^(l m)_x (k)+i d^(l m)_y (k)) $
ギャップ関数の大元が$Delta^(l m)_(sigma,sigma')(k)=expval(c^l_(k,sigma)c^m_(-k,sigma'))$であることを思い出すと, 
$ Delta^(l m)_(sigma,sigma')(k)&=expval(c^l_(k,sigma)c^m_(-k,sigma'))=-expval(c^m_(-k,sigma')c^l_(k,sigma))=-Delta^(m l)_(sigma',sigma)(-k) $
よって$bold(d)$ベクトルは
$ d_0(k)=d^T_0(-k),\
bold(d)(k)=-bold(d)^T (-k) $
を満たす。

異常グリーン関数と超伝導ギャップ(異常自己エネルギー)はBCSハミルトニアンの非対角項に対応するので, エルミート性から自身のみの対称性を議論するのは困難, 
ギャップのスピン,パリティ, 松原の対称性によって, これらの対称性は決まる。
一方共役函数同士に関してはある程度の対称性がある。
異常グリーン関数(とその複素共役)は定義より
$ F_(alpha,beta)(k,tau)=-expval(T_tau [c_(k,alpha)(tau),c_(-k,beta)(0)])\
F^dagger_(alpha,beta)(k,tau)=-expval(T_tau [c^dagger_(-k,alpha)(tau),c^dagger_(k,beta)(0)]) $
この形から分かるように異常グリーン関数は演算子の入れ替えに対する反対称性をもつそのため波数, 松原周波数表記の場合でも,
$ F_(alpha,beta)(k,i epsilon_n)=-F_(beta,alpha)(-k,-i epsilon_n)\
F^dagger_(alpha,beta)(k,i epsilon_n)=-F^dagger_(beta,alpha)(-k,-i epsilon_n) $<psym>
シングレット, トリプレットに対応する異常グリーン関数を考えると, 転置時のスピンの入れ替わりの影響によって
$ F_s (k,i epsilon_n)=F^T_s (-k,-i epsilon_n)\
F_t (k,i epsilon_n)=-F^T_t (-k,-i epsilon_n) $
と言う対称性を最低限持つことが分かる。ギャップ関数が異常自己エネルギーであることから自明ではあるが, これは上の$bold(d)$の対称性を松原周波数にまで拡張した物に対応する。

== 時間反転対称性から来る対称性
=== ハミルトニアン
時間反転演算子$Theta=-i hat(sigma)_y hat(K)$を作用させると, $bold(r)arrow bold(r), bold(p)arrow -bold(p) (bold(k)arrow -bold(k)), sigma arrow macron(sigma)(-sigma)$と変換される。
ここで$hat(K)$は複素共役演算子で作用させると複素共役を返す。スピンを考えない場合$-i hat(sigma)_y$は無視できるため, 時間反転操作は複素共役を求めることに対応する。
基底が軌道(実球面調和関数)やスピンを含まないWannier関数なら時間反転による基底由来の符号の変化は起こらない(複素球面調和関数などの場合, 非対角項の符号が変わる可能性有)。
そのため先のハミルトニアンが時間反転対称性を持つとき以下の対称性を持つ。
$ hat(H)_(sigma,sigma')(k)=sigma sigma'hat(H)^*_(macron(sigma),macron(sigma)')(-k) $
これとエルミート性から
$ hat(H)_(sigma,sigma')(k)=sigma sigma'hat(H)^T_(macron(sigma)',macron(sigma))(-k) $
も成り立つ。$sigma,sigma'$はスピンの符号に対応する。
またこの時固有値を与えるユニタリー行列も
$ hat(U)(k)=-i sigma_y K hat(U)(-k) $
となるので, その成分は$l,sigma,n$を軌道,スピン,バンドの足とすると,
$ U^(l,sigma)_(n)(k)=-sigma [U^(l,macron(sigma))_(n)(-k)]^* $
これはスピン対称性が丸いとき(スピンの依存性が無視できる時)にハミルトニアンは
$ hat(H)(k)=hat(H)^dagger (k)=hat(H)^* (-k)=hat(H)^T (-k) $
となることを示している。

時間反転対称性があるとき, ホッピングは次の関係を満たす。
$ hat(t)_(sigma,sigma')(r)=sigma sigma'hat(t)^*_(macron(sigma),macron(sigma)')(r) $
そのため, スピンが丸い系では$hat(t)(r)=hat(t)^* (r)$を満たすので, ホッピングは実数になる。
=== 松原グリーン関数
先のハミルトニアンの関係を用いると, 松原グリーン関数は次の対称性を持つ。
$ G^(l,l')_(sigma,sigma')(k,i epsilon_m)&=sum_n (U^(l,sigma)_(n)(k)[U^(l',sigma')_(n)(k)]^*)/(i epsilon_m -xi_n (k))\
&=sum_n (-sigma [U^(l,macron(sigma))_(n)(-k)]^* (-sigma')U^(l',macron(sigma)')_(n)(-k))/(i epsilon_m -xi_n (-k))\
&=sigma sigma'[sum_n ( U^(l,macron(sigma))_(n)(-k)[U^(l',macron(sigma)')_(n)(-k)]^*)/(-i epsilon_m -xi_n (-k))]^*\
&=sigma sigma'[G^(l,l')_(macron(sigma),macron(sigma)')(-k,-i epsilon_m)]^*\
&=sigma sigma'G^(l',l)_(macron(sigma)',macron(sigma))(-k,i epsilon_m) $ 
よってスピンが丸いときには
$ hat(G)(k,i epsilon_m)=hat(G)^dagger (k,-i epsilon_m)=hat(G)^* (-k,-i epsilon_m)=hat(G)^T (-k,i epsilon_m) $
という対称性を松原グリーン関数は持つ。
=== 既約感受率
既約感受率に関しても同じ議論を行うと, 
$ (chi_(0))^(l_1,l_2,l_3,l_4)_(sigma_1,sigma_2,sigma_3,sigma_4)(q,i omega_n)&=sum_(k,i epsilon_m) G^(l_1,l_3)_(sigma_1,sigma_3)(k+q,i(epsilon_m+omega_n))G^(l_4,l_2)_(sigma_4,sigma_2)(k,i epsilon_m)\
&=sigma_1 sigma_2 sigma_3 sigma_4sum_(k,i epsilon_m) [G^(l_1,l_3)_(macron(sigma)_1,macron(sigma)_3)(-k-q,-i(epsilon_m+omega_n))]^* [G^(l_4,l_2)_(macron(sigma)_4,macron(sigma)_2)(-k,-i epsilon_m)]^* \
&=sigma_1 sigma_2 sigma_3 sigma_4 [(chi_(0))^(l_1,l_2,l_3,l_4)_(macron(sigma)_1,macron(sigma)_2,macron(sigma)_3,macron(sigma)_4)(-q,-i omega_n)]^* \
&=sigma_1 sigma_2 sigma_3 sigma_4 (chi_(0))^(l_3,l_4,l_1,l_2)_(macron(sigma)_3,macron(sigma)_4,macron(sigma)_1,macron(sigma)_2)(-q,i omega_n) \
&=sigma_1 sigma_2 sigma_3 sigma_4 [(chi_(0))^(l_4,l_3,l_2,l_1)_(macron(sigma)_4,macron(sigma)_3,macron(sigma)_2,macron(sigma)_1)(q,i omega_n)]^* \
&=sigma_1 sigma_2 sigma_3 sigma_4 (chi_(0))^(l_2,l_1,l_4,l_3)_(macron(sigma)_2,macron(sigma)_1,macron(sigma)_4,macron(sigma)_3)(q,-i omega_n) $
となる。スピンが丸いときにはグリーン関数同様
$ hat(chi)_0 (q,i omega_n)=hat(chi)^dagger_0 (q,-i omega_n)=hat(chi)^*_0 (-q,-i omega_n)=hat(chi)^T_0 (-q,i omega_n) $
を満たす。同様の議論は自己エネルギーやスピン/電荷感受率でも行える。
より正確に対称性を考えると, スピンが丸いときの既約感受率は
$ chi_(0)^(l_1,l_2,l_3,l_4)(q,i omega_n)&=(chi^*_(0))^(l_3,l_4,l_1,l_2)(q,-i omega_n)=(chi^*_(0))^(l_1,l_2,l_3,l_4)(-q,-i omega_n)=chi_(0)^(l_3,l_4,l_1,l_2)(-q,i omega_n)\
&=chi_(0)^(l_4,l_3,l_2,l_1)(-q,-i omega_n)=(chi^*_(0))^(l_2,l_1,l_4,l_3)(-q,i omega_n)\
&=(chi^*_(0))^(l_4,l_3,l_2,l_1)(q,i omega_n)=chi_(0)^(l_2,l_1,l_4,l_3)(q,-i omega_n) $

時間反転対称性は超伝導ギャップや異常グリーン関数の対称性にも影響を与える。
通常の超伝導の議論を行う際は時間反転対称性を考慮する。(時間反転対称性を破った議論はFFLOの時などは必要)

=== ギャップ関数・異常グリーン関数
BCSギャップ関数
$ Delta^(l m)(k)=(d^(l m)_0 (k) sigma_0+bold(d)^(l m)(k)dot bold(sigma))i sigma_y $
に時間反転演算子を作用させると, 
$ hat(Theta)^(-1)Delta^(l m)(k)hat(Theta)&=hat(K)i hat(sigma)_y (d^(l m)_0 (k) 
sigma_0+bold(d)^(l m)(k)dot bold(sigma))i sigma_y (-i hat(sigma)_y) hat(K)\
&=i hat(sigma)_y hat(K)(d^(l m)_0 (k) sigma_0+bold(d)^(l m)(k)dot bold(sigma))hat(K) $
これを行列で表すと
$ hat(Theta)^(-1)Delta^(l m)(k)hat(Theta)=mat(d^(* l m)_x (k)-i d^(* l m)_y (k),
d^(* l m)_0 (k)-d^(* l m)_z (k);-d^(* l m)_0 (k)-d^(*l m)_z (k),
-d^(*l m)_x (k)-i d^(* l m)_y (k)) $
と表せる。時間反転対称性があるとき$hat(Theta)^(-1)Delta^(l m)(k)hat(Theta)=Delta^(l m)(-k)$を満たすので,
$bold(d)$ベクトルは
$ d_0(k)=d^*_0(-k),\
bold(d)(k)=-bold(d)^* (-k) $
を満たす。これとギャップの持つ対称性$Delta (k)=-Delta^T (-k)$を用いると, 
$ d_0(k)&=d^T_0(-k)=d^*_0(-k),\
bold(d)(k)&=-bold(d)^T (-k)=-bold(d)^* (-k) $
となるので, 
$ d_0(k)&=d_0^dagger (k),\
bold(d)(k)&=bold(d)^dagger (k) $
を満たすことがわかる。よってギャップ関数は軌道対角部分で実数となる。
ギャップが時間反転を持たないとき, $bold(d)$ベクトルは対角成分に複素成分を持ってもよい。
この状態がいわゆるカイラル超伝導体に対応する。

松原周波数依存性を考慮した場合, 時間反転演算子によって$(k,i epsilon_m)arrow (-k,-i epsilon_m)$と変換されるため, ギャップ関数は
$ Delta^(l m)_(sigma, sigma ')(k,i epsilon_m)=sigma sigma' Delta^(* l m)_(macron(sigma), macron(sigma)')(-k,-i epsilon_m) $<TSRgap>
と言う対称性をみたす。これと交換関係からギャップ関数は
$ Delta^(l m)_(sigma, sigma ')(k,i epsilon_m)=-sigma sigma' Delta^(* m l)_(macron(sigma)', macron(sigma))(k,i epsilon_m) $
を満たす。この条件を用いるとsingletのギャップ関数は
$ Delta^(l m)_(s)(k,i epsilon_m)&=Delta^(* m l)_(s)(k,i epsilon_m) $
を満たし, 各$(k,i epsilon_m)$成分毎では軌道に対してエルミート対称性を持つ。
そのため時間反転を破らない限り対角成分は実数, つまり各バンド上のギャップは実数となる。
tripletにおいては$d$ベクトルが同様の性質を持つようになる。また異常グリーン関数も同じ性質を持つ。

ここで気をつけるべき事として, 実は時間反転対称性はそれだけではギャップ関数のパリティを決定しない。
あくまで満たすべき対称性は @psym,@TSRgap だけである。そのため, 周波数に対しての偶奇性が, 決まらないと, パリティの偶奇性も決まらない。
一般的には静的な安定性に作用するため, フェルミ準位で大きなギャップが開く偶周波数超伝導体が発現する。

線形Eliashberg方程式では, 異常グリーン関数は
$ F_(alpha,beta,s,s')(k,i epsilon)=sum_(a,b,sigma,sigma')G_(alpha,a,s,sigma)(k,i epsilon)Delta_(a,b,sigma,sigma')(k,i epsilon)G_(beta,b,s',sigma')(-k,-i epsilon) $
と書き表せる。この時, 時間反転に対するグリーン関数の対称性を用いる事で,
$ F_(alpha,beta,s,s')(k,i epsilon)=sum_(a,b,sigma,sigma')sigma' s'G_(alpha,a,s,sigma)(k,i epsilon)Delta_(a,b,sigma,sigma')(k,i epsilon)G^*_(beta,b,macron(s'),macron(sigma'))(k,i epsilon) $
と書き直すことが出来る。特にスピンの丸い系では
$ F_(alpha,beta)(k,i epsilon)=-sum_(a,b)G_(alpha,a)(k,i epsilon)Delta_(a,b)(k,i epsilon)G^*_(beta,b)(k,i epsilon) $
となる。
== 時間反転もうちょっと詳しい話
スピン非対角まで考慮したハミルトニアンは, パウリ行列で分解することができる。$hat(sigma)_x,hat(sigma)_y,hat(sigma)_z$に対応する成分を$hat(g)_x,hat(g)_y,hat(g)_z$と表すと, ハミルトニアンは
$ hat(H)(k)&=hat(H)_0(k)hat(sigma)_0+hat(bold(g))(k)dot hat(bold(sigma))\
&=hat(H)_0(k)hat(sigma)_0+hat(g)_x (k)hat(sigma)_x+hat(g)_y (k)hat(sigma)_y+hat(g)_z (k)hat(sigma)_z $
スピンを$2times 2$の行列で表すと上のハミルトニアンは
$ hat(H)(k)=mat(hat(H)_0 (k)+hat(g)_z (k), hat(g)_x (k)-i hat(g)_y (k); hat(g)_x (k)+i hat(g)_y (k),hat(H)_0 (k)-hat(g)_z (k)) $
これに時間反転演算子を作用させると
$ hat(Theta)^(-1) hat(H)(k) hat(Theta)&=hat(K)mat(0,1;-1,0)hat(H)(k)mat(0,-1;1,0)hat(K)\
&=mat(hat(H)^*_0 (k)-hat(g)^*_z (k), -(hat(g)^*_x (k)-i hat(g)^*_y (k)); -(hat(g)^*_x (k)+i hat(g)^*_y (k)),hat(H)^*_0 (k)+hat(g)^*_z (k)) $
となる。時間反転対称性があるときは, これが
$ hat(H)(-k)=mat(hat(H)_0 (-k)+hat(g)_z (-k), hat(g)_x (-k)-i hat(g)_y (-k); hat(g)_x (-k)+i hat(g)_y (-k),hat(H)_0 (-k)-hat(g)_z (-k)) $
と等価なので各成分は
$ hat(H)_0(k)&=hat(H)^*_0 (-k)=hat(H)^T_0 (-k)\
hat(bold(g))(k)&=-hat(bold(g))^* (-k)=-hat(bold(g))^T (-k) $
を満たす。
以前のスピンが丸い系での議論は$H_0 (k)$のみを考慮する場合に対応する。

スピンを考慮した系で感受率を求める場合, 既約感受率から, 相互作用バーテックス$V$を用いて, 相互作用を考慮した応答関数, 
$ hat(chi)&=(hat(I)-hat(chi)_0hat(V))^(-1)hat(chi)_0=hat(chi)_0(hat(I)-hat(V)hat(chi)_0)^(-1) $
を求めた後に, スピン/軌道感受率を以下の式から求める。
$ hat(chi)^(z z)_s&=1/4sum_sigma (hat(chi)_(sigma,sigma,sigma,sigma)-hat(chi)_(sigma,sigma,macron(sigma),macron(sigma)))\
hat(chi)^(plus.minus)_s&=hat(chi)_(arrow.t,arrow.b,arrow.b,arrow.t)\
hat(chi)^(minus.plus)_s&=hat(chi)_(arrow.b,arrow.t,arrow.t,arrow.b)\
hat(chi)_c&=1/2sum_sigma (hat(chi)_(sigma,sigma,sigma,sigma)+hat(chi)_(sigma,sigma,macron(sigma),macron(sigma))) $

== 空間の対称性
空間の対称性は, 原子の内部座標や軌道にも影響を与えるので, 色々面倒くさい。
一般的に同一サイト内の$d$-$d$間や$p$-$p$間の寄与のみを考慮する場合は空間反転対称性から来る符号の変化が打ち消されるので気にしなくて良いが, 
$d$-$p$間を考慮したり, 多サイトやnonsymmorphic(原子が反転中心に居ない)な系(+SOC)の模型だと複雑になるので注意。
空間反転の演算子に対して, $bold(r)arrow -bold(r), bold(p)arrow -bold(p) (bold(k)arrow -bold(k)), sigma arrow sigma$と変換される。
空間反転演算子を作用させると, 基底のパリティによって符号の変化が起こる。そのため先に述べたように偶奇の積になる$d-p$遷移の成分等は符号が変化する。
またユニットセル内の複数原子サイトを考慮した模型の場合, 同一サイト間のホッピングは, 基本並進ベクトルの空間対称性と整合するため, 問題無いが, 
別サイト間のホッピングは基本並進ベクトルだけで見ると, 空間反転が破れた状態に見える形になるので, 遠方のホッピングに関して整合性を取る必要がある。
同様のことは回転対称性や鏡面対称性を作用させた場合にも言え, 基底関数の$x<->y$の入れ替え等に対応した行列成分の入れ替えなどを考慮しないといけないので注意。

実際超伝導ギャップに関しては, パリティ混成が無い場合, 交換則の関係上軌道対角部分では, 
singletでは偶周波数偶パリティ or 奇周波数奇パリティ, tripletでは偶周波数奇パリティ or 奇周波数偶パリティが成り立つ。
非対角成分に関しては, 基底に取った軌道のパリティの積に依って符号が変化する。
例えば, $d$軌道内若しくは$p$軌道内のみを考える場合, 各軌道のパリティは同じなので, 
符号反転は起こらないため, ギャップの交換則から来る対称性を用いて, 
$ Delta(k,i epsilon_n)=-eta_k Delta^T (k,-i epsilon_n) $
が成り立つ。$eta_k$はギャップが持つパリティに対応する符号。
さらにスピンsinglet, triplet別々に見ると,
$ Delta_s (k,i epsilon_n)&=Delta_s^T (k,-i epsilon_n)\
Delta_t (k,i epsilon_n)&=Delta_t^T (k,-i epsilon_n) $
が成り立つ。ただしWannier関数ベースの模型の場合, 複数の軌道が混成しているため, パリティが綺麗に決まらない場合があるので注意。

== 周期性が無い場合の対称性
不純物を含む系や合金などは周期性が無いためハミルトニアンは$k$に対して非対角な項を持つ。そのため, ハミルトニアンは,
$ hat(H)=sum_(k,k')sum_(alpha,beta)H_(alpha,beta)(k,k')hat(c)^dagger_(k,alpha)hat(c)_(k',beta) $
のような形となる。これは, ハミルトニアンのエルミート性から
$ H_(alpha,beta)(k,k')=H^*_(beta,alpha)(k',k) $
と言う対称性を持つ(非エルミートはしらん)。$H$のエルミート性から来ているので,この対称性は波数以外が基底の場合でも成り立つ。
例えば周期性が無い時のホッピングが持つ対称性は
$ t_(alpha,beta)(r_i,r_j)=t^*_(beta,alpha)(r_j,r_i) $
となる。通常は周期性があるので$r=r_j-r_i$を用いてホッピングを表す事になる。すると, 先に述べたような対称性が成り立つ。
二体の相互作用$V$が持つ対称性も, 交換関係とエルミート性から与えられるので, 二体の相互作用項を
$ H'=sum_(1,2,3,4)V^(1234)c^dagger_1 c^dagger_4 c_3 c_2 $
と置くと 
$ V^(1234)=V^(4321)=(V^*)^(2143)=(V^*)^(3412) $
は周期性の有無は関係なく成り立つ。

普通の計算では, 格子サイトの周期性(並進対称性)と境界条件の周期性が一致するので, 軌道の自由度の分だけ(対称性に依っては縮退するが),独立した固有値が求まる。
不純物を考慮した計算の不純物抜きの様な, 格子の周期性と境界の周期性が異なる模型の計算では, 
(数値計算なので誤差は含むが)格子の周期性と境界の周期性の比に対応して固有値が縮退する。