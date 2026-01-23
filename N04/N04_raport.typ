#import "@preview/physica:0.9.7": *
#set heading(numbering: "1.")
#show title: set text(size: 20pt)
#show title: set align(center)
#set par(justify: true)
#set page(numbering: "1")


#let custom_link(href, body, color: purple) = {
  link(href)[#text(color)[_*#body*_]]
}

#title[
  Metody Numeryczne N04
]
#align(center)[Jakub Kurek]

= Wstęp
Celem zadania jest znalezienie interpolacji Lagrange'a $f_"approx" (x)$ funkcji:
$ f(x)=1/(1+x^2) $
na przedziale $x in [-5, 5]$.

Program liczący rozwiązania został napisany w C++23. Wykresy zostały stworzone w języku Python przy użyciu biblioteki matplotlib. Cały kod wykorzystywany do obliczeń oraz generowania wykresów znajduje się w repozytorium na #custom_link("https://github.com/Eg-Nickname/num_methods/tree/main/N04")[GitHub].

= Interpolacja Lagrange'a
Wzór interpolacji Lagrange ma postać:
$ f(x) = sum^n_(j=1) l_j (x) f_j $
gdzie:
$ l_(j)(x)=product_(j!=i)(x-x_i)/(x_j-x_i) $
== Liczenie wartości bezpośrednio z wzoru interpolacyjnego
Normlane obliczanie wartości $M$ punktów z wzoru interpolacyjnego Lagrange zajmowałoby $O(M dot N^2)$ operacji. Można jednak zastosować pewne optymalizację. Najepierw wyliczamy mianowniki wielomianów $l_j (x)$ i oznaczamy je $m_j$ gdyż nie zależą one od $x$ i są identyczne w każdym punkcie. Dla każdego punktu nastepnie obliczamy:
$ omega(x) = product_i^N x - x_i $
Wówczas wartość w punkcie liczymy ze wzoru:
$ y(x) = sum_(j=1)^n omega(x)/((x-x_j) dot m_j) $

Matematycznie jest to poprawne działanie, lecz w progrmaie należy sprawdzić czy nasz $x$ nie jest węzłem interpolacji, gdyż wówczas $omega = 0$ oraz $(x-x_j)=0$ co daje nam $0/0 = "NaN"$.

Koszt obliczenia punktów tą metoda to $O(N^2 + M*N)$.

== Liczenie współczynników wielomianu interpolacyjnego
Zakładamy, że żaden z węzłów interpolacji nie jest zerem. Jeśli, któryś z węzłów interpolacyjnych jest równy 0. To tworzymy wielomian Lagrange i wyliczamy jego wartość w punkcie $tilde(x) != 0$ w tym przypadku $tilde(x) = x_(i+1)/2$. Następnie usuwamy węzeł $x=0$ i zastępujemy go przez $tilde(x)$ oraz $tilde(f) = y(tilde(x))$.

#pagebreak()
Wówczas wielomian ma postać:
$ y(x) = sum^n_(j=1) l_(j)(x)f_j $

Współczynnik $a_0$ wielomianu obliczamy, licząc wartość w zerze:
$ a_0 = y(0) = sum^n_(j=1) l_(j)(0)f_j $
Koszt obliczenia $O(N^2)$, przy czym zapamiętujemy wielomiany $l_(j)(0)$.

Można zauważyć:
$ (y(x)-a_0)/x = dots + a_1 $
Nie możemy jednak podzielić przez x, gdyż obliczamy wartości w $0$. Możemy jednak skorzystać z interpolacji:
$
  (y(x_i) - a_0) / x_i
  = 1 / x_i (sum_(j=1)^n l_j (x_i) f_j - a_0)
  = 1 / x_i (sum_(j=1)^n delta_(i j) f_j - a_0)
  = (f_i - a_0) / x_i
$

Wówczas budujemy tabelę
#align(center)[
  #table(
    columns: (auto, auto, auto, auto, auto, auto),
    align: center + horizon,
    stroke: (x, y) => if x == 0 {
      (right: 1pt)
    } else {
      (left: 0.5pt)
    },

    $x_i$,
    $x_1$,
    $x_2$,
    $x_3$,
    $dots$,
    $x_n$,

    table.hline(),

    $f_i^((1))$,
    $(f_1 - a_0) / x_1$,
    $(f_2 - a_0) / x_2$,
    $(f_3 - a_0) / x_3$,
    $dots$,
    $(f_n - a_0) / x_n$,
  )
]

a na jej podstawie wielomian interpolacyjny:
$ y^((1))(x)= sum^n_(j=1) l_(j)(x)f^((1))_j $
Korzystamy z faktu iż znamy $l_(j)(0)$, ponieważ zależą one tylko od węzłów. Możemy obliczyć wartość:
$ a_1 = y^((1))(0) = sum^n_(j=1) l_(j)(0)f^((1))_j $

Uogólniając ten algorytm otrzymujemy:
$ a_k = sum^n_(j=1) l_(j)(0)f^((k))_j quad k = 0, dots, n-1 $
gdzie:
$
  f^((k))_j = cases(
    f_j & #h(3em) k=0,
    (f_j^((k-1))-a_(k-1))/(x_j) & #h(3em) k!= 0
  )
$

Obliczenie każdego kolejnego $a_k$ dokonujemy w czasie $O(N)$, więc obliczenie współczynników wielomianu zajmuje w sumie $O(N^2)$.

#pagebreak()
== Obliczanie wartości wielomianu
Dla efektywnego obliczania wartości wielominau korzystamy z algorytu Hornera:
#block(
  inset: 10pt,
  radius: 4pt,
  [
    #set enum(indent: 10pt)
    + $P = a_n$
    + $k = n$
    + *while* $k > 0$:
      - $k = k - 1$
      - $P = P dot z + a_k$
    + *return* $P$
  ],
)

= Interpolacja funkcji
Celem zadania było znalezienie optymalnej interpolacji, która zminimalizuje błąd przybliżenia $Delta f = max{abs(f - f_"approx")}$ dla gęstej siatki punktów, które będą znajdować się równiez pomiędzy węzłami interpolacji. Dla zadanej funkcji $f(x)=1/(1+x^2)$ dla $x in [-5,5]$. Przeanalizowane były dwie metody dobierania węzłów:
- węzły równoodległe: $x_n =-5 + 10n/(N-1)$ dla $n in 0, dots, N$
- węzły Czebyszewa: $x_n =- 5 cos((n pi)/(N-1))$ dla $n in 0, dots, N$

== Znajdowanie optymalnej ilości węzłów interpolacyjnych
Dla każdej metody doboru węzłów wygenerowane zostały interpolacje dla $n in {2,3, dots 256}$ oraz została policzona wartość $Delta f(x)$ korzystając z siatki $1024$ równorozmieszczonych punktów na przedziale $[-5, 5]$, gdyż w węzłach interpolacyjnych błąd zawsze wynosi $0$, natomiast różnice w wartościach funkcji i interpolacji między węzłami determinujeą jakość interpolacji.

== Porównanie metod doboru wezłów interpolacyjnych przy użyciu współczynników wielomianu
Dla węzłów równoodległych najbardziej optymalny był wielomian dla $N=8$ z $Delta f(x) = 0.247339$.

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/equal_plot.jpg"),

  caption: [Interpolacji wielomianem, którego współczynniki zostały policzone z równoodległych węzłów
  ],
) <polynomial_equal_points_plot>

#pagebreak()
Dla węzłów Czebyszewa najbardziej optymalny był wielomian dla $N=30$ z $Delta f(x) = 0.00616472$

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/czebyszew_plot.jpg"),

  caption: [Interpolacji wielomianem, którego współczynniki zostały policzone z węzłów Czebyszewa
  ],
) <polynomial_czebyszew_points_plot>

== Porównanie metod doboru wezłów interpolacyjnych, korzystając bezpośrednio ze wzoru interpolacyjnego
Dla węzłów równoodległych również najbardziej optymalny był wielomian z $N=8$ i odpowiadający mu @polynomial_equal_points_plot. Dużo lepszą dokładność jednak dostaliśmy korzystając z węzłów Czebyszewa gdzie najbardziej optymalna aproksymacja była dla $N=222$ węzłów a $Delta f(x) = 1.30104 dot 10^(-18)$.

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/direct_czebyszew_plot.jpg"),

  caption: [Interpolacja wielomianem Lagrange'a
  ],
) <direct_czebyszew_points_plot>



== Problemy z ilością węzłów
Interpolację wielomianem Lagrange'a dla równoodległych punktów bardzo szybko zabijają oscylacjie Rungego, czyli gwałtowne zmiany między węzłami interpolacyjnymi.

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/runge_plot.jpg"),

  caption: [Wizualizacja oscylacji Rungego
  ],
) <runge_points_plot>
Jak widać na @runge_points_plot dla $x in [-3, 3]$ funkcja ta dużo lepiej przybliża funkcję $f$ lecz oscylacje dla $x in [-5, -4] union [4, 5]$ powodują, że nasz błąd $Delta f =0.556723$ urósł ponad dwukrotnie. Dalsze zwiększanie ilości węzłów jeszcze bardziej potęguje siłę oddziaływania oscylacji Rungego co skutkuje jeszcze większym błędem $Delta f$.

Rozwiązaniem dla szybko pojawiających się oscylacji Rungego dla równoodlęgłych węzłów jest skorzystanie z węzłów Czebyszewa, których zagęszczenie zwiększa się na granicach przedziału i opóźnia się ich pojawienie. Pozwala to zwiększyć liczbę węzłów i lepiej dopasować interpolowaną funkcję.

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/errors.jpg"),

  caption: [Wizualizacja błędów interpolacji dla różncyh ilości węzłów
  ],
) <errors_plot>

@errors_plot pokazuje jak zmieniają się błędy $Delta f$ w zależności od ilości węzłów interpolacji. Dla obu metod z współczynnikami wielomianu najpierw błąd maleje, aż w pewnym momencie zaczyna rosnąć logarytmicznie do nieskończoności, odpowiednio są to $N=10$ i $N=30$ dla węzłów równoodległych oraz węzłów Czebyszewa. Dla węzłów równoodległych brak zwiększenia dokładności  powodują oscylacje Rungego, natomiast dla węzłów Czebyszewa problemem staje się dokładność numeryczna bezpośredniego wyliczania współczynników wielomianu. Dopiero wyliczanie punktów bezpośrednio z wzoru Lagrange dla każdego punktu rozwiązuje ten problem. Zmniejszanie dokładności kończy się dla $N=222$ i nastepuje stagnacja spowodowana skończoną dokładnością numeryczną, od tego punktu wartości oscylują wokoło wartości $2 dot 10^(-18)$.

= Podsumowanie
Jeśli to tylko możliwe nie powinniśmy używać węzłów równo rozmieszonych i zastąpić je węzłami Czebyszewa, które zmniejszają oscylacje Rungego. Dla dużych wielomianów znaczenie przejmują błędy numeryczne i to one psują interpolację w przyadku węzłów Czebyszewa. Najlepszym rozwiązaniem jest liczenie interpolowanych punktów bezpośrednio ze wzoru Lagrange'a wykorzystując jego właściwości i zmniejszając koszt obliczeń do $O(N^2 + M dot N)$ gdzie dla dostateczine dużych ilość interpolowanych punktów metoda ta ma złożoność $O(M dot N)$, gdyż człon z $N^2$ staje się zaniedbywalnie mały.
