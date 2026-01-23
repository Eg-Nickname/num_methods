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
  Metody Numeryczne N07
]
#align(center)[Jakub Kurek]

= Wstęp
Celem zadania jest znalezienie rozwiązania układu równań:
$
  cases(
    - (u_(n-1) - 2u_n + u_(n+1)) / h^2 + 2u_n (u_n^2 - 1) & = 0 #h(3em) n = 1 dots N-1,
    u_0 & = 0,
    u_N & = 1,
  )
$
gdzie $h = 20/(N-1)$.

Program liczący rozwiązania został napisany w C++23 przy użyciu biblioteki #custom_link("https://libeigen.gitlab.io/")[Eigen] w wersji 5.0.1. Wykresy zostały stworzone w języku Python przy użyciu biblioteki matplotlib. Cały kod wykorzystywany do obliczeń i rysowania wykresów znajduje się w repozytorium na #custom_link("https://github.com/Eg-Nickname/num_methods/tree/main/N07")[GitHub].

= Wielowymiarowa metoda Newtona
Niech funkcja $g: RR^N -> RR^N$ będzie funkcją klasy co najmniej $C^1$ (wszystkie pochodne cząstkowe są funkcjami ciągłymi zmiennej $x$). Rozważamy równanie:
$ g(x) = 0 $
W zależności od parametrów układ równań może mieć różną ilość rozwiązań.

Rozwijając funkcję $g$ w szereg Taylora do pierwszego rzędu otrzymujemy:
$ g(x + delta x) tilde.eq g(x) J delta x $
gdzie J jest jakobianem funkcji g:
$ J(x)_(i j) = partial g_i/ partial x_j bar_x $
Aby znaleźć się w punkcie spełniającym równanie $g(x + delta x) = 0$:
$ delta x = - J^(-1) g(x) $

Prowadzi to do iteracji:
$ x_(k+1) = x_k - J^(-1)(x_k)g(x_k) $

Jakobian trzeba obliczać w każdym kroku iteracji. (rozwiązanie innego układu równań liniowych).

Zamiast rozwiązywać układ równań algebraicznych można rozwiązać funkcję $G: RR^n -> RR$, ponieważ minimalizacja jest łatwiejsza:
$ G(x)=1/2 norm(g(x))^2 = 1/2 (g(x))^TT g(x) $

Globalne minimum funkcji $G=0$ odpowiada rozwiązaniu układu równań, jednak funkcja $G$ może mieć wiele minimów lokalnych (lub globalne minimum może nie istnieć)

Rozwiązując równanie $g(x)=0$ metodą Newtona krok iteracji wynosi:
$ delta x = -J^(-1) g $
oraz mamy:
$ (partial G)/ (partial x_i) = sum_j J_(j i)g_j $
więc:
$ nabla G = J^TT g $

Funkcja $G$ po wykonaniu kroku Newtona:
$ (nabla G)^T delta x = g^T J (-J^(-1)) g = -g^T g < 0 $
Co pokazuje, że kierunek kroku Newtona jest lokalnym kierunkiem spadku G.
Jednak przesunięcie się o pełną długość kroku Newtona nie musi prowadzić do spadku G. Postępujemy więc:

Obliczamy $delta x$ oraz ustawiamy wagę początkową $w = 1$.

Dopóki nie znajdziemy poprawy ($G(x_"test") < G(x_i)$):
- $x_"test" = x_i + w delta x$
- Jeśli $G(x_"test") < G(x_i)$:
  - $x_(i+1) = x_"test"$
  - Przechodzimy do kolejnej iteracji
- W przeciwnym przypadu:
  - $w = w * alpha quad "gdzie" alpha < 1$
  - Wracamy do obliczania $x_"test"$.

Metoda ta jest zawsze zbieżna do jakiegoś minimum funkcji $G$, ale niekoniecznie do jej minimum globalnego. Jeśli znajdziemy minimum lokalne $G_"min" > epsilon$, rozpoczynamy metodę z innym warunkiem początkowym. Jeśli kilka warunków początkowych nie daje pożądanego rozwiązania należy się poddać.

Im lepszy warunek początkowy tym szansa na znalezienie rozwiązania jest większa. Należy zainwestować naszą wiedzę o funkcji $g$ w jego znalezienie.
== Jakobian
W naszym przyadku Jakobian pochodnych cząstkowych ma postać trójdiagonalną, gdyż nieliniowy składnik zależy tylko od kolejnych 3 zmiennych (inne pochodne cząstkowe będą się zerować).
#let h2 = $h^2$
$
  J = mat(
    1, 0, 0, dots, 0;
    -1/h2, (partial f_1) / (partial u_1), -1/h2, , dots.v;
    0, dots.down, dots.down, dots.down, 0;
    dots.v, , -1/h2, (partial f_(N-1)) / (partial u_(N-1)), -1/h2;
    0, dots, 0, 0, 1
  )
$
Gdzie dla elementów na diagonali
$
  (partial f_n) / (partial u_n) = frac(partial, partial u_n) [ - (u_(n-1) - 2u_n + u_(n+1)) / h^2 + 2u_n (u_n^2 - 1) ] = 2/h^2 + 6u_n^2 - 2
$
Struktura trójdiagonalna macierzy pozwala nam rozwiązać układ równań $J delta u = g(u)$ algorytmem Thomasa w czasie $O(N)$.

= Rozwiązania otrzymane tłumioną wielowymiarową metoda Newtona
#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/solution_plot.jpg"),

  caption: [Rozwiązanie równania
  ],
) <nonlinear_eq_sol_plot>

= Podsumowanie
Tłumiona wielowymiarowa Newtona w zależności od warunku początkowego, może osiągać różne minima lokalne, ze względu na nieliniowość członu $- (u_(n-1) - 2u_n + u_(n+1)) / h^2 + 2u_n (u_n^2 - 1)$. W zależność od tego jak wyglądał początkowy wektor przybliżeń $u$ metoda ta dała inne wyniki. Dla poczatkowego wektora $u$ będacego zapoczatkowanego jedynkami, metoda dośc szybko zbiega do ustabilizowania wyników kolejnych wartości $u_n$ jako jeden. Podobnie jest dla wektora losowego, który najpier wykonuje oscylację do wartości $-1$, żeby następnie ustabilizować się w $1$. Dla liniowo zapoczątkowanego wektora $u_n = 1/(N-1)$, wartości rozwiązania przez dośc długi czas znajdują się w $-1$, żeby na sam koniec zbiec do wartości $1$. Najciekawszy jest wektor $u$ będący początkowo zerami oprócz ostatniego $u_N=1$, wykonuje on 2 oscylacje do momentu ustabilizowania się końcowych wartości $u_n=1$. @nonlinear_eq_sol_plot dobrze pokazuje, że nasze równanie posiada różne rozwiązania, które spełniają założoną dokładnośc $10^(-6)$.
