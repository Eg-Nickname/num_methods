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
  Metody Numeryczne N03
]
#align(center)[Jakub Kurek]

= Wstęp
Celem zadania jest znalezienie dwóch największych oraz czterech najmniejszych co do modułu wartości własnych oraz ich wektorów własnych macierzy:
#grid(
  columns: (1fr, 1fr),
  align: horizon,
  [
    $
      A_(i j) = cases(
        -2/h^2 + V_i & "dla" i = j,
        1/h^2 & "dla" i - j = plus.minus 1,
        0 & "w przeciwnym przypadku"
      )
    $
  ],
  $
    h = 20/(N-1) quad V_i=(i h - 10)^2
  $,
)


Program liczący rozwiązania został napisany w C++23 przy użyciu biblioteki #custom_link("https://libeigen.gitlab.io/")[Eigen] w wersji 5.0.1. Wykresy zostały stworzone w języku Python przy użyciu biblioteki matplotlib. Cały kod wykorzystywany do obliczeń i rysowania wykresów znajduje się w repozytorium na #custom_link("https://github.com/Eg-Nickname/num_methods/tree/main/N03")[GitHub].

= Metoda potęgowa
Dla naszej macierzy $A$  wiemy, że jest ona symetryczna $A^TT=A in RR^(N times N)$. Dzieki temu wiemy, że jest ona diagonalizowalna, ma rzeczywiste wartości własne oraz jej unormowane wektory własne tworzą bazę ortonormalną w $RR^N$. Oznaczmy tą bazę ${e_i}^N_(i=1)$ wówczas $A e_i = lambda_i e_i$. Przyjmujemy dodatkowo, że wartości własne są dodatnie i uporządkowane $lambda_1 > lambda_2 > dots > lambda_N > 0$.
Dla dowolnego wektora $y in R^N$ istnieje rozkład w bazie:
$ y = sum_i^N Beta_i e_i $
Mamy wówczas:
$
  A y = A sum_i^N Beta_i e_i = sum_i^N Beta_i A e_i = sum_i^N Beta_i lambda_i e_i \
  A^k y = sum_i^N Beta_i lambda_i^k e_i
$
Dla dostatecznie dużych $k$ wyraz z $lambda_1^k$ będzie dominował w sumie, a reszta współczynników będzie zaniedbywalna. Dla dostatecznie dużych $k$ prawa strona będzie dążyć do wektora proporcjonalnego do $e_i$ czyli wektora własnego największej wartości własnej.

Dla $norm(y_1)=1$ iteracja:
$
      z_k & = A y_k \
  y_(k+1) & = z_k / norm(z_k)
$
zbiega się do unormowanego wektora własnego $A$ odpowiadającego największej wartości własnej.
Dla $y_k tilde.eq e_1$ (wartości przestają się zauważalnie zmieniać). Unormowaną wartość własną obliczamy $norm(lambda_1) = norm(z_k) = norm(A y_k)$.

== Kolejne wartości własne
Wektor $y_k$ nie będzie zbiegał do największej wartości własnej (każdej wartości własnej większej niż $lambda_i$). Wtedy i tylko wtedy gdy współczynniki $Beta_j=0$ dla $j < i$. Dlatego trzeba zapewnić, żeby $y_1$ oraz kolejne $y_k$ były prostopadłe do wcześniej obliczonych wektorów własnych.

W arytmetyce dokładnej ortogonalizacja musi nastąpić tylko dla $y_1$, lecz ze względu na błędy arytmetyki nie dokładnej. Trzeba ortogonalizować każdy wektor $y_k$.

Iteracja wówczas dla i-tej wartości własnej wygląda:\
Założenia: $norm(y_1)=1 and forall_(j<i)e^TT_j y_1 = 0$

$
      z_k & = A y_k \
      z_k & = z_k - sum_(j<i)e_j(e_j^TT z_k) \
  y_(k+1) & = z_k / norm(z_k) \
$

== Najmniejsze wartości własne
Dla macierzy $A^TT=A in RR^(N times N)$ najmniejsza wartość własna odpowiada $1 / max{mu_i}$ macierzy $A^(-1)$.
Metoda potęgowa wówczas przybiera postać:
$
      z_k & = A^(-1) y_k \
  y_(k+1) & = z_k / norm(z_k)
$

Zapis $A^(-1)y_k = z_k$ rozwiązujemy jako $A z_k = y$. Dla naszej macierzy $A$ ze względu na jej trójdiagonalna postać możemy skorzystać z algorytmu Thomasa, przez co obliczanie kolejnych wektorów $z_n$ wykonujemy w czasie $O(N)$.

= Metoda Rayleigha
Metoda Rayleigha jest modyfikacją metody potęgowej, która wykorzystuje fakt iż wartości własne macierzy $A - sigma II$ wynoszą $lambda_i - sigma$, gdzie $lambda_i$ są wartościami własnymi macierzy $A$ dla odpowiadających wektorów własnych. Szukając największej wartości własnej $lambda_i$ tak naprawdę szukamy najmniejszych wartości własnych macierzy $A - sigma II$, ponieważ $max{lambda_i}-sigma tilde.eq 0$

Dla najlepszych efektów najlepiej zacząć rozwiązywanie zwykła metoda potęgowa na pare iteracji, żeby później aktualizować lepiej przybliżoną wartość $sigma$.

Po zmianie metody iteracja przybiera postać:
$
    sigma & = y_k^TT A y_k \
      z_k & = (A- sigma II)^(-1) y_k \
      z_k & = z_k - sum_(j<i)e_j(e_j^TT z_k) \
  y_(k+1) & = z_k / norm(z_k) \
$

W zależności od tego czy chcemy znaleźć największa czy najmniejsza wartość własną pre-iterujemy metoda potęgowa lub odwrotną metoda potęgowa.

#pagebreak()
= Rozwiązanie
Przyjmując $N=500$ metody dają rozwiązania
#grid(
  columns: (1fr, 1fr),
  align: center,
  [
    #figure(
      kind: "table",
      supplement: [Tabela],
      table(
        columns: (auto, auto),
        inset: 8pt,
        align: horizon + center,
        table.header([*$min{abs(lambda_i)}$*], [*Wartość*]),
        $1/mu_1$, [$0.490098$],
        $1/mu_2$, [$0.713051$],
        $1/mu_3$, [$-0.806307$],
        $1/mu_4$, [$-1.75069$],
      ),

      caption: [Najmniejsze wartości własne
      ],
    ) <smallest_evals>

  ],
  [
    #figure(
      kind: "table",
      supplement: [Tabela],
      table(
        columns: (auto, auto),
        inset: 8pt,
        align: horizon + center,
        table.header([*$max{abs(lambda_i)}$*], [*Wartość*]),
        $lambda_1$, [$-2489.01$],
        $lambda_2$, [$-2487.01$],
      ),

      caption: [Największe wartości własne
      ],
    ) <biggest_evals>
  ],
)

Kluczową różnica w obu metodach jest szybkośc zbieżności:
#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (1fr, 2fr, 2fr, 2fr),
    inset: 8pt,
    align: horizon + center,
    table.header(
      [*$min{abs(lambda_i)}$*],
      [*Ilośc iteracji potęgowa*],
      [*Ilośc iteracji Rayleigha*],
      [*Ilośc iteracji Rayleigha minus preiteracje*],
    ),
    $1/mu_1$, [$21$], [$6$], [$3$],
    $1/mu_2$, [$89$], [$9$], [$6$],
    $1/mu_3$, [$13$], [$5$], [$2$],
    $1/mu_4$, [$22$], [$5$], [$2$],
  ),

  caption: [Zbieżność dla najmniejszych wartości własnych
  ],
) <smallest_iters>

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/smallest_ev_cmp_plot.jpg"),

  caption: [Zbieżność dla najmniejszych wartości własnych
  ],
) <smallest_iters_plot>


#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (1fr, 2fr, 2fr, 2fr),
    inset: 8pt,
    align: horizon + center,
    table.header(
      [*$max{abs(lambda_i)}$*],
      [*Ilośc iteracji potęgowa*],
      [*Ilośc iteracji Rayleigha*],
      [*Ilośc iteracji Rayleigha minus preiteracje*],
    ),
    $lambda_2$, [$1941$], [$129$], [$126$],
    $lambda_1$, [$2127$], [$111$], [$108$],
  ),

  caption: [Zbieżność dla największych wartości własnych
  ],
) <biggest_iters>

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/biggest_ev_cmp_plot.jpg"),

  caption: [Zbieżność dla największych wartości własnych
  ],
) <biggest_iters_plot>

W zależności od tego jak dobrze jest zbieżna metoda dla danej macierzy i kolejnych wartości własnych  metoda Rayleigha potrzebowała 15-20 razy mniej iteracji niż metoda potęgowa @biggest_iters oraz 10-15 razy mniej iteracji niż odwrotna metoda potęgowa @smallest_iters. Ze względu na dość szybką zbieżność odwrotnej metody potegowej dla zadanej macierzy połowa wszystkich iteracji w metodzie Rayleigha zostaje wykorzystana do naprowadzenia metodą potęgową. Dopiero porównanie ich od momemtu, w którym iteracje się róznią daje nam rzeczywiste porównanie różnic w zbiezności.

= Podsumowanie
Metoda Rayleigha jest zdecydowania bardziej efektywna niż metoda potęgowa pod względem wykonywanych iteracji, chociaż trzeba zwrócić uwagę, że dla niej rozkład macierzy jest wykonywany w kazdym kroku iteracji, a nie jak w metodzie potęgowej tylko na początku. W tym przypadku rozkład i rozwiązywanie przy uzyciu algorytmu Thomasa jest szybkie, lecz nie zawsze zachodzi to w ogólnosci.

Drugim ważnym problem metody Rayleigha jest to, że jeśli nie wykonamy dostatecznie dużo iteracji metodą potęgową lub dla szukanej wartości wartości własnej $lambda_i$ zachodzi $lambda_(i+1) / lambda_i<< 1$ to metoda ta, może zbiec do fałszywego maksimum (lub minimum), gdyż najlepiej ona znajduje wartość własną najbliższą $sigma$, co w przypadku złego uwarunkowania może prowadzić do złych szacowań $sigma$.




