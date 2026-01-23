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
  Metody Numeryczne N06
]
#align(center)[Jakub Kurek]

= Wstęp
Celem zadania jest znalezienie miejsc zerowych wielomianu:
$ 243x^7 - 486x^6 + 783x^5 - 990x^4 + 558x^3 - 28x^2 - 72x + 16 $

Program liczący rozwiązania został napisany w C++23. Cały kod wykorzystywany do obliczeń znajduje się w repozytorium na #custom_link("https://github.com/Eg-Nickname/num_methods/tree/main/N06")[GitHub].

= Metoda Laguerra
Dla wielomianów $P_n(z)$ metoda znajdywania miejsc zerowych jest metoda Laguerre'a, która zadana jest iteracją:
$
  z_(i+1) = z_i - frac(n P_n (z_i), P'_n (z_i) plus.minus sqrt((n-1) ((n-1) [P'_n (z_i)]^2 - n P_n (z_i) P''_n (z_i))))
$
gdzie znak w mianowniku wybieramy tak, aby moduł mianownika był większy.

== Obniżanie stopnia wielomianu i wygładzanie
Aby unikać zbiegania metody do poprzednio znalezionego miejsca zerowego obniżamy stopień wielomianu, czyli znajdujemy faktoryzacje $P_n (z) = (z - z_1)P_(n-1) (z)$, gdzie $z_1$ jest znalezionym miejscem zerowym.

Następnie szukamy miejsca zerowego wielomianu $P_(n-1) (z)$. Ze względu na istnienie zaburzeń znalezione miejsce zerowe poprawiamy za pomocą pełnego, nie wydzielonego wielomianu. Dla znalezionego miejsca zerowego $tilde(z)_2$ wielomianu $P_(n-1)$, wykonujemy wygładzenie poprzez zapoczątkowanie metody Laguerre'a miejscem $tilde(z)_2$ jako warunkiem początkowym dla wielominau P_n.

Dopiero dla wygładzonego miejsca $z_2$ znajdujemy faktoryzację $P_(n-1) = (z-z_2) P_(n-2) (z)$

Aby znaleźć faktoryzację należy rozwiązać układ:
$
                 b_(n-1) & = a_n \
  -z_0 b_(n-1) + b_(n-2) & = a_(n-1) \
  -z_0 b_(n-2) + b_(n-3) & = a_(n-2) \
                         & dots.v \
          -z_0 b_2 + b_1 & = a_2 \
          -z_0 b_1 + b_0 & = a_1 \
                -z_0 b_0 & = a_0
$
Po uproszczeniu i korzystaniu z znanych nam informacji otrzymujemy układ równań, który można rozwiązać metoda forward substitution:
$
  mat(
    1, , , , , ;
    -z_0, 1, , , , ;
    , -z_0, 1, , , ;
    , , dots.down, dots.down, , ;
    , , , -z_0, 1, ;
    , , , , -z_0, 1;
  ) mat(b_(n-1); b_(n-2); b_(n-3); dots.v; b_1; b_0;) = mat(a_n; a_(n-1); a_(n-2); dots.v; a_2; a_1;)
$

Obniżamy kolejne wielomiany $P_(n-k)$ aż do uzyskania wielomianu $P_2$, który rozwiązujemy ze znanych nam wzorów analitycnzych.

= Rozwiązanie
#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto),
    inset: 8pt,
    align: horizon + center,
    // table.header([*$i$*], [*z_i*]),
    [*i*], $z_1$, $z_2$, $z_3$, $z_4$, $z_5$, $z_6$, $z_7$,
    [*wartość*], [$0.333333$], [$0.666662$], [$0.666665$], [$-0.333333$], [$0.666665$], [$1.41421 i$], [$-1.41421 i$],
  ),

  caption: [Miejsca zerowe wielomianu
  ],
) <polynomial_solutions>

= Podsumowanie
Metoda Laguerra jest efektywną metodą znajdowania miejsc zerowych wielomianów. Posiada zbieżność sześcienna dla pojedyńczych pierwiastków oraz liniowa dla wielokrotnych. Znajduje ona równiez zespolone pierwiaski, dzięki czemu dla każdegowielomianu możemy znaleźć jego pełną faktoryzację. Metoda ta jest równiez praktycznie zawsze zbieżna do pierwiastku co praktycznie niweluje potrzebe specjalnego dobierania miejsca, z którego zaczynamy iterację.
