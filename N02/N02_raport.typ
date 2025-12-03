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
  Metody Numeryczne N02
]
#align(center)[Jakub Kurek]

= Wstęp
Celem zadania jest sprawdzenie róznych metod rozwiązania układu równań $A u=b$

#let diag = $(-2) / h^2$
#let other = $1 / h^2$
#grid(
  columns: (1fr, 0.5fr, 1fr),
  align(center)[
    $
      A^(N times N) = mat(
        #diag, #other, , , #other ;
        #other, #diag, #other, , ;
        , #other, #diag, #other, , ;
        , , , dots.down, ;
        #other, , , #other, #diag
      )
    $],
  align(center + horizon)[
    $ h=2/(N-1) $
  ],
  align(center + horizon)[
    $b_n = cos((4 pi (n-1)) / N)$
  ],
)
Program liczący rozwiązania został napisany w C++23 przy użyciu biblioteki #custom_link("https://libeigen.gitlab.io/")[Eigen] w wersji 5.0.1. Wykresy zostały stworzone w języku Python przy użyciu biblioteki matplotlib. Cały kod wykorzystywany do obliczeń oraz generowania wykresów znajduje się w repozytorium na #custom_link("https://github.com/Eg-Nickname/num_methods/tree/main/N02")[GitHub].

== Konfiguracja sprzętowa
- Procesor: Intel i7-8650U (8) @ 4.200GHz
- Pamięć RAM: 32GB DDR4 2133MHz
- System operacyjny: Arch Linux 6.17.5-arch1-1

== Metodyka testowania
Wszystkie testy były wykonywane na odciążonym procesorze bezpośrednio z powłoki. Dla obliczeń trwających poniżej 20ms testy były wykonywane 10 krotnie a czas był uśredniany. Pomiary czasu były wykonywane w mikrosekundach w celu lepszego zobrazowania czasu obliczeń dla małych macierzy oraz optymalnych metod. Po wprowadzeniu optymalizacje wszystkie dalsze testy były kompilowane przy użyciu kompilatora clang++ w wersji 21.1.4 z flagami ```-Wall -Wextra -std=c++23 -pedantic -I ./../eigen-5.0.1/ -O3 -march=native -DNDEBUG```.
= Wstępne optymalizacje
Najważniejszą optymalizacją poprawiające wyniki obliczeń było zastosowanie odpowiednich flag kompilatora:

- ```-O3``` - poziom optymalizacji kodu przez kompilator.  @o3cmp pokazuje róznice czasu wykonania zwiększoną około 10-krotnie. W przypadku funkcji bibliotecznych ta różnica zmniejsza się do około 2-krotności.
#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/O3_comparison_100k.jpg", width: 70%),
  caption: [Porównanie algorytmu Thomasa z wzorem Shermana-Morrisona z optymalizacją -O3 i bez
  ],
) <o3cmp>

#box[
  - ```-march=native``` - możliwośc wykonywania instrukcji wektoryzujących AVX. @march_cmp pokazuje róznice czasu obliczeń wybranych algorytmów z bilbioteki Eigen po włączeniu nstrukcji AVX.
  #figure(
    kind: "chart",
    supplement: [Wykres],
    image("./figures/march_comparison_4k.jpg", width: 75%),
    caption: [Porównanie algorytmów z biblioteki Eigen z użyciem instrukcji AVX i bez
    ],
  ) <march_cmp>
]

#box[
  - ```-DNDEBUG``` - wyłączenie assercji. Eigen dla macierzy o dynamicznych rozmiarach sprawdza poprawność ich rozmiaru przy użyciu assercji. Nie ma wielkiego wpływu na pojedyńcze działanie programu dla dużych macierzy, lecz może wpływać kiedy wykonujemy wiele działań na wektorach lub macierzach. Zalecane przez dokumentacje.
]

\
= Rozwiązania metodami gęstymi
Biblioteka Eigen oferuje moduł służący do rozwiązywania układów równań przy użyciu macierzy gęstych. Wybrane zostały metody:
- Full Pivot LU
- Partial Pivot LU
- Householder Full Pivot QR
- Householder Partial Pivot QR

Wszystkie metody uzywały wspólnej funkcji ```cpp gen_desnse_A()``` (@gen_dense_mat_code)  do generowania gęstej macierzy oraz funkcji ```cpp gen_b_vector()``` (@gen_b_vec_code) do tworzenia wektora rowiązania.
Funkcja ```cpp gen_dense_B``` była wykorzystywana do sprawdzania działania własnej implementacji alorytmu Thomasa, co wymusiło rozbicie tworzenia macierzy gęstej na dwa etapy.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::MatrixXd gen_dense_B(long N) {
    // matrix must be greater than 1
    // assert(N != 0);
    float64_t h = 2 / ((float64_t)N - 1.0);
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(N, N);

    for (long x = 1; x < N - 1; x++) {
      mat(x, x - 1) = 1.0 / (h * h);
      mat(x, x + 1) = 1.0 / (h * h);
      mat(x, x) = -2.0 / (h * h);
    }
    mat(0, 0) = -3.0 / (h * h);
    mat(N - 1, N - 1) = -3.0 / (h * h);

    // Diagonal fix
    mat(0, 1) = 1.0 / (h * h);
    mat(N - 1, N - 2) = 1.0 / (h * h);

    return mat;
  }

  // A = B + uvT
  Eigen::MatrixXd gen_dense_A(long N) {
    auto mat = gen_dense_B(N);
    float64_t h = 2 / ((float64_t)N - 1.0);

    // Distortion from 3 diagonal
    mat(0, N - 1) = 1.0 / (h * h);
    mat(N - 1, 0) = 1.0 / (h * h);

    mat(0, 0) += 1.0 / (h * h);
    mat(N - 1, N - 1) += 1.0 / (h * h);
    return mat;
  }

  ```,
  caption: [Generowanie gęstej macierzy pełnej oraz samej wersji diagonalnej
  ],
) <gen_dense_mat_code>

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd gen_b_vector(long N) {
    Eigen::VectorXd b = Eigen::VectorXd(N);
    b.setZero();
    for (long n = 1; n <= N; n++) {
      float64_t val = std::cos((4.0 * std::numbers::pi * (n - 1)) / (double)N);
      b(n - 1) = val;
    }

    return b;
  }
  ```,
  caption: [Generowanie wektora rozwiązania
  ],
) <gen_b_vec_code>


== Full Pivot LU
Pierwszą z omawianych metod jest rozwiązanie układu równań przy użyciu dekompozycji LU z pełnym pivotingiem. Analizując wybrane dane z @full_piv_lu_data  metoda rośnie w czasie $O(N^3)$ co jest zgodne z teoretyczna złozonością algorytmu. Jest to bazowy wynik, do którego kolejne metody będą się odnosić.
#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd solve_d_mat_fullpiv_lu(long N) {
    Eigen::MatrixXd d_mat = gen_dense_A(N);
    Eigen::VectorXd b = gen_b_vector(N);

    Eigen::VectorXd u = d_mat.fullPivLu().solve(b);
    return u;
  }
  ```,
  caption: [Rozwiązanie równania przy użyciu rozkładu LU z pełnym pivotem dla podanego N],
) <solve_full_lu_code>

#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/dense_full_pivot_lu.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N przy użyciu fullPivLU
  ],
) <full_piv_lu_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*], $324μ s$, $254 m s$, $2947 m s$, $10s$, $25s$, $49s$, $86s$, $134s$, $200s$, $283s$, $349s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu fullPivLU
  ],
) <full_piv_lu_data>

== Partial Pivot LU
Kolejną z omawianych metod jest rozwiązanie układu równań przy użyciu dekompozycji LU z częściowym pivotingiem. Kosztem precyzji numerycznej zyskujemy krótszy czas wykonania. Analizując wybrane dane z @par_piv_lu_data  metoda rośnie w czasie $O(N^3)$ , lecz w porównaniu z czasami dla pełnego pivotingu (@full_piv_lu_data) czas rozwiązywania układu równań jest ponad 12-krotnie mniejszy.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd solve_d_mat_partialpiv_lu(long N) {
    Eigen::MatrixXd d_mat = gen_dense_A(N);
    Eigen::VectorXd b = gen_b_vector(N);

    Eigen::VectorXd u = d_mat.partialPivLu().solve(b);
    return u;
  }
  ```,
  caption: [Rozwiązanie równania przy użyciu rozkładu LU z częściowym pivotingiem],
) <solve_par_lu_code>

#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/dense_partial_pivot_lu.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N z częściowym pivotem dla podanego N przy użyciu partialPivLU
  ],
) <par_piv_lu_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*], $451μ s$, $102m s$, $395m s$, $999m s$, $2138m s$, $4090m s$, $6753m s$, $11s$, $16s$, $22s$, $29s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu partialPivLU
  ],
) <par_piv_lu_data>

== Householder Full Pivot QR
Następną omawianą metodą jest dekompozycja QR z pełnym pivotingiem. Analizując wybrane dane z @full_piv_qr_data metoda rośnie w czasie $O(N^3)$ i potwierdza to teoretyczną złożoność algorytmu. W porównaniu z czasami rozkładu LU z pełnym pivotingiem (@full_piv_lu_data) jest około $1.4$ razy wolniejsza. W naszym przypadku rozwiązywania pojedyńczego rówania, rozkład QR nie oferuje wystarczających korzyści w odniesieniu do poniesionych kosztów.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd solve_d_mat_fullpiv_qr(long N) {
    Eigen::MatrixXd d_mat = gen_dense_A(N);
    Eigen::VectorXd b = gen_b_vector(N);

    Eigen::VectorXd u = d_mat.fullPivHouseholderQr().solve(b);
    return u;
  }

  ```,
  caption: [Rozwiązanie równania przy użyciu rozkładu QR z pełnym pivotem dla podanego N],
) <solve_full_qr_code>

#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/dense_full_pivot_qr.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N przy użyciu fullPivHouseholderQr
  ],
) <full_piv_qr_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*], $363μ s$, $512m s$, $5169m s$, $17s$, $39s$, $74s$, $122s$, $191s$, $282s$, $399s$, $499s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu fullPivHouseholderQr
  ],
) <full_piv_qr_data>


== Householder Partial Pivot QR
Ostatnia omawianą metoda dla macierzy gęstych jest dekompozycja QR z częściowym pivotingiem. Analizując czasy z @par_piv_qr_data metoda rośnie w czasie $O(N^3)$. Prównując ją z jej odpowiednikiem, czyli rozkładem LU z częściowym pivotingiem (@par_piv_lu_data) jest ona około $1.8$ razy wolniejsza. Jak w przypadku dekompozycji z pełnym pivotingiem nie oferuje ona wystarczających zysków w stosunku odpowiadającej jej metodzie z rozkładem LU.
#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd solve_d_mat_partialpiv_qr(long N) {
    Eigen::MatrixXd d_mat = gen_dense_A(N);
    Eigen::VectorXd b = gen_b_vector(N);

    Eigen::VectorXd u = d_mat.householderQr().solve(b);
    return u;
  }
  ```,
  caption: [Rozwiązanie równania przy użyciu rozkładu QR z częściowym pivotem dla podanego N],
) <solve_full_qr_code>

#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/dense_partial_qr.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N przy użyciu householderQr
  ],
) <par_piv_lu_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*], $208μ s$, $87m s$, $506 m s$, $1663m s$, $3869m s$, $7588m s$, $11s$, $19s$, $26s$, $43s$, $53s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu partialPivQR
  ],
) <par_piv_qr_data>

= Rozwiązania metodami dla macierzy rzadkich
Biblioteka Eigen oferuje również moduł służący do rozwiązywania układów równań przy użyciu macierzy rzadkich. Nasza macierz jest trójdagonalna z 2 elementami poza nimi. Dla dużych N jest ona w większości wypełniona zerami, co może zostać wykorzystane poprzez przedstawienie jej w reprezentacji rzadkiej. Wybrane zostały metody:
- Sparse LU
- Sparse QR

Obie metody uźywały wspólnej funkcji ```cpp gen_sparse_A()``` (@gen_sparse_mat_code)  do generowania rzadkiej macierzy. Metody uzywają gęstego wektora generowanego przez @gen_b_vec_code, gdyż jest on całkowicie zapełniony.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
   Eigen::SparseMatrix<double> gen_sparse_A(long N) {
    float64_t h = 2 / ((float64_t)N - 1.0);

    Eigen::SparseMatrix<double> sp_mat = Eigen::SparseMatrix<double>(N, N);
    std::vector<Eigen::Triplet<double>> tripletList =
        std::vector<Eigen::Triplet<double>>(N);

    for (long x = 1; x < N - 1; x++) {
      tripletList.push_back(Eigen::Triplet<double>(x, x, -2.0 / (h * h)));
      tripletList.push_back(Eigen::Triplet<double>(x, x + 1, 1.0 / (h * h)));
      tripletList.push_back(Eigen::Triplet<double>(x, x - 1, 1.0 / (h * h)));
    }

    // Diagonal fix
    tripletList.push_back(Eigen::Triplet<double>(0, 1, 1.0 / (h * h)));
    tripletList.push_back(Eigen::Triplet<double>(N - 1, N - 2, 1.0 / (h * h)));
    tripletList.push_back(Eigen::Triplet<double>(0, 0, -2.0 / (h * h)));
    tripletList.push_back(Eigen::Triplet<double>(N - 1, N - 1, -2.0 / (h * h)));

    // Corners
    tripletList.push_back(Eigen::Triplet<double>(N - 1, 0, 1.0 / (h * h)));
    tripletList.push_back(Eigen::Triplet<double>(0, N - 1, 1.0 / (h * h)));

    sp_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    sp_mat.makeCompressed();

    return sp_mat;
  }

  ```,
  caption: [Generowanie rzadkiej macierzy
  ],
) <gen_sparse_mat_code>


== Sparse LU
Pierwszą z opisywanych metod, które próbują wykorzystać wcześniej znaną strukture macierzy jest rozkład LU działający tylko na nie zerowych elementach macierzy. Dla $N < 10000$ rozkład skaluje się w przybliżeniu liniowo $O(N)$ po analizie czasów z @sparse_lu_data. Dokładniejsza analiza dla $N>10000$ zostanie przeprowadzona wspólnie z resztą względnie efektywnych metod. Dla $N<10000$ mimo uśredniania czasów dalej duży wpływ na efektywnośc metody posiada system operacyjny i przełączanie zadań. Wstępnie porównując czasy dla obecnie najbardziej efektywnej metody, czyli rozkładu LU z częściowym pivotingiem dla macierzy gęstej (@par_piv_lu_data), wybranie metody rzadkiej powoduje zwiększenie wydajności dla macierzy $N=9010$ o około 2200 razy.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd solve_sp_mat_lu(long N) {
    auto sp_mat = gen_sparse_A(N);
    auto b = gen_b_vector(N);

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
        sparse_lu_solver;
    sparse_lu_solver.analyzePattern(sp_mat);
    sparse_lu_solver.factorize(sp_mat);

    auto u = sparse_lu_solver.solve(b);
    return u;
  }
  ```,
  caption: [Rozwiązanie równania przy użyciu rozkładu LU wykorzystującego rzadką postac macierzy],
) <solve_sparse_lu_code>

#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/sparse_lu.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N przy użyciu sparse LU
  ],
) <sparse_lu_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*],
    $292μ s$,
    $3925μ s$,
    $3484μ s$,
    $4324μ s$,
    $5113μ s$,
    $6554μ s$,
    $7048μ s$,
    $7856μ s$,
    $8193μ s$,
    $9594μ s$,
    $9959μ s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu sparse LU
  ],
) <sparse_lu_data>

== Sparse QR
Drugą z oferowanych przez bibliotekę Eigen metod dekompozycji macierzy rzadkich jest zmodyfikowany rozkład QR. Nasza macierz mimo rzadkiej reprezentacji źle wpływa na rzadki rozkład QR i powoduje, że czasy rozwiązań z @sparse_qr_data są gorsze niż dla gęstych rozkładów QR i LU z częściowym pivotingiem. W porównaniu z gęstym QR z częściowym pivotingiem (@par_piv_qr_data) jest on około 1.5 razy wolniejsza. Jest to dobry przykład, że mimo rzadkiej struktury nie zawsze rzadki algorytm będzie lepszym rozwiązaniem dla czasu rozwiązania.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd solve_sp_mat_lu(long N) {
    auto sp_mat = gen_sparse_A(N);
    auto b = gen_b_vector(N);

    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
        sparse_lu_solver;
    sparse_lu_solver.analyzePattern(sp_mat);
    sparse_lu_solver.factorize(sp_mat);

    auto u = sparse_lu_solver.solve(b);
    return u;
  }
  ```,
  caption: [Rozwiązanie równania przy użyciu rozkładu QR wykorzystującego rzadką postac macierzy],
) <solve_sparse_qr_code>

#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/sparse_qr.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N przy użyciu sparse QR
  ],
) <sparse_qr_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*], $173μ s$, $112m s$, $722m s$, $2290m s$, $5425m s$, $10s$, $18s$, $29s$, $43s$, $61s$, $75s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu sparse QR
  ],
) <sparse_qr_data>

= Rozwiązanie wykorzystujące strukture macierzy
== Obserwacje

Analizując postać macierzy $A$ można zauważyć, że jest ona trójdiagonalna z dodatkowymi elementami w pozycjach $(1, N)$ oraz $(N, 1)$. Korzystając z obserwacji można zapisać macierz $A$ jako $A = B + u v^TT$, wówczas:
#grid(
  columns: (2fr, 0.5fr, 0.5fr, 1fr),
  align(center)[
    $
      B^(N times N) = mat(
        2*#diag, #other, , , ;
        #other, #diag, #other, , ;
        , #other, #diag, #other, , ;
        , , , dots.down, ;
        , , , #other, #diag - (1/(2h^2));
      )
    $
  ],
  align(center)[
    $ u = mat(2/h^2; 0; dots.v; 0; 1/h^2) $
  ],

  align(center)[
    $ v = mat(1; 0; dots.v; 0; 1/2) $
  ],

  align(center)[
    $
      u v^TT = mat(
        2/h^2, 0, ..., 0, #other ;
        0, 0, ..., 0, 0;
        dots.v, , dots.down, , dots.v;
        0, 0, ..., 0, 0;
        #other, 0, ..., 0, 1/(2h^2);
      )
    $
  ],
)

Pozwoli to rozwiązać układ równań przy użyciu wzoru Shermana-Morrisona:
$ (B + u v^TT)^(-1) = B^(-1) - (B^-1 u v^TT B^(-1)) / (1 + v^TT B^(-1) u) $
Dzięki zastosowaniu tego wzoru i postaci macierzy będziemy mogli rozwiązać układ równań $A u=b$ w czasie $O(N)$. Nasza macież $B$ jest trójdiagonalna co pozwala rozwiązywać układy $y = B^(-1)b$ oraz $z=B^(-1)u$ w czasie $O(N)$ poprzez użycie algorytmu Thomasa.

== Implementacja algorytmu Thomasa oraz rozwiązania wzorem Shermana-Morrisona

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  class ThomasSolver {
    Eigen::VectorXd y;
    Eigen::VectorXd z;
    Eigen::VectorXd c;
    /// Creates Thomas Solver for three diag matrix
    /// x upper diag 1 - N-1 (a_N is ignored)
    /// y diag 1 - N
    /// z lower diag 2 - N (c_1 is ignored)
  public:
    // Based on: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    ThomasSolver(Eigen::VectorXd &&x, Eigen::VectorXd &&y, Eigen::VectorXd &&z)
        : y(std::move(y)), z(std::move(z)) {

      long N = this->y.size();
      this->c = Eigen::VectorXd(N);

      this->c(0) = x(0) / this->y(0);
      for (long i = 1; i < this->y.size(); i++) {
        this->c(i) = x(i) / (this->y(i) - this->z(i) * this->c(i - 1));
      }
    }

    Eigen::VectorXd solve(const Eigen::VectorXd &b) const {
      auto N = this->y.size();
      assert(b.size() == N);
      Eigen::VectorXd d = Eigen::VectorXd::Zero(N);

      d(0) = b(0) / y(0);
      // Calculate coefficients in forwoard sweep
      for (long i = 1; i < this->y.size(); i++) {
        d(i) = (b(i) - this->z(i) * d(i - 1)) /
               (this->y(i) - this->z(i) * this->c(i - 1));
      }

      auto x = Eigen::VectorXd(N);
      // Backsubstitute for x_n
      x(this->y.size() - 1) = d(N - 1);
      for (long i = N - 2; i >= 0; i--) {
        x(i) = d(i) - c(i) * x(i + 1);
      }

      return x;
    }
  };
  ```,
  caption: [Własna implementacja algorytmu Thomasa],
) <thomas_solver_code>

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  // Create specific solver for our B matrix (A = B + uuT)
  // tr - top right corner
  // bl - bottom left corner
  ThomasSolver create_solver_for_B(long N, double tr, double bl) {
    Eigen::VectorXd a = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
    Eigen::VectorXd c = Eigen::VectorXd::Zero(N);

    float64_t h = 2 / ((float64_t)N - 1.0);

    for (long n = 0; n < N; n++) {
      a(n) = 1.0 / (h * h);
      b(n) = -2.0 / (h * h);
      c(n) = 1.0 / (h * h);
    }
    // Fix B_(1,1) and B_(N, N)
    float64_t gamma = -b(0);
    b(0) -= gamma;
    b(N - 1) -= (tr * bl) / gamma;

    return ThomasSolver(std::move(a), std::move(b), std::move(c));
  }
  ```,
  caption: [Solver dla macierzy B],
) <thomas_b_solver_code>

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd solve_mat_sherman_morrison(long N) {
    // Generate B matrix as threediagonal
    // Use thomas algorithm to solve intermidied equations for sherman morrison
    // solve for u using Sherman Morrison

    // A = B + uvT
    float64_t h = 2.0 / (N - 1);
    float64_t gamma = 2.0 / (h * h);

    // top left corner val
    float64_t tr = 1.0 / (h * h);
    // bototm lef corner val
    float64_t bl = 1.0 / (h * h);

    Eigen::VectorXd u = Eigen::VectorXd::Zero(N);
    u(0) = gamma;
    u(N - 1) = 1.0 / (h * h);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(N);
    v(0) = 1;
    v(N - 1) = (bl * tr) / gamma;

    auto B = create_solver_for_B(N, tr, bl);

    auto b = gen_b_vector(N);
    // y = B^-1 b => By = b
    auto y = B.solve(b);
    // z = B^-1 u = Bz = u
    auto z = B.solve(u);

    float64_t p = (v.transpose() * z);
    Eigen::VectorXd x = y - (z * (v.transpose() * y)) / (1.0 + p);

    return x;
  }
  ```,
  caption: [Rozwiązanie rówania przy użyciu wzoru Shermana-Morrisona],
) <sherman_solver_code>

== Analiza rozwiazania
Wykorzystanie rozwiązania wykorzystującego specyfikę naszej macierzy jesteśmy w stanie znajdować rozwiązanie układu w czasie $O(N)$ (@sherman_data). W prównaniu dla rzadkiego rozkładu LU (@sparse_lu_data) czasy rozwiązań są około 20 razy szybsze. Dla macierzy $N=9670$ wykorzystanie alorytmu Thomasa i wzoru Shermana-Morrisona powoduje wykonanie obliczeń w czasie mniejszym niż $1m s$ co dla wcześniej omawianych metod mogło by zostać uznane jako wariacje czasów wykonania miedzy próbami.
#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/sherman_morrison.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N przy użyciu wzoru Shermana-Morrisona i  algorytmu Thomasa
  ],
) <sherman_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*], $3μ s$, $31μ s$, $76μ s$, $128μ s$, $181μ s$, $240μ s$, $296μ s$, $365μ s$, $448μ s$, $498μ s$, $557μ s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu wzoru Shermana-Morrisona i  algorytmu Thomasa
  ],
) <sherman_data>

= Prównanie wszystkich metod
#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/all_combined_10k.jpg", width: 85%),
  caption: [Prównanie czasów rozwiązania $A u=b$ dla różnych N dla wszystkich omawianych metod
  ],
) <10k_methods_chart>
Na @10k_methods_chart widać zgrupowane metody, które układają się w cztery rozróżnialne grupy. Zaczynając od najwolniejszych są to metody używające rozkładów gęstych z pełnym pivotingiem. Są one o około 6 rzędów wielkości wolniejsze niż nasza najardziej efektywna metoda. Kolejną grupą są metody używające częściowego pivotingu oraz rozkład QR dla macierzy rzadkich. Są one o rząd wielkości szybsze niż metody używające pełengo pivotingu. Metody z tych dwóch grup chrakteryzują się czasem obliczeń $O(N^3)$ i róznią się jedynie stałym współczynnikiem niezależnym od wielkości macierzy. Pozostałe dwie metody cechują się czasem obliczeń $O(N)$ i zostaną omówione dogłebniej w kolejnej sekcji, gdyż one jedyne są rozsądnym wyborem do szukania rozwiązań dla macierzy $N>10000$.

== Porównanie Sparse LU oraz wzoru Shermana-Morrisona z algorytmem Thomasa
#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/big_mat_100k.jpg", width: 90%),
  caption: [Prównanie czasów rozwiązania $A u=b$ dla dużych N
  ],
) <100k_methods_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center + horizon,
    [*N*], $1099$, $10090$, $20080$, $30070$, $40060$, $50050$, $60040$, $70030$, $80020$, $90010$, $99001$,

    [*Sherman*],
    $36μ s$,
    $642μ s$,
    $1426μ s$,
    $2158μ s$,
    $2446μ s$,
    $3467μ s$,
    $4169μ s$,
    $4889μ s$,
    $5609μ s$,
    $6376μ s$,
    $7452μ s$,

    [*Sparse LU*], $754μ s$, $7630μ s$, $15m s$, $23m s$, $30m s$, $37m s$, $44m s$, $47m s$, $59m s$, $63m s$, $68m s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla dużych N efektywnymi metodami
  ],
) <100k_data>

Analizując dane z @100k_methods_chart oraz @100k_data rozwiązywanie równania przy użyciu najbardziej optymalnej metody daje nam około 10-krotnie szybsze wykonanie obliczeń. Rzadki rozkład LU dla macierzy $N=99001$ posiada porównywalny rząd czasu wykonania jak gęste LU z cześciowym pivotingiem (@par_piv_lu_data) lecz dla macierzy $N=1090$, więc w porównywalym czasie jesteśmy w stanie policzyć macierz 98 razy większą. Ekstrapolując czas wykonania przy założeniu, że gęste LU z częściowym pivotem posiada złożoność obliczeniowa $O(N^3)$ policzenie macierzy $N=96700$ zajeło by $29000s tilde.eq 8h$, pomijając fakt pamięci potrzebnej do jej przechowania (około 70GB).

== Analiza wyników
#grid(
  columns: (1fr, 1fr),
  align(center)[
    #figure(
      kind: "chart",
      supplement: [Wykres],
      image("./figures/sol_chart.jpg", width: 90%),
      caption: [Wykres zależności $u_n$ i wartości dla \ $N=1000$
      ],
    ) <solution_chart>
  ],
  align(center)[
    #figure(
      kind: "chart",
      supplement: [Wykres],
      image("./figures/sol_chart_with_estimation.jpg", width: 90%),
      caption: [Wykres przybliżenia wartości $u_n$ dla \ $N=1000$ przy użyciu $-cos(6.27x)/39.5 +0.0248$
      ],
    ) <solution2_chart>
  ],
)

@solution_chart pokazuje rozkład wartości wektora $arrow(u)$ unormowane na przedziale $(0, 2)$. Wartości te układują się w funkcje przypominająca $-cos(x)$. Na @solution2_chart została nałożona przybliżenie stworzone przu użyciu funkcji $-cos(6.27x)/39.5 +0.0248$. Jest to oczywiście przybliżenie, lecz można o użyć w metodzie iteracyjnej jako wektor początkowy co zagwarantuje nam szybką zbierzność i przyśpieszy czas rozwiązania równania. Głebsza analiza wyników i próba znalezienia dokładniejszego przybliżenia mogła by się okazać dokładnym rozwiązaniem naszego układu równań.

= Podsumowanie
Po przeanalizowaniu wszystkich metod można zauważyć, że dobór odpowiedniej metody jest kluczową rzeczą przy rozwiązywaniu układów równań. Znajomość struktury macierzy pozwala nam zejść z czasem wykonania o 6 rzędów wielkości dla średnich macierzy. Dla dużych macierzy różnice te są jeszcze bardziej widoczne i dobór i implementacja odpowiedniej metody jest bardziej efektywna czasowo niż siłowe rozwiązanie metodami gęstymi. (Opis czasów dla macierzy 100000x100000). Analiza wynikowego wektora może umożliwić nam znalezienie funkcji przybliżającej rozwiązanie. Może być użyteczne w przypadku macierzy, których unormowany rozkład wyników będzie posiadał podobną strukture.
