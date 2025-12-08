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
  Metody Numeryczne N01
]
#align(center)[Jakub Kurek]

= Wstęp
Celem zadania jest sprawdzenie różnych metod rozwiązania układu równań $A u=b$

#let diag = $(-2) / h^2$
#let other = $1 / h^2$
#grid(
  columns: (1fr, 1fr, 1fr),
  align(center)[
    $
      A^(N times N) = mat(
        1, , , , ;
        #other, #diag, #other, , ;
        , #other, #diag, #other, , ;
        , , , dots.down, ;
        , , , , 1
      )
    $],
  align(center + horizon)[
    $ h=2/(N-1) $
  ],
  align(center + horizon)[
    $ b = mat(1; 0; dots.v; 0; 1) $
  ],
)
Program liczący rozwiązania został napisany w C++23 przy użyciu biblioteki #custom_link("https://libeigen.gitlab.io/")[Eigen] w wersji 5.0.1. Wykresy zostały stworzone w języku Python przy użyciu biblioteki matplotlib. Cały kod wykorzystywany do obliczeń oraz generowania wykresów znajduje się w repozytorium na #custom_link("https://github.com/Eg-Nickname/num_methods/tree/main/N01")[GitHub].

== Konfiguracja sprzętowa
- Procesor: Intel i7-8650U (8) @ 4.200GHz
- Pamięć RAM: 32GB DDR4 2133MHz
- System operacyjny: Arch Linux 6.17.5-arch1-1

== Metodyka testowania
Wszystkie testy były wykonywane na odciążonym procesorze bezpośrednio z powłoki. Dla obliczeń trwających poniżej 20ms testy były wykonywane 10 krotnie a czas był uśredniany. Pomiary czasu były wykonywane w mikrosekundach w celu lepszego zobrazowania czasu obliczeń dla małych macierzy oraz optymalnych metod. Po wprowadzeniu optymalizacje wszystkie dalsze testy były kompilowane przy użyciu kompilatora clang++ w wersji 21.1.4 z flagami ```-Wall -Wextra -std=c++23 -pedantic -I ./../eigen-5.0.1/ -O3 -march=native -DNDEBUG```.

= Wstępne optymalizacje
Najważniejszą optymalizacją poprawiające wyniki obliczeń było zastosowanie odpowiednich flag kompilatora:

- ```-O3``` - poziom optymalizacji kodu przez kompilator.  @o3cmp pokazuje róznice czasu wykonania zwiększoną około 10-krotnie.
#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/O3_comparison_100k.jpg", width: 70%),
  caption: [Porównanie algorytmu Thomasa z optymalizacją -O3 i bez
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
  - ```-DNDEBUG``` - wyłączenie assercji. Eigen dla macierzy o dynamicznych rozmiarach sprawdza poprawność ich rozmiaru przy użyciu assercji. Nie ma wielkiego wpływu na pojedyncze działanie programu dla dużych macierzy, lecz może wpływać kiedy wykonujemy wiele działań na wektorach lub macierzach. Zalecane przez dokumentacje.
]

\
= Rozwiązania metodami gęstymi
Biblioteka Eigen oferuje moduł służący do rozwiązywania układów równań przy użyciu macierzy gęstych. Wybrane zostały metody:
- Full Pivot LU
- Partial Pivot LU
- Householder Full Pivot QR
- Householder Partial Pivot QR

Wszystkie metody uzywały wspólnej funkcji ```cpp gen_desnse_A()``` (@gen_dense_mat_code)  do generowania gęstej macierzy oraz funkcji ```cpp gen_b_vector()``` (@gen_b_vec_code) do tworzenia wektora rowiązania.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::MatrixXd gen_dense_A(long N) {
    // matrix must be greater than 1
    // assert(N != 0);
    float64_t h = 2 / ((float64_t)N - 1.0);
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(N, N);

    for (long x = 1; x < N - 1; x++) {
      mat(x, x - 1) = 1.0 / (h * h);
      mat(x, x + 1) = 1.0 / (h * h);
      mat(x, x) = -2.0 / (h * h);
    }

    mat(0, 0) = 1.0;
    mat(N - 1, N - 1) = 1.0;

    return mat;
  }
  ```,
  caption: [Generowanie gęstej macierzy pełnej
  ],
) <gen_dense_mat_code>

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::VectorXd gen_b_vector(long N) {
    Eigen::VectorXd b = Eigen::VectorXd(N);
    b.setZero();
    b(0) = 1;
    b(N - 1) = 1;
    return b;
  }
  ```,
  caption: [Generowanie wektora rozwiązania
  ],
) <gen_b_vec_code>


== Full Pivot LU
Pierwszą z omawianych metod jest rozwiązanie układu równań przy użyciu dekompozycji LU z pełnym pivotingiem. Analizując wybrane dane z @full_piv_lu_data  metoda rośnie w czasie $O(N^3)$ co jest zgodne z teoretyczna złozonością algorytmu. Jest to bazowy wynik, do którego kolejne metody będą się odnosić.
// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   Eigen::VectorXd solve_d_mat_fullpiv_lu(long N) {
//     Eigen::MatrixXd d_mat = gen_dense_A(N);
//     Eigen::VectorXd b = gen_b_vector(N);

//     Eigen::VectorXd u = d_mat.fullPivLu().solve(b);
//     return u;
//   }
//   ```,
//   caption: [Rozwiązanie równania przy użyciu rozkładu LU z pełnym pivotem dla podanego N],
// ) <solve_full_lu_code>

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

    [*Czas*], $164μ s$, $272m s$, $2998m s$, $10s$, $25s$, $50s$, $91s$, $144s$, $228s$, $336s$, $436s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu fullPivLU
  ],
) <full_piv_lu_data>

== Partial Pivot LU
Kolejną z omawianych metod jest rozwiązanie układu równań przy użyciu dekompozycji LU z częściowym pivotingiem. Kosztem precyzji numerycznej zyskujemy krótszy czas wykonania. Analizując wybrane dane z @par_piv_lu_data  metoda rośnie w czasie $O(N^3)$ , lecz w porównaniu z czasami dla pełnego pivotingu (@full_piv_lu_data) czas rozwiązywania układu równań jest ponad 14-krotnie mniejszy.

// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   Eigen::VectorXd solve_d_mat_partialpiv_lu(long N) {
//     Eigen::MatrixXd d_mat = gen_dense_A(N);
//     Eigen::VectorXd b = gen_b_vector(N);

//     Eigen::VectorXd u = d_mat.partialPivLu().solve(b);
//     return u;
//   }
//   ```,
//   caption: [Rozwiązanie równania przy użyciu rozkładu LU z częściowym pivotingiem],
// ) <solve_par_lu_code>

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

    [*Czas*], $108μ s$, $100m s$, $591m s$, $1206m s$, $2530m s$, $4762m s$, $7967m s$, $12s$, $17s$, $25s$, $30s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu partialPivLU
  ],
) <par_piv_lu_data>

== Householder Full Pivot QR
Następną omawianą metodą jest dekompozycja QR z pełnym pivotingiem. Analizując wybrane dane z @full_piv_qr_data metoda rośnie w czasie $O(N^3)$ i potwierdza to teoretyczną złożoność algorytmu. W porównaniu z czasami rozkładu LU z pełnym pivotingiem (@full_piv_lu_data) jest około $1.4$ razy wolniejsza. W naszym przypadku rozwiązywania pojedynczego rówania, rozkład QR nie oferuje wystarczających korzyści w odniesieniu do poniesionych kosztów.

// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   Eigen::VectorXd solve_d_mat_fullpiv_qr(long N) {
//     Eigen::MatrixXd d_mat = gen_dense_A(N);
//     Eigen::VectorXd b = gen_b_vector(N);

//     Eigen::VectorXd u = d_mat.fullPivHouseholderQr().solve(b);
//     return u;
//   }

//   ```,
//   caption: [Rozwiązanie równania przy użyciu rozkładu QR z pełnym pivotem dla podanego N],
// ) <solve_full_qr_code>

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

    // TODO update
    [*Czas*], $393μ s$, $807m s$, $7225m s$, $22s$, $48s$, $91s$, $153s$, $242s$, $360s$, $508s$, $621s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu fullPivHouseholderQr
  ],
) <full_piv_qr_data>


== Householder Partial Pivot QR
Ostatnia omawianą metoda dla macierzy gęstych jest dekompozycja QR z częściowym pivotingiem. Analizując czasy z @par_piv_qr_data metoda rośnie w czasie $O(N^3)$. Prównując ją z jej odpowiednikiem, czyli rozkładem LU z częściowym pivotingiem (@par_piv_lu_data) jest ona około $2$ razy wolniejsza. Jak w przypadku dekompozycji z pełnym pivotingiem nie oferuje ona wystarczających zysków w stosunku odpowiadającej jej metodzie z rozkładem LU.
// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   Eigen::VectorXd solve_d_mat_partialpiv_qr(long N) {
//     Eigen::MatrixXd d_mat = gen_dense_A(N);
//     Eigen::VectorXd b = gen_b_vector(N);

//     Eigen::VectorXd u = d_mat.householderQr().solve(b);
//     return u;
//   }
//   ```,
//   caption: [Rozwiązanie równania przy użyciu rozkładu QR z częściowym pivotem dla podanego N],
// ) <solve_full_qr_code>

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

    // TODO update
    [*Czas*], $195μ s$, $59m s$, $405m s$, $1349m s$, $3184m s$, $6780m s$, $12s$, $21s$, $34s$, $51s$, $64s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu partialPivQR
  ],
) <par_piv_qr_data>

== Porównanie metod gęstych
#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/dense_combined_10k.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla $2000 < N < 10000$ metody gęste
  ],
) <dense_chart>

= Rozwiązania metodami dla macierzy rzadkich
Biblioteka Eigen oferuje również moduł służący do rozwiązywania układów równań przy użyciu macierzy rzadkich. Nasza macierz jest trójdagonalna z 2 elementami poza nimi. Dla dużych N jest ona w większości wypełniona zerami, co może zostać wykorzystane poprzez przedstawienie jej w reprezentacji rzadkiej. Wybrane zostały metody:
- Sparse LU
- Sparse QR

Obie metody uźywały wspólnej funkcji ```cpp gen_sparse_A()``` (@gen_sparse_mat_code)  do generowania rzadkiej macierzy. Metody uzywają gęstego wektora generowanego przez @gen_b_vec_code, gdyż zmiana go na sparse vector nie oferuje zauważalnych zysków.

#figure(
  kind: "code",
  supplement: [Kod],
  ```cpp
  Eigen::SparseMatrix<double> gen_sparse_A(long N) {
    // matrix must be greater than 1
    // assert(N != 0);
    float64_t h = 2 / ((float64_t)N - 1.0);

    Eigen::SparseMatrix<double> sp_mat = Eigen::SparseMatrix<double>(N, N);
    std::vector<Eigen::Triplet<double>> tripletList =
        std::vector<Eigen::Triplet<double>>(N);

    for (int i = 1; i < N - 1; i++) {
      tripletList.push_back(Eigen::Triplet<double>(i, i, -2.0 / (h * h)));
      if (i < N - 1)
        tripletList.push_back(Eigen::Triplet<double>(i, i + 1, 1.0 / (h * h)));
      if (i > 0)
        tripletList.push_back(Eigen::Triplet<double>(i, i - 1, 1.0 / (h * h)));
    }
    sp_mat.setFromTriplets(tripletList.begin(), tripletList.end());

    sp_mat.coeffRef(0, 0) = 1;
    sp_mat.coeffRef(N - 1, N - 1) = 1;
    sp_mat.makeCompressed();
    return sp_mat;
  }

  ```,
  caption: [Generowanie rzadkiej macierzy
  ],
) <gen_sparse_mat_code>


== Sparse LU
Pierwszą z opisywanych metod, które próbują wykorzystać wcześniej znaną strukture macierzy jest rozkład LU działający tylko na nie zerowych elementach macierzy. Dla $N < 10000$ rozkład skaluje się w przybliżeniu liniowo $O(N)$ po analizie czasów z @sparse_lu_data. Dokładniejsza analiza dla $N>10000$ zostanie przeprowadzona wspólnie z resztą względnie efektywnych metod. Dla $N<10000$ mimo uśredniania czasów dalej niewielki wpływ na efektywnośc metody posiada system operacyjny i przełączanie zadań (ostatnie pomiary nie są znacząco rosnące). Wstępnie porównując czasy dla obecnie najbardziej efektywnej metody, czyli rozkładu LU z częściowym pivotingiem dla macierzy gęstej (@par_piv_lu_data), wybranie metody rzadkiej powoduje zwiększenie wydajności dla macierzy $N=9010$ o około 5000 razy.

// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   Eigen::VectorXd solve_sp_mat_lu(long N) {
//     auto sp_mat = gen_sparse_A(N);
//     auto b = gen_b_vector(N);

//     Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
//         sparse_lu_solver;
//     sparse_lu_solver.analyzePattern(sp_mat);
//     sparse_lu_solver.factorize(sp_mat);

//     auto u = sparse_lu_solver.solve(b);
//     return u;
//   }
//   ```,
//   caption: [Rozwiązanie równania przy użyciu rozkładu LU wykorzystującego rzadką postac macierzy],
// ) <solve_sparse_lu_code>

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
    $44μ s$,
    $563μ s$,
    $1152μ s$,
    $1756μ s$,
    $2348μ s$,
    $2902μ s$,
    $3467μ s$,
    $4057μ s$,
    $4642μ s$,
    $5378μ s$,
    $5860μ s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu sparse LU
  ],
) <sparse_lu_data>

== Sparse QR
Drugą z oferowanych przez bibliotekę Eigen metod dekompozycji macierzy rzadkich jest zmodyfikowany rozkład QR. Analizując czasy z @sparse_lu_data zauważamy, że metoda jest około 5 razy wolniejsza niż odpowiadający jej rozkład LU. Z samych danych z tabeli ciężko stwierdzić złożonośc algorytmu, gdyż dla $N=1090$ i $N=2080$, dla 2-krotnego wzrostu, czas rośnie 2.5-krotnie. Natomiast dla macierzy $N=4060$ i $N=8020$ prawie 3-krotnie. Dokładniejsza analiza czasów wykonania i złożoności czasowej dla $N>10000$ zostanie przeprowadzona wspólnie z resztą względnie efektywnych metod.

// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   Eigen::VectorXd solve_sp_mat_lu(long N) {
//     auto sp_mat = gen_sparse_A(N);
//     auto b = gen_b_vector(N);

//     Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
//         sparse_lu_solver;
//     sparse_lu_solver.analyzePattern(sp_mat);
//     sparse_lu_solver.factorize(sp_mat);

//     auto u = sparse_lu_solver.solve(b);
//     return u;
//   }
//   ```,
//   caption: [Rozwiązanie równania przy użyciu rozkładu QR wykorzystującego rzadką postac macierzy],
// ) <solve_sparse_qr_code>

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

    // TODO update
    [*Czas*],
    $29μ s$,
    $502μ s$,
    $1385μ s$,
    $3507μ s$,
    $5687μ s$,
    $8354μ s$,
    $11m s$,
    $14m s$,
    $16m s$,
    $21m s$,
    $24m s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu sparse QR
  ],
) <sparse_qr_data>

= Rozwiązanie wykorzystujące strukture macierzy
== Obserwacje
Nasza macierz $A$ posiada struktóre trójdagonalną co pozwala nam skorzystać z algorytmu Thomasa (brak zer na głównej diagonali), który sprowadza liczbę operacji do minimum i pozwala rozwiązać układ rówań w czasie $O(N)$.

== Implementacja algorytmu
// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   class ThomasSolver {
//     Eigen::VectorXd y;
//     Eigen::VectorXd z;
//     Eigen::VectorXd c;
//     /// Creates Thomas Solver for three diag matrix
//     /// x upper diag 1 - N-1 (a_N is ignored)
//     /// y diag 1 - N
//     /// z lower diag 2 - N (c_1 is ignored)
//   public:
//     // Based on: https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
//     ThomasSolver(Eigen::VectorXd &&x, Eigen::VectorXd &&y, Eigen::VectorXd &&z)
//         : y(std::move(y)), z(std::move(z)) {

//       long N = this->y.size();
//       this->c = Eigen::VectorXd(N);

//       this->c(0) = x(0) / this->y(0);
//       for (long i = 1; i < this->y.size(); i++) {
//         this->c(i) = x(i) / (this->y(i) - this->z(i) * this->c(i - 1));
//       }
//     }

//     Eigen::VectorXd solve(const Eigen::VectorXd &b) const {
//       auto N = this->y.size();
//       assert(b.size() == N);
//       Eigen::VectorXd d = Eigen::VectorXd::Zero(N);

//       d(0) = b(0) / y(0);
//       // Calculate coefficients in forwoard sweep
//       for (long i = 1; i < this->y.size(); i++) {
//         d(i) = (b(i) - this->z(i) * d(i - 1)) /
//                (this->y(i) - this->z(i) * this->c(i - 1));
//       }

//       auto x = Eigen::VectorXd(N);
//       // Backsubstitute for x_n
//       x(this->y.size() - 1) = d(N - 1);
//       for (long i = N - 2; i >= 0; i--) {
//         x(i) = d(i) - c(i) * x(i + 1);
//       }

//       return x;
//     }
//   };
//   ```,
//   caption: [Własna implementacja algorytmu Thomasa],
// ) <thomas_solver_code>

// #figure(
//   kind: "code",
//   supplement: [Kod],
//   ```cpp
//   ThomasSolver create_solver_for_A(long N) {
//     Eigen::VectorXd a = Eigen::VectorXd::Zero(N);
//     Eigen::VectorXd b = Eigen::VectorXd::Zero(N);
//     Eigen::VectorXd c = Eigen::VectorXd::Zero(N);

//     float64_t h = 2 / ((float64_t)N - 1.0);

//     for (long n = 0; n < N; n++) {
//       a(n) = 1.0 / (h * h);
//       b(n) = -2.0 / (h * h);
//       c(n) = 1.0 / (h * h);
//     }
//     a(0) = 0;
//     c(N - 1) = 0;
//     b(0) = 1;
//     b(N - 1) = 1;

//     return ThomasSolver(std::move(a), std::move(b), std::move(c));
//   }

//   ```,
//   caption: [Solver dla macierzy A],
// ) <thomas_a_solver_code>

== Analiza wyników algorytmu Thomasa
Wykorzystanie rozwiązania wykorzystującego specyfikę naszej macierzy jesteśmy w stanie znajdować rozwiązanie układu w czasie $O(N)$ (@thomas_data). W prównaniu dla rzadkiego rozkładu LU (@sparse_lu_data) czasy rozwiązań są około 20 razy krótsze. Dla macierzy $N=9670$ wykorzystanie alorytmu Thomasa. Mimo powtarzania i uśredniania czasów rozwązań, dla macierzy $4000<N<8000$ osylują wokół jednego czasu, co wskazuje na duży wpływ systemu peracyjnego na rozwiazania, w tak krótkim czasie.

#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/10k_sep/thomas_algorithm.jpg", width: 85%),
  caption: [Czas rozwiązania $A u=b$ dla różnych N przy algorytmu Thomasa
  ],
) <thomas_chart>

#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto, auto),
    align: center,
    [*N*], $100$, $1090$, $2080$, $3070$, $4060$, $5050$, $6040$, $7030$, $8020$, $9010$, $9670$,

    [*Czas*], $2μ s$, $29μ s$, $58μ s$, $113μ s$, $177μ s$, $173μ s$, $161μ s$, $176μ s$, $189μ s$, $220μ s$, $241μ s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla wybranych N przy użyciu wzoru Shermana-Morrisona i  algorytmu Thomasa
  ],
) <thomas_data>

= Prównanie wszystkich metod
#figure(
  kind: "chart",
  supplement: [Wykres],
  image("./figures/all_combined_10k.jpg", width: 85%),
  caption: [Prównanie czasów rozwiązania $A u=b$ dla różnych N dla wszystkich omawianych metod
  ],
) <10k_methods_chart>
Na @10k_methods_chart widać zgrupowane metody, które układają się w cztery rozróżnialne grupy. Zaczynając od najwolniejszych są to metody używające rozkładów gęstych z pełnym pivotingiem. Są one o około 6 rzędów wielkości wolniejsze niż algorytm Thomasa. Kolejną grupą są metody używające częściowego pivotingu, sa około 10-krotnie szybsze niż ich odpowiedniki z pełnym pivotingiem. Metody z tych dwóch grup chrakteryzują się czasem obliczeń $O(N^3)$ i róznią się jedynie stałym współczynnikiem niezależnym od wielkości macierzy. Nastepną grupą sa metody rzadkie. Najszybszy w osobnej grupie jest algorytm Thomasa. Metody z dwóch ostatnich grup zostaną omówione dogłebniej w kolejnej sekcji, gdyż one jedyne są rozsądnym wyborem do szukania rozwiązań dla macierzy $N>10000$.

== Porównanie metod efektywnych
Analizując dane z @100k_methods_chart dobrze widać, że nawet dla efektywnych metod czasy mogą się różnić o ponad 2-rzędy wielkosci (Alorytm Thomasa i Sparse QR). Dla rzadkiego rozkładu QR dopiero widać, iż nie rośnie on w czasie liniowym $O(N)$ tak jak rzadkie LU i algorytm Thomasa. Z danych @100k_data można zauwazyć, że złożoność czasowa rośnie w czasie wielomianowym $O(N^2)$. Dla pojedyńczych macierzy $N tilde.eq 100000$, wszystkie metody dają względnie rozsądne czasy rozwiązań, chociaż bezdyskusyjnie widać przewagę metod $O(N)$. Obie metody mimo liniowego czasu wykonania różnią się od siebie o rząd wielkości poprzez stały współczynik wykonywanyhch operacji.
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
    columns: (1.2fr, 0.85fr, 0.95fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr, 1fr),
    align: center + horizon,
    [*N*], $1099$, $10090$, $20080$, $30070$, $40060$, $50050$, $60040$, $70030$, $80020$, $90010$, $99001$,

    [*Thomas*],
    $20μ s$,
    $420μ s$,
    $1135μ s$,
    $1878μ s$,
    $2479μ s$,
    $3369μ s$,
    $4048μ s$,
    $4233μ s$,
    $5369μ s$,
    $5328μ s$,
    $6034μ s$,

    [*Sparse LU*], $977μ s$, $7310μ s$, $12m s$, $24m s$, $31m s$, $37m s$, $38m s$, $46m s$, $49m s$, $64m s$, $67m s$,

    [*Sparse QR*],
    $671μ s$,
    $36m s$,
    $119m s$,
    $261m s$,
    $472m s$,
    $767m s$,
    $1087m s$,
    $1540m s$,
    $2080m s$,
    $2986m s$,
    $3580m s$,
  ),

  caption: [Czasy rozwiązań $A u=b$ dla dużych N efektywnymi metodami
  ],
) <100k_data>

== Analiza wyników
#grid(
  columns: (1fr, 1fr),
  align(center)[
    $
      A^(N times N) = mat(
        1, , , , ;
        #other, #diag, #other, , ;
        , #other, #diag, #other, , ;
        , , , dots.down, ;
        , , , , 1
      )
    $],
  align(center + horizon)[
    $ b = mat(1; 0; dots.v; 0; 1) $
  ],
)
Analizując dokładniej strukturę naszej macierzy oraz wektora rozwiązania należy zwrócić uwagę na to, że $u_1$ oraz $u_n$ muszą być równe $1$. Dla $1<n<N$ nasza macierz jest symetryczna i każdy wiersz musi się zerować. Elementy na lewo i lewo od diagonali po dodaniu towrzą element diagonalny z przeciwnym znakiem. Pierwsze równanie oraz ostatnie. Dyktują nam współczyniki przy $u$, które propagują się na cały wektor rozwiązania, który zapełnia się $1$. @solution_chart i @solution2_chart pokazują, że oczekiwania zgadzają się z rzeczywistością. Znając analizę tego wyniku, można zauwazyć, że dla naszego $b$ rozwiązanie $u$ zawsze będzie mieć postać $u = mat(1, dots.h.c, 1)^TT$. Co neguje jakąkolwiek potrzebę przeprowadzania rozkładów macierzy $A$.

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
      caption: [Wykres przybliżenia wartości $u_n$ dla \ $N=1000$ przy użyciu $u_n = 1$
      ],
    ) <solution2_chart>
  ],
)

= Podsumowanie
Po przeanalizowaniu wszystkich metod można zauważyć, że dobór odpowiedniej metody jest kluczową rzeczą przy rozwiązywaniu układów równań. Znajomość struktury macierzy pozwala nam skrócić czas wykonania o 6 rzędów wielkości dla średnich macierzy. Dla dużych macierzy różnice te są jeszcze bardziej widoczne i dobór i implementacja odpowiedniej metody jest bardziej efektywna czasowo niż siłowe rozwiązanie metodami gęstymi. Wstepna analiza rozwiązania pozwala nam, również zaoszczędzić czas przy liczeniu macierzy, lecz zyski czasowe będą nie wystarczające, gdyż wykorzystanie samejj postaci trójdiagonalnej będzie wystarczająca optymalizacją dla większości przypadków.
