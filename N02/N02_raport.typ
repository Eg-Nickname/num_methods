#import "@preview/physica:0.9.7": *

Rozwiązujemy układ równań $A u=b$

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
= Rozwiązania metodami gęstymi
== LUP
== Householder QR

= Rozwiązania metodami dla macierzy rzadkich
== Sparse LU
== Sparse QR

= Rozwiązanie wykorzystujące strukture macierzy
Macierz $A$ można zapisac jako $A = B + u v^TT$, wówczas:
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
$ (A + u v^TT)^(-1) = A^(-1) - (A^-1 u v^TT A^(-1)) / (1 + v^TT A^(-1) u) $
