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
Macierz $A$ można zapisac jako $A = B + u u^TT$, wówczas:
#grid(
  columns: (1fr, 0.5fr, 1fr),
  align(center)[
    $
      B^(N times N) = mat(
        (-3)/h^2, #other, , , ;
        #other, #diag, #other, , ;
        , #other, #diag, #other, , ;
        , , , dots.down, ;
        , , , #other, (-3)/h^2
      )
    $
  ],
  align(center)[
    $ u = mat(1/h; 0; dots.v; 0; 1/h) $
  ],
  align(center)[
    $
      u u^TT = mat(
        #other, 0, ..., 0, #other ;
        0, 0, ..., 0, 0;
        dots.v, , dots.down, , dots.v;
        0, 0, ..., 0, 0;
        #other, 0, ..., 0, #other ;
      )
    $
  ],
)

Pozwoli to rozwiązać układ równań przy użyciu wzoru Shermana-Morrisona ($v =u$):
$ (A + u v^TT)^(-1) = A^(-1) - (A^-1 u v^TT A^(-1)) / (1 + v^TT A^(-1) u) $
