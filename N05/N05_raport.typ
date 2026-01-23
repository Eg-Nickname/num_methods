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
  Metody Numeryczne N05
]
#align(center)[Jakub Kurek]

= Wstęp
Celem zadania jest numerycznie rozwiązanie całki:
$ I = integral_(-1)^1 exp(x^2)/sqrt(1-x^2) dif x $

Program liczący rozwiązania został napisany w C++23. Cały kod wykorzystywany do obliczeń oraz generowania wykresów znajduje się w repozytorium na #custom_link("https://github.com/Eg-Nickname/num_methods/tree/main/N05")[GitHub].

= Przekształcenie całki
Aby całka nie posiadała nieskończoności na granicach przedziału, możemy dokonać podstawienia:
$
      x & = sin(t) \
  dif x & = cos(t) dif t
$
Obliczamy nowe przedziały całkowania:
$
  x = -1 & => sin(t) = -1 arrow.r t = -pi/2 \
   x = 1 & => sin(t) = 1 arrow.r t = pi/2
$
Podstawiamy do naszej całki i upraszczamy:
$
  I = integral_(-pi/2)^(pi/2) exp(sin(t)^2)/sqrt(1-sin(x)^2) cos(t) dif t = integral_(-pi/2)^(pi/2) exp(sin(t)^2)/sqrt(cos(x)^2) cos(t) dif t = integral_(-pi/2)^(pi/2) exp(sin(t)^2) dif t
$

= Całkowanie numeryczne
Całki były liczone przy użyciu zamkniętych kwadratur Newtona-Cotensa z różnymi stopniami wielomianów interpolacyjnych, wraz z iteracyjnym zagęszczaniem przedziałów. Wiele przybliżeń całek na małych przedziałach daje lepsza dokładność niż jeden duży wielomina interpolacyjny.

== Zagęszczanie przedziałów
W każdej iteracji powstawało $2 dot (N-1)$ dodatkowych punktów, które były umieszczone pomiędzy punktami poprzedniej iteracji. Wartości obliczonej funkcji całkowej były przechowywane co pozwoliło na ponowne wykorzystanie raz obliczonej wartości w każdej późniejszej iteracji.

== Metody
Użyte metody:
- Metoda trapezów
$ integral_a^b tilde.eq (b-a)/2 dot (f_0 + f_1) $
- Metoda Simpsona
$ integral_a^b tilde.eq (b-a)/6 dot (f_0 + 4 f_1 + f_2) $
#pagebreak()
- Metoda 3/8
$ integral_a^b tilde.eq (b-a)/8 dot (f_0 + 3 f_1 + 3 f_2 + f_3) $

== Porównanie metod
#figure(
  kind: "table",
  supplement: [Tabela],
  table(
    columns: (auto, auto),
    inset: 8pt,
    align: horizon + center,
    table.header([*Metoda*], [*Ilość obliczeń funkcji*]),
    "trapezów", [$17$],
    "Simpsona", [$33$],
    "3/8", [$49$],
  ),

  caption: [Ilość obliczonych wartości funkcji dla różnych metod
  ],
) <biggest_evals>

Obliczenia były prowadzone do momentu spełnienia warunku stopu:
$ |I_n - I_(n-1)| < epsilon, quad "gdzie" epsilon = 10^(-6) $

Koniec obliczen następował dla:
$ I = 5.508429773886 $

#figure(
  kind: "plot",
  supplement: [Wykres],
  image("./figures/trapezoid_plot.jpg"),

  caption: [Całkowana funkcja z metody trapezów
  ],
) <trapezoid_plot>

= Podsumowanie
Najbardziej efektywna okazała się metoda trapezów, co pokazuje @trapezoid_plot funkcja ta dla 17 punktów bardzo dobrze przybliża całkowaną funkcję. Wykonała ona najmniejszą liczbę obliczeń funkcji całkowanej, co stanowi największy numeryczny koszt we wszystkich  metodach. Podstawienie $x=sin(t)$ pozwoliło nam na usunięcie nieskończoności na krancach przedziałów i zastosowanie zamkniętych kwadratur Newtona-Cotensa.
