# Jak uruchomić

Program można uruchomić przy użyciu poleceń po sklonowaniu całego repozytorium (https://github.com/Eg-Nickname/num_methods) lub bezpośrednio z kodu.

```bash
make run # generuje dane do wykresów oraz znajduje najbardziej optymalną interpolacje
python ./scripts/prrety_plot_gen.py # generowanie wykresów
typst compile N04_raport.typ # wymagany typst. Gotowy pdf w repozytorium
make compile && ./target/exe/main.x -p # analogicznie do run
```
