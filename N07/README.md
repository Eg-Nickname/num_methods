# Jak uruchomić

Program można uruchomić przy użyciu poleceń po sklonowaniu całego repozytorium (https://github.com/Eg-Nickname/num_methods), lub po pobraniu fragmentu kodu i edycji ściezki do biliteki eigen w pliku Makefile

```bash
make run # generuje dane do wykresów oraz liczy rozwiązania układu równań
python ./scripts/plot.py # generowanie wykresów
typst compile N07_raport.typ # wymagany typst. Gotowy pdf w repozytorium
make compile && ./target/exe/main.x #  analogicznie do run
```
