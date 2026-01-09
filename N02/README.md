# Jak uruchomić

Program można uruchomić przy użyciu poleceń po sklonowaniu całego repozytorium (https://github.com/Eg-Nickname/num_methods), lub po pobraniu fragmentu kodu i edycji ściezki do biliteki eigen w pliku Makefile

```bash
make run_methods_10k # generuje dane do wykresów dla macierzy 100 - 10000
make big_mat # generuje dane do wykresów ``dla macierzy 100 - 100 000
make gen_figures # create plots with generated data
typst compile N02cd_raport.typ # wymagany typst. Gotowy pdf w repozytorium
make compile && ./target/exe/main.x -p # Gen data for cmp plot
```
