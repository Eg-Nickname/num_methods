import os
import pandas as pd
import matplotlib.pyplot as plt


df_zeroed = pd.read_csv("./data/zeroed_u_n.csv", delimiter=";")
df_oned = pd.read_csv("./data/oned_u_n.csv", delimiter=";")
df_random = pd.read_csv("./data/random_u_n.csv", delimiter=";")
df_linear = pd.read_csv("./data/linear_u_n.csv", delimiter=";")

plt.figure(figsize=(10, 6))
plt.plot(df_oned["u_n"], label="$u_n$ jeden", color="#6929c4")
plt.plot(df_linear["u_n"], label="$u_n$ losowe", color="#3ddbd9")
plt.plot(df_zeroed["u_n"], label="$u_n$ zero", color="#0072c3")
plt.plot(
    df_random["u_n"],
    label="$u_n$ liniowe",
    color="#005d5d",
)


plt.xlabel("N")
plt.ylabel("$u_n$")


ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)


plt.xlim(0.0, 1000.05)


plt.grid(True, linestyle="--", linewidth=0.5)
plt.legend()

plt.subplots_adjust(
    left=0.10,  # lewy margines (0–1)
    right=0.90,  # prawy margines
    top=0.9,  # górny margines
    bottom=0.1,  # dolny margines
)

plt.grid(True)
plt.legend()
plt.savefig("./figures/solution_plot.jpg", dpi=300)
