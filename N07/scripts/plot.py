import os
import pandas as pd
import matplotlib.pyplot as plt


df = pd.read_csv("./data/u_n.csv", delimiter=";")

plt.figure(figsize=(10, 6))
plt.plot(df["u_n"], label="$u_n$", color="#6929c4")
plt.xlabel("N")
plt.ylabel("$u_n$")


ax = plt.gca()
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)


plt.ylim(0.0, 1.05)
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
