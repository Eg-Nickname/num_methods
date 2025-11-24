import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

plt.figure(figsize=(12, 6))
# Pallete https://carbondesignsystem.com/data-visualization/color-palettes/#categorical-palettes
color_palete = [
    # TEAL
    "#005d5d",  # TBD
    "#007d79",
    "#009d9a",
    "#08bdba",
    "#3ddbd9",
    # PURPLE
    "#8a3ffc",
    "#a56eff",
]

# GRAY
GRAY = "#666666"

# Execute for all datapoints in benchdata
ax = plt.gca()
files = [
    "./bench_data/methods_10k/optimised_fullpiv_dense_qr_10k.csv",
    "./bench_data/methods_10k/optimised_fullpiv_dense_lu_10k.csv",
    "./bench_data/methods_10k/optimised_sparse_qr_10k.csv",
    "./bench_data/methods_10k/optimised_par_dense_qr_10k.csv",
    "./bench_data/methods_10k/optimised_par_dense_lu_10k.csv",
    "./bench_data/methods_10k/optimised_sparse_lu_10k.csv",
    "./bench_data/methods_10k/optimised_sherman_10k.csv",
]

i = 0
y_offsets = [70000000, -70000000, 7000000, -8000000, -4000000, 0, 0]
for csv_file in files:
    df = pd.read_csv(csv_file, sep=";")

    x = df["n"].values
    y_columns = [col for col in df.columns if col != "n"]
    y = df[y_columns[0]].values
    y_pos = df[y_columns[0]].iloc[-1] + y_offsets[i]
    ax.text(
        df["n"].iloc[-1] + 100,
        y_pos,
        y_columns[0],
        color=color_palete[i],
        fontsize=10,
        va="center",
    )
    plt.scatter(x, y, s=15, color=color_palete[i], label=f"{y_columns[0]}")
    i += 1


plt.xlim(0, 9700)
plt.xticks(np.arange(100, 9341, 660))

ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)


plt.xlabel("N")
plt.ylabel("Time μs")
plt.yscale("log")

plt.grid(True, linestyle="--", linewidth=0.5)
ax.axhline(10**3, linewidth=0.7, linestyle="--", color="#222222")
ax.axhline(10**6, linewidth=0.7, linestyle="--", color="#222222")
ax.axhline(60 * 10**6, linewidth=0.7, linestyle="--", color="#222222")

ax.text(10, 10**3 + 3 * 10**2, "1ms", va="center", ha="left")
ax.text(10, 10**6 + 3 * 10**5, "1s", va="center", ha="left")
ax.text(10, 60 * 10**6 + 20 * 10**6, "1m", va="center", ha="left")


# plt.legend()
plt.subplots_adjust(
    left=0.05,  # lewy margines (0–1)
    right=0.85,  # prawy margines
    top=0.95,  # górny margines
    bottom=0.1,  # dolny margines
)

plt.savefig("figures/combined_10k.jpg", dpi=300)
