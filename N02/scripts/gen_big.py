import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# GRAY
GRAY = "#666666"

# Execute for all datapoints in benchdata


color_palete = {
    "Dense Full Pivot QR": "#005d5d",
    "Dense Full Pivot LU": "#007d79",
    "Sparse QR": "#6929c4",
    "Dense Partial QR": "#08bdba",
    "Dense Partial Pivot LU": "#3ddbd9",
    "Sparse LU": "#8a3ffc",
    "Sherman Morrison": "#a56eff",
}


def gen_combined(input_files, output_file, offsets, x_start=0):
    plt.figure(figsize=(12, 6))
    ax = plt.gca()
    i = 0
    y_min = 10**10
    for csv_file in input_files:
        df = pd.read_csv(csv_file, sep=";")

        x = df["n"].values
        y_columns = [col for col in df.columns if col != "n"]
        y = df[y_columns[0]].values
        y_pos = df[y_columns[0]].iloc[-1] + offsets[i]
        ax.text(
            df["n"].iloc[-1] + 1000,
            y_pos,
            y_columns[0],
            color=color_palete[y_columns[0]],
            fontsize=10,
            va="center",
        )
        mask = x > x_start
        x_cut = x[mask]
        y_cut = y[mask]
        if y_cut[0] < y_min:
            y_min = y_cut[0]
        plt.scatter(
            x, y, s=15, color=color_palete[y_columns[0]], label=f"{y_columns[0]}"
        )
        i += 1

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("N")
    plt.ylabel("Time μs")
    plt.yscale("log")
    plt.xlim(
        x_start,
    )
    plt.xlim(0, 98500)
    plt.xticks(np.arange(100, 100000, 999 * 6))

    plt.grid(True, linestyle="--", linewidth=0.5)
    ax.axhline(10**3, linewidth=0.7, linestyle="--", color="#222222")
    ax.axhline(10**6, linewidth=0.7, linestyle="--", color="#222222")

    ax.text(x_start + 10, 10**3 + 3 * 10**2, "1ms", va="center", ha="left")
    ax.text(x_start + 10, 10**6 + 3 * 10**5, "1s", va="center", ha="left")

    # plt.legend()
    plt.subplots_adjust(
        left=0.05,  # lewy margines (0–1)
        right=0.85,  # prawy margines
        top=0.95,  # górny margines
        bottom=0.1,  # dolny margines
    )
    plt.savefig(output_file, dpi=300)


big_files = [
    "./bench_data/optimised_sparse_lu_100k.csv",
    "./bench_data/optimised_sherman_100k.csv",
]
gen_combined(big_files, "figures/big_mat_100k.jpg", [0, 0])
