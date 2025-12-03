import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# GRAY
GRAY = "#666666"

# Execute for all datapoints in benchdata


color_palete = {
    "Dense Full Pivot LU No AVX": "#007d79",
    "Dense Partial QR No AVX": "#08bdba",
    "Dense Full Pivot LU AVX": "#6929c4",
    "Dense Partial QR AVX": "#8a3ffc",
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
            df["n"].iloc[-1] + 100,
            y_pos,
            y_columns[0],
            color=color_palete[y_columns[0]],
            fontsize=14,
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

    plt.xlim(x_start, 4000)
    plt.ylim(y_min * 0.8)
    plt.xticks(np.arange(x_start, 4000, 400) + 100)

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("N")
    plt.ylabel("Time μs")
    plt.yscale("log")

    plt.grid(True, linestyle="--", linewidth=0.5)
    ax.axhline(10**3, linewidth=0.7, linestyle="--", color="#222222")
    ax.axhline(10**6, linewidth=0.7, linestyle="--", color="#222222")
    ax.axhline(60 * 10**6, linewidth=0.7, linestyle="--", color="#222222")
    ax.axhline(600 * 10**6, linewidth=0.7, linestyle="--", color="#222222")

    ax.text(x_start + 10, 10**3 + 3 * 10**2, "1ms", va="center", ha="left", fontsize=14)
    ax.text(x_start + 10, 10**6 + 3 * 10**5, "1s", va="center", ha="left", fontsize=14)
    ax.text(
        x_start + 10,
        60 * 10**6 + 20 * 10**6,
        "1m",
        va="center",
        ha="left",
        fontsize=14,
    )
    ax.text(
        x_start + 10,
        600 * 10**6 + 20 * 10**7,
        "10m",
        va="center",
        ha="left",
        fontsize=14,
    )

    # plt.legend()
    plt.subplots_adjust(
        left=0.05,  # lewy margines (0–1)
        right=0.77,  # prawy margines
        top=0.95,  # górny margines
        bottom=0.1,  # dolny margines
    )
    plt.savefig(output_file, dpi=300)


# Generate chart for all data
all_y_offsets = [0, 0, 0, -50000000]
all_files = [
    "./bench_data/optimisation_cmp/march_fullpiv_dense_fullpiv_lu_4k.csv",
    "./bench_data/optimisation_cmp/march_fullpiv_dense_partialpiv_qr_4k.csv",
    "./bench_data/optimisation_cmp/no_march_dense_fullpiv_lu_4000.csv",
    "./bench_data/optimisation_cmp/no_march_dense_partialpiv_qr_4000.csv",
]
gen_combined(all_files, "figures/march_comparison_4k.jpg", all_y_offsets)
