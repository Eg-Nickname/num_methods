import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

color_palete = {
    "Dense Full Pivot QR": "#005d5d",
    "Dense Full Pivot LU": "#007d79",
    "Sparse QR": "#6929c4",
    "Dense Partial QR": "#08bdba",
    "Dense Partial Pivot LU": "#3ddbd9",
    "Sparse LU": "#8a3ffc",
    "Sherman Morrison": "#a56eff",
}


def plot_csv(csv_path):
    df = pd.read_csv(csv_path, sep=";")

    x = df["n"].values
    y_columns = [col for col in df.columns if col != "n" and col != "Unnamed"]
    y = df[y_columns[0]].values

    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, s=20, color=color_palete[y_columns[0]])

    plt.xlim(0, 9700)
    plt.xticks(np.arange(100, 9341, 660))

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.title(f"{y_columns[0]}")

    plt.xlabel("N")
    plt.ylabel("Time μs")
    plt.yscale("log")

    plt.grid(True, linestyle="--", linewidth=0.5)
    # Horizontal lines
    ax.axhline(10**3, linewidth=0.7, linestyle="--", color="#222222")
    if y[-1] < 10**6:
        ax.text(10, 10**3 + 1 * 10**2, "1ms", va="center", ha="left")
    else:
        ax.text(10, 10**3 + 3 * 10**2, "1ms", va="center", ha="left")
    if y[-1] > 10**6:
        ax.text(10, 10**6 + 3 * 10**5, "1s", va="center", ha="left")
        ax.axhline(10**6, linewidth=0.7, linestyle="--", color="#222222")
    if y[-1] > 60 * 10**6:
        ax.axhline(60 * 10**6, linewidth=0.7, linestyle="--", color="#222222")
        ax.text(10, 60 * 10**6 + 20 * 10**6, "1m", va="center", ha="left")

    filename = y_columns[0].replace(" ", "_").lower()
    plt.savefig(f"figures/10k_sep/{filename}.jpg", dpi=300)


# Execute for all datapoints in benchdata
for csv_file in glob.glob("./bench_data/methods_10k/*.csv"):
    print(csv_file)
    plot_csv(csv_file)
