import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob

color_palete = {
    "y_n": "#6929c4",
}


def plot_csv(csv_path):
    df = pd.read_csv(csv_path, sep=";")

    x = df["x_n"].values
    y_columns = [col for col in df.columns if col != "n" and col != "Unnamed"]
    y = df[y_columns[1]].values

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, color="#6929c4")

    # x = np.linspace(0, 2, 500)
    # y = -np.cos(6.27 * x) / 39.5 - 0.0188

    # plt.plot(x, y, label="-cos(6.27x)/39.5 - 0.0188", color="#08bdba")
    plt.xlim(0, 2)
    # plt.xticks(np.arange(100, 9341, 660))

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # plt.title(f"{y_columns[1]}")

    plt.xlabel("x_n")
    plt.ylabel("y_n")
    # plt.yscale("log")

    plt.grid(True, linestyle="--", linewidth=0.5)

    filename = y_columns[0].replace(" ", "_").lower()
    plt.savefig(f"figures/sol_chart.jpg", dpi=300)

    x = np.linspace(0, 2, 1000)
    y = -np.cos(6.27 * x) / 39.5 + 0.0248

    plt.plot(x, y, label="-cos(6.27x)/39.5 - 0.0188", color="#08bdba")
    # plt.legend()
    # plt.show()
    plt.savefig(f"figures/sol_chart_with_estimation.jpg", dpi=300)


# Execute for all datapoints in benchdata
plot_csv("./bench_data/solution_data.csv")
