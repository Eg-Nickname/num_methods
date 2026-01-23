import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

GRAY = "#666666"
pd.options.mode.use_inf_as_na = True


def gen_plot():
    plt.figure(figsize=(12, 6))

    x_smooth = np.linspace(-math.pi / 2, math.pi / 2, 500)
    s = np.sin(x_smooth)
    y_smooth = np.exp(s * s)

    plt.plot(
        x_smooth,
        y_smooth,
        label=r"$f(x) = exp(sin(x)^2)$",
        color="#333333",
        linestyle="-",
        alpha=0.4,
        zorder=1,
    )

    df = pd.read_csv("./data/trapezoid_points.csv", sep=";")

    plt.plot(
        df["x"],
        df["y"],
        label="Interpolacja",
        color="#6929c4",
        alpha=0.7,
        zorder=2,
        linewidth=1,
    )

    plt.scatter(
        df["x"],
        df["y"],
        color="#6929c4",
        s=8,
        label="Obliczone wartości funkcji",
        zorder=3,
    )

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("x")
    plt.ylabel("f(x)")

    plt.xlim(-math.pi / 2 - 0.1, math.pi / 2 + 0.1)
    # plt.ylim(-0.2, 1.2)

    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.legend()

    plt.subplots_adjust(
        left=0.10,  # lewy margines (0–1)
        right=0.90,  # prawy margines
        top=0.9,  # górny margines
        bottom=0.12,  # dolny margines
    )

    plt.savefig("./figures/trapezoid_plot.jpg", dpi=300)


if __name__ == "__main__":
    gen_plot()
