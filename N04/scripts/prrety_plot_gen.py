import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# GRAY
GRAY = "#666666"

# Execute for all datapoints in benchdata


color_palete = {
    "Dark Teal": "#08bdba",
    "Light Teal": "#3ddbd9",
    "Dark Purple": "#6929c4",
    "Light Purple": "#8a3ffc",
}


def gen_equal_inter_plot():
    plt.figure(figsize=(12, 6))

    x_smooth = np.linspace(-5, 5, 500)
    y_smooth = 1 / (1 + x_smooth**2)

    plt.plot(
        x_smooth,
        y_smooth,
        label=r"$f(x) = \frac{1}{1+x^2}$",
        color=GRAY,
        linestyle="-",
        alpha=0.4,
        zorder=1,
    )

    # Plot interpolation nodes
    df_nodes = pd.read_csv("./data/nodes_equal_polynomial.csv", sep=";")

    # Plot dense polynomial points
    df_polynomial = pd.read_csv("./data/equal_polynomial.csv", sep=";")

    plt.scatter(
        df_nodes["x"],
        df_nodes["y"],
        color=color_palete["Dark Teal"],
        label="Węzły Interpolacji",
        zorder=3,
    )

    plt.plot(
        df_polynomial["x"],
        df_polynomial["y"],
        label="Wielomian Interpolacji",
        color=color_palete["Light Teal"],
        alpha=0.5,
        zorder=2,
        linewidth=2,
    )

    plt.xticks(df_nodes["x"])
    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("x")
    plt.ylabel("f(x)")

    plt.xlim(-5.5, 5.5)
    plt.ylim(-0.2, 1.2)

    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.title("Interpolacja równoodległymi węzłami")
    plt.legend()

    plt.subplots_adjust(
        left=0.05,  # lewy margines (0–1)
        right=0.85,  # prawy margines
        top=0.95,  # górny margines
        bottom=0.1,  # dolny margines
    )
    plt.savefig("./figures/equal_plot.jpg", dpi=300)


def gen_czybyszew_inter_plot():
    plt.figure(figsize=(12, 6))

    x_smooth = np.linspace(-5, 5, 500)
    y_smooth = 1 / (1 + x_smooth**2)

    plt.plot(
        x_smooth,
        y_smooth,
        label=r"$f(x) = \frac{1}{1+x^2}$",
        color=GRAY,
        linestyle="-",
        alpha=0.4,
        zorder=1,
    )

    # Plot interpolation nodes
    df_nodes = pd.read_csv("./data/nodes_czybyszew_polynomial.csv", sep=";")

    # Plot dense polynomial points
    df_polynomial = pd.read_csv("./data/czybyszew_polynomial.csv", sep=";")

    plt.scatter(
        df_nodes["x"],
        df_nodes["y"],
        color=color_palete["Dark Purple"],
        label="Węzły Interpolacji",
        zorder=3,
    )

    plt.plot(
        df_polynomial["x"],
        df_polynomial["y"],
        label="Wielomian Interpolacji",
        color=color_palete["Light Purple"],
        alpha=0.5,
        zorder=2,
        linewidth=2,
    )

    # Ugly bot working ticks selection
    ticks = []
    prev_tick = -10.0
    for x in df_nodes["x"]:
        if x > 0:
            break
        if abs(prev_tick - x) > 0.3:
            ticks.append(x)
            ticks.append(-x)
            prev_tick = x
    plt.xticks(ticks, rotation=45)

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("x")
    plt.ylabel("f(x)")

    plt.xlim(-5.5, 5.5)
    plt.ylim(-0.2, 1.2)

    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.title("Interpolacja węzłami Czybyszewa")
    plt.legend()

    plt.subplots_adjust(
        left=0.10,  # lewy margines (0–1)
        right=0.90,  # prawy margines
        top=0.9,  # górny margines
        bottom=0.1,  # dolny margines
    )
    plt.savefig("./figures/czybyszew_plot.jpg", dpi=300)


if __name__ == "__main__":
    gen_equal_inter_plot()
    gen_czybyszew_inter_plot()
