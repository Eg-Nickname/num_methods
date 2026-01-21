import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

GRAY = "#666666"
pd.options.mode.use_inf_as_na = True


def gen_inter_plot(nodes_file, polynomial_file, node_color, polynomial_color, filename):
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
    df_nodes = pd.read_csv(nodes_file, sep=";")

    # Plot dense polynomial points
    df_polynomial = pd.read_csv(polynomial_file, sep=";")

    plt.scatter(
        df_nodes["x"],
        df_nodes["y"],
        color=node_color,
        label="Węzły Interpolacji",
        zorder=3,
    )

    plt.plot(
        df_polynomial["x"],
        df_polynomial["y"],
        label="Wielomian Interpolacji",
        color=polynomial_color,
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
        if abs(prev_tick - x) > 0.30 and abs(x) > 0.1:
            ticks.append(x)
            ticks.append(-x)
            prev_tick = x
    plt.xticks(ticks, rotation=45)

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    plt.xlabel("x")
    plt.ylabel("y(x)")

    plt.xlim(-5.5, 5.5)
    plt.ylim(-0.2, 1.2)

    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.legend()

    plt.subplots_adjust(
        left=0.10,  # lewy margines (0–1)
        right=0.90,  # prawy margines
        top=0.9,  # górny margines
        bottom=0.12,  # dolny margines
    )

    plt.savefig(filename, dpi=300)


def gen_error_plot():
    plt.figure(figsize=(12, 6))

    # Plot interpolation nodes
    df_czebyszew_polynomial_err = pd.read_csv(
        "./data/errors_czybyszew_polynomial.csv", sep=";"
    )
    df_eq_czebyszew_polynomial_err = df_czebyszew_polynomial_err[
        df_czebyszew_polynomial_err["N"] % 2 == 0
    ]
    df_equal_polynomial_err = pd.read_csv("./data/errors_equal_polynomial.csv", sep=";")

    df_czebyszew_direct_err = pd.read_csv("./data/errors_czebyszew_direct.csv", sep=";")

    plt.plot(
        df_eq_czebyszew_polynomial_err["N"],
        df_eq_czebyszew_polynomial_err["Error"],
        label="Błędy dla wielomianu z węzłów Czebyszewa",
        color="#6929c4",
        alpha=1,
        zorder=2,
        linewidth=2,
    )
    plt.plot(
        df_equal_polynomial_err["N"],
        df_equal_polynomial_err["Error"],
        label="Błędy dla wielominau z równoodległych wezłów",
        color="#3ddbd9",
        alpha=1,
        zorder=2,
        linewidth=2,
    )
    plt.plot(
        df_czebyszew_direct_err["N"],
        df_czebyszew_direct_err["Error"],
        label="Błędy dla bezpośrednioliczonych węzłów czebyszewa",
        color="#0072c3",
        alpha=1,
        zorder=2,
        linewidth=2,
    )

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlim(1)

    plt.xlabel("N")
    plt.ylabel("$\Delta f$")

    plt.yscale("log")
    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.legend()

    plt.subplots_adjust(
        left=0.10,  # lewy margines (0–1)
        right=0.90,  # prawy margines
        top=0.9,  # górny margines
        bottom=0.1,  # dolny margines
    )

    plt.savefig("./figures/errors.jpg", dpi=300)


if __name__ == "__main__":
    gen_inter_plot(
        "./data/nodes_equal_polynomial.csv",
        "./data/equal_polynomial.csv",
        "#08bdba",
        "#3ddbd9",
        "./figures/equal_plot.jpg",
    )

    gen_inter_plot(
        "./data/nodes_czybyszew_polynomial.csv",
        "./data/czybyszew_polynomial.csv",
        "#6929c4",
        "#8a3ffc",
        "./figures/czebyszew_plot.jpg",
    )
    gen_inter_plot(
        "./data/nodes_runge_polynomial.csv",
        "./data/runge_polynomial.csv",
        "#005d5d",
        "#007d79",
        "./figures/runge_plot.jpg",
    )
    gen_error_plot()

    gen_inter_plot(
        "./data/nodes_czebyszew_direct.csv",
        "./data/czebyszew_direct.csv",
        "#0072c3",
        "#0072c3",
        "./figures/direct_czebyszew_plot.jpg",
    )
