import os
import pandas as pd
import matplotlib.pyplot as plt

data_path = "./data/"


power_color_palete = [
    "#005d5d",
    "#007d79",
    "#08bdba",
    "#3ddbd9",
]

ray_color_palete = ["#491d8b", "#6929c4", "#8a3ffc", "#a56eff"]


def gen_biggest_plot():
    biggest_ray_data = []
    biggest_power_data = []
    for file in os.listdir(data_path):
        if file.startswith("biggest") and file.endswith(".csv"):
            file_path = os.path.join(data_path, file)
            df = pd.read_csv(file_path, delimiter=";")

            if "ray" in file.lower():
                biggest_ray_data.append(df)
            elif "power" in file.lower():
                biggest_power_data.append(df)

    plt.figure(figsize=(10, 6))

    for i, df in enumerate(biggest_ray_data):
        label = "Ray $" + "\lambda_" + str(i + 1) + "$"
        plt.plot(
            df["N"], df["Approx"], color=ray_color_palete[i], alpha=0.6, label=label
        )

    for i, df in enumerate(biggest_power_data):
        label = "Power $" + "\lambda_" + str(i + 1) + "$"
        plt.plot(
            df["N"], df["Approx"], color=power_color_palete[i], alpha=0.6, label=label
        )

    plt.title(
        "Porównanie Metofy Rayleigha i Metody Potęgowej dla największych wartości własnych"
    )
    plt.xlabel("N")
    plt.ylabel("Approx")
    plt.yscale("log")
    # plt.xscale("log")

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlim(0)

    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.legend()

    plt.subplots_adjust(
        left=0.1,  # lewy margines (0–1)
        right=0.9,  # prawy margines
        top=0.95,  # górny margines
        bottom=0.1,  # dolny margines
    )

    plt.savefig("./figures/biggest_ev_cmp_plot.jpg", dpi=300)


def gen_smallest_plot():
    smallest_ray_data = []
    smallest_power_data = []
    for file in os.listdir(data_path):
        if file.startswith("smallest") and file.endswith(".csv"):
            file_path = os.path.join(data_path, file)
            df = pd.read_csv(file_path, delimiter=";")

            if "ray" in file.lower():
                smallest_ray_data.append(df)
            elif "power" in file.lower():
                smallest_power_data.append(df)

    plt.figure(figsize=(10, 6))

    for i, df in enumerate(smallest_ray_data):
        label = "Ray $" + "\lambda_" + str(i + 1) + "$"
        plt.plot(
            df["N"], df["Approx"], color=ray_color_palete[i], alpha=0.6, label=label
        )

    for i, df in enumerate(smallest_power_data):
        label = "Power $" + "\lambda_" + str(i + 1) + "$"
        plt.plot(
            df["N"], df["Approx"], color=power_color_palete[i], alpha=0.6, label=label
        )

    plt.title(
        "Porównanie Metofy Rayleigha i Metody Potęgowej dla najmniejszych wartości własnych"
    )
    plt.xlabel("N")
    plt.ylabel("Approx")
    plt.yscale("log")
    # plt.xscale("log")

    ax = plt.gca()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    plt.xlim(0)

    plt.grid(True, linestyle="--", linewidth=0.5)
    plt.legend()

    plt.subplots_adjust(
        left=0.1,  # lewy margines (0–1)
        right=0.9,  # prawy margines
        top=0.95,  # górny margines
        bottom=0.1,  # dolny margines
    )

    plt.savefig("./figures/smallest_ev_cmp_plot.jpg", dpi=300)


if __name__ == "__main__":
    gen_biggest_plot()
    gen_smallest_plot()
