import numpy as np
import matplotlib.pyplot as plt

def read_data(path):
    res = []
    with open(path, "r") as file:
        for line in file.readlines()[1:]:
            res.append(tuple(map(int, line.strip().split(","))))
    return res

def create_plot(output_path, plot_type="semi-log", input_path="solution_counting-full.csv"):
    # Read data from .csv provided
    raw_data = read_data(input_path)

    # Get type-1 and type-2 values separately
    xs, ys_t1 = zip(*sorted([(n, t1) for n, _, t1, _, _ in raw_data]))
    xs, ys_t2 = zip(*sorted([(n, t2) for n, _, _, t2, _ in raw_data]))

    # Make the argument None (or something else) to get a normal plot
    if plot_type == "semi-log":
        xs = np.log10(xs)

    fig, ax = plt.subplots(figsize=(12, 12), layout="constrained", dpi=300)

    # Title and axis labels
    ax.set_title("Solution-count by value and type", fontsize=16)
    ax.set_xlabel("x", fontsize=12)
    ax.set_ylabel("f(x)", fontsize=12)

    # Scatter type-1 and type-2 solution data
    ax.scatter(xs, ys_t1, color="blue", alpha=0.06, s=1)
    ax.scatter(xs, ys_t2, color="orange", alpha=0.06, s=1)

    # Dummies with high alpha for the legend
    t1_dummy = ax.scatter([], [], color="blue", alpha=0.7, s=14, label="Type-1")
    t2_dummy = ax.scatter([], [], color="orange", alpha=0.7, s=14, label="Type-2")
    ax.legend()

    # Grid setup
    ax.minorticks_on()
    plt.grid(which='major', alpha=0.3)
    plt.grid(which='minor', alpha=0.15)

    # Save the file
    plt.savefig(output_path)


create_plot("Figure1.png", "semi-log")
create_plot("Figure2.png", "linear")
