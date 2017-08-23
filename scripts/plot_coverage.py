import matplotlib
matplotlib.use("agg")

import pandas
from itertools import islice
from collections import deque

import matplotlib.pyplot as plt
plt.style.use("ggplot")


def _sliding_window(iterable, size=3, step=1, fillvalue=None):

    if size < 0 or step < 1:
        raise ValueError
    it = iter(iterable)
    q = deque(islice(it, size), maxlen=size)
    if not q:
        return
    q.extend(fillvalue for _ in range(size - len(q)))
    while True:
        yield iter(q)
        try:
            q.append(next(it))
        except StopIteration:
            return
        q.extend(next(it, fillvalue) for _ in range(step - 1))


def _smooth_coverage(df, size=10000, mode="median"):

    coverage_summary = []

    for win in _sliding_window(df["Depth"], size=size, step=size):
        win = pandas.Series(list(win))

        if mode == "median":
            value = win.median()
        elif mode == "mean":
            value = win.mean()
        else:
            raise ValueError("Mode not supported.")

        coverage_summary.append(value)

    return pandas.DataFrame({"Coverage": coverage_summary,
                             "Reference": [i * size for i in range(len(coverage_summary))]})


def plot_coverage(file, output, window_size=20000, mode="median"):

    df = pandas.read_csv(file, sep="\t", header=None, names=["Chromosome", "Position", "Depth"])

    if mode == "mean":
        coverage = df.groupby("Chromosome").mean()
    else:
        coverage = df.groupby("Chromosome").median()

    print("{mode} depth of coverage per chromosome is:\n\n".format(mode=mode.capitalize()), coverage["Depth"], "\n")

    chromosomes = []
    for chromosome in df["Chromosome"].unique():
        chr_df = df[df["Chromosome"] == chromosome]
        chromosomes.append(_smooth_coverage(chr_df, size=window_size, mode=mode))

    n_chr = len(chromosomes)

    fig, axs = plt.subplots(ncols=n_chr)

    colors = ["#762a83", "#1b7837"]

    for i in range(n_chr):
        chrom = "Chromosome {chrom}".format(chrom=i+1)

        ax = chromosomes[i].plot(x="Reference", y="Coverage", ax=axs[i], legend=False,
                                 title=chrom, color=colors[i], rot=30)

        ax.set(xlabel="Reference position (bp)", ylabel="Depth of coverage (x)")
        ax.set_ylim(0)

    fig.set_size_inches(12.7, 8.7)
    fig.savefig(output)

plot_coverage(file=snakemake.input[0], output=snakemake.output[0], window_size=snakemake.params[0], mode="mean")
plot_coverage(file=snakemake.input[0], output=snakemake.output[1], window_size=snakemake.params[0], mode="median")