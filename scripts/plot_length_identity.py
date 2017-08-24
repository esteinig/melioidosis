import matplotlib

matplotlib.use('agg')

import pandas
import seaborn as sns

sns.set_style("white")

for i in range(len(snakemake.input)):
    file = snakemake.input[i]
    df = pandas.read_csv(file, sep="\t")

    df = df[df["Length"] > 0]

    p = sns.jointplot(x="Length", y="Identity", kind="scatter", data=df, stat_func=None, space=0, s=0.1)

    p.fig.set_size_inches(12.7, 8.7)
    p.fig.savefig(snakemake.output[i])
