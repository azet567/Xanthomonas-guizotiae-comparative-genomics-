import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


# ============================
# Read FastANI pairwise output
# ============================

df = pd.read_csv(
    "ANI_matrix.tsv",
    sep="\t",
    header=None,
    names=["query", "reference", "ANI", "fragments", "total"]
)

print("ANI comparisons:", df.shape)


# ============================
# Build ANI matrix
# ============================

genomes = sorted(
    set(df["query"]).union(set(df["reference"]))
)

ani = pd.DataFrame(
    100.0,
    index=genomes,
    columns=genomes
)


for _, row in df.iterrows():
    ani.loc[row["query"], row["reference"]] = row["ANI"]
    ani.loc[row["reference"], row["query"]] = row["ANI"]


# Remove .fna
ani.index = ani.index.str.replace(".fna", "", regex=False)
ani.columns = ani.columns.str.replace(".fna", "", regex=False)


print("Final ANI matrix:", ani.shape)


# ============================
# Clustered heatmap
# ============================

sns.set(font_scale=0.7)

g = sns.clustermap(
    ani,
    method="average",
    metric="euclidean",
    cmap="RdYlBu_r",
    vmin=85,
    vmax=100,
    figsize=(16,14),
    linewidths=0.2,
    xticklabels=True,
    yticklabels=True,
    dendrogram_ratio=(0.15,0.15),
    cbar_kws={
        "label":"Average Nucleotide Identity (%)"
    }
)


g.ax_heatmap.set_xticklabels(
    g.ax_heatmap.get_xticklabels(),
    rotation=90,
    fontsize=7
)

g.ax_heatmap.set_yticklabels(
    g.ax_heatmap.get_yticklabels(),
    fontsize=7
)


# ============================
# Save figures
# ============================

g.savefig(
    "ANI_heatmap_Xanthomonas.png",
    dpi=600,
    bbox_inches="tight"
)

g.savefig(
    "ANI_heatmap_Xanthomonas.pdf",
    bbox_inches="tight"
)

g.savefig(
    "ANI_heatmap_Xanthomonas.svg",
    bbox_inches="tight"
)


print("✅ ANI heatmap generated")
