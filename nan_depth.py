import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import ast
import sklearn as sk
from collections import Counter
import seaborn as sns

with open("go_cache.pkl", "rb") as f:
    cache = pickle.load(f)

plt.rcParams.update({
    'font.size': 24,
    'axes.titlesize': 24,
    'axes.labelsize': 24,
    'xtick.labelsize': 20,
    'ytick.labelsize': 20,
})

# would it make more sense to make just a thing of subplots with :
# depth/level/Ufreq/Afreq vs P/R/f1 

# maybe here we can have a subplot with depth and level and f1? rn

def get_vecs():
    all_terms = set(cache.keys())

    f1_vec = []
    depth_vec = []
    level_vec = []
    
    for term in all_terms:
        if "stats" not in cache[term]:
            continue
        f1 = cache[term]["stats"]["f1"]
        depth = cache[term]["depth"]
        level = cache[term]["level"]

        f1_vec.append(f1)
        depth_vec.append(depth)
        level_vec.append(level)

    return f1_vec, depth_vec, level_vec

def DL_percentages():

    all_terms = set(cache.keys())

    depths = [cache[term]["depth"] for term in all_terms]
    levels = [cache[term]["level"] for term in all_terms]

    depth_pct = pd.Series(Counter(depths)).sort_index() / len(depths) * 100
    level_pct = pd.Series(Counter(levels)).sort_index() / len(levels) * 100

    return depth_pct, level_pct

def main():
    f1_vec, depth_vec, level_vec = get_vecs()

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    df_plot = pd.DataFrame({"depth": depth_vec, "level": level_vec, "f1": f1_vec})

    sns.boxplot(data=df_plot, x="depth", y="f1", ax=axes[0])
    axes[0].set_xlabel('Depth of GoTerm')
    axes[0].set_ylabel('F1 Value')
    axes[0].set_title('Depth of GoTerm vs F1')

    sns.boxplot(data=df_plot, x="level", y="f1", ax=axes[1])
    axes[1].set_xlabel('Level of GoTerm')
    axes[1].set_ylabel('F1 Value')
    axes[1].set_title('Level of GoTerm vs F1')

    plt.tight_layout()
    plt.savefig("DepthDensity_big.png")
    plt.show()

# ok just the level boxplots

def just_level():
    f1_vec, depth_vec, level_vec = get_vecs()
    fig, ax = plt.subplots(figsize=(10, 8))

    df_plot = pd.DataFrame({"depth": depth_vec, "level": level_vec, "f1": f1_vec})
    depth_pct, level_pct = DL_percentages()

    sns.boxplot(data=df_plot, x="level", y="f1", ax=ax)
    ax.set_xlabel('Level of GoTerm')
    ax.set_ylabel('F1 Value')
    ax.set_title('Level of GoTerm vs F1')
    y_max = df_plot["f1"].max()
    y_min = df_plot["f1"].min()
    offset = (y_max - y_min) * 0.03

    # categories = sorted(df_plot["level"].unique())
    # for i, cat in enumerate(categories):
    #     box_top = df_plot.loc[df_plot["level"] == cat, "f1"].max()
    #     pct = level_pct.get(cat, 0)
    #     ax.text(i, box_top + offset, f"{pct:.2f}%",
    #             ha="center", va="bottom", fontsize=12, color="dimgray")

    # ax.set_ylim(top=y_max + offset * 6)

    categories = sorted(df_plot["level"].unique())
    for i, cat in enumerate(categories):
        group = df_plot.loc[df_plot["level"] == cat, "f1"]
        q1 = group.quantile(0.25)
        q2 = group.quantile(0.5)
        box_center = (q1 + q2) / 2
        y_pos = q1 + 0.09
        pct = level_pct.get(cat, 0)
        ax.text(i, y_pos, f"{pct:.3f}%",
                ha="center", va="center", fontsize=14, color="white", fontweight="bold",
                rotation=90)
                # bbox=dict(facecolor="white", alpha=0.6, edgecolor="none", pad=1))


    plt.savefig("LevelDensity_brightbig.png")

just_level()