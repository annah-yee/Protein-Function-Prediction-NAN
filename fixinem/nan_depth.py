import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import ast
import sklearn as sk
from collections import Counter
from goatools.obo_parser import GODag
import seaborn as sns

with open("go_cache.pkl", "rb") as f:
    cache = pickle.load(f)

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

def main():
    f1_vec, depth_vec, level_vec = get_vecs()

    fig, axes = plt.subplots(1, 2, figsize=(20, 8))
    df_plot = pd.DataFrame({"depth": depth_vec, "level": level_vec, "f1": f1_vec})

    # print(df_plot.describe())

    sns.violinplot(data=df_plot, x="depth", y="f1", ax=axes[0])
    axes[0].set_xlabel('Depth of GoTerm')
    axes[0].set_ylabel('F1 Value')
    axes[0].set_title('Depth of GoTerm vs F1')

    sns.violinplot(data=df_plot, x="level", y="f1", ax=axes[1], cut=0)
    axes[1].set_xlabel('Level of GoTerm')
    axes[1].set_ylabel('F1 Value')
    axes[1].set_title('Level of GoTerm vs F1')

    plt.tight_layout()
    plt.savefig("DepthDensity.png")
    plt.show()

main()