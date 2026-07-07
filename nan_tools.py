import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

def plot_confusion(metrics, comparison, ax=None):

    green = "#c9e6b0"
    red = "#f2b6b0"
    orange = "#f7d39d"
    cell_font = 18

    # unpack metrics
    tp = metrics["tp"]
    fp = metrics["fp"]
    fn = metrics["fn"]
    precision = metrics["precision"]
    recall = metrics["recall"]
    f1 = metrics["f1"]

    cells = [
    # (row, col, label, sublabel, color)
    (0, 0, "True\nPositive", f"{tp:,}", green),
    (0, 1, "False\nPositive", f"{fp:,}", orange),
    (1, 0, "False\nNegative", f"{fn:,}", red),
    (1, 1, "True\nNegative", "N/A", green)]

    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6, 6))

    ax.set_xlim(0, 3.55)
    ax.set_ylim(0, 3.75)
    ax.set_aspect("equal")
    ax.axis("off")

    # layout
    left, top = 0.5, 3.25 
    cell_size = 1.5

    # make cells
    for r, c, label, sub, color in cells:
        x = left + c * cell_size
        y = top - (r + 1) * cell_size
        ax.text(x + cell_size / 2, y + cell_size / 2 + 0.15, label,
                ha="center", va="center", fontsize=cell_font, fontweight="bold")
        ax.text(x + cell_size / 2, y + cell_size / 2 - 0.2, sub,
                ha="center", va="center", fontsize=cell_font)
        ax.add_patch(Rectangle((x, y), cell_size, cell_size,
                                facecolor=color, edgecolor="none"))

    # labels
    ax.text(left + cell_size, top + 0.55, f"{comparison} Annotation",
            ha="center", va="center", fontsize=cell_font +2, fontweight="bold")


    ax.text(1.75, top - 2*cell_size - 0.2, f"Precision {precision:.4f} | Recall: {recall:.4f} | F1: {f1:.4f}",
            ha="center", va="center", fontsize=12)

    ax.text(left + cell_size / 2, top + 0.2, "Positive", ha="center", va="center", fontsize=cell_font)
    ax.text(left + 1.5 * cell_size, top + 0.2, "Negative", ha="center", va="center", fontsize=cell_font)

    # row labels
    ax.text(left - 0.55, top - cell_size, "AlphaFunctor Prediction",
            ha="center", va="center", fontsize=cell_font +2, fontweight="bold", rotation=90)

    ax.text(left - 0.2, top - cell_size / 2, "Positive", ha="center", va="center", fontsize=cell_font, rotation=90)
    ax.text(left - 0.2, top - 1.5 * cell_size, "Negative", ha="center", va="center", fontsize=cell_font, rotation=90)

    # any and all black lines

    # outline
    lw = 2
    # vertical divider
    ax.plot([left + cell_size, left + cell_size], [top - 2*cell_size, top], color="black", linewidth=lw)
    # horizontal divider
    ax.plot([left, left + 2*cell_size], [top - cell_size, top - cell_size], color="black", linewidth=lw)
    # outer rectangle
    ax.add_patch(Rectangle((left, top - 2 * cell_size), 2 * cell_size, 2 * cell_size,
                            fill=False, edgecolor="black", linewidth=lw))

    # colunm separator
    ax.plot([left + cell_size, left + cell_size], [top, top + 0.4], color="black", linewidth=1)
    ax.plot([left, left], [top, top + 0.4], color="black", linewidth=1)
    ax.plot([left + 2*cell_size, left + 2*cell_size], [top, top + 0.4], color="black", linewidth=1)


    # row separator
    ax.plot([left - 0.4, left], [top - cell_size , top- cell_size], color="black", linewidth=1)
    ax.plot([left - 0.4, left], [top , top], color="black", linewidth=1)
    ax.plot([left - 0.4, left], [top - 2*cell_size, top - 2*cell_size], color="black", linewidth=1)

    # cap row and column
    ax.plot([left, left+2*cell_size], [top + 0.4, top + 0.4], color="black", linewidth=1)
    ax.plot([left - 0.4, left - 0.4], [top, top - 2*cell_size], color="black", linewidth=1)

    if standalone:
        plt.tight_layout()
        plt.savefig(f"{comparison}_confusion.png", bbox_inches="tight")
        plt.savefig(f"{comparison}_confusion.pdf", bbox_inches="tight")

import os
from collections import Counter

def find_missing():
    # this will hopefully find what outputs we're missing
    files = os.listdir('/outputs')
    numbers = sorted(int(re.search(r'_(\d+)\.csv$', f).group(1)) for f in files if re.search(r'_(\d+)\.csv$', f))

    full_range = set(range(numbers[0], numbers[-1] + 1))
    missing = sorted(full_range - set(numbers))
    counts = Counter(numbers)
    duplicates = sorted(n for n, c in counts.items() if c > 1)
    
    print(f"Missing: {missing}")
    print(len(missing))
    print(f"Duplicates: {duplicates}")

import pandas as pd

def quick_vs_uni(quick_data, uni_data):
    q_df = pd.read_csv(quick_data)
    u_df = pd.read_csv(uni_data)

    q_df = q_df.rename(columns={'go_terms': 'quick_terms'})
    u_df = u_df.rename(columns={'go_terms': 'uni_terms'})
    vs_data = pd.merge(u_df, q_df[['uniprot_id', 'quick_terms']], on='uniprot_id', how='left')
    vs_data.to_csv('quick_vs_uni.csv', index=False)

quick_vs_uni("quickgo.csv", "annotated_swissprotkb.csv")

import pickle 
with open("go_cache.pkl", "rb") as f:       # holds information about each go_term
     cache = pickle.load(f)


def go_term_lookup(go_id):
    print(cache[go_id])

# go_term_lookup("GO:0046686")


       
