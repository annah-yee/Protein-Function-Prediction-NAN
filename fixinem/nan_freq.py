import pickle
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
import ast
import sklearn as sk
from collections import Counter

tqdm.pandas()

def retrieve(row):
    '''Return basic stats per row'''
    try:
        quick_set = set(ast.literal_eval(row["go_terms"]))
        pred_set = set(ast.literal_eval(row["predictions"]))

        return pd.Series({
            "tp": quick_set & pred_set,
            "fn": quick_set - pred_set,
            "fp": pred_set - quick_set
        })
    except (ValueError, SyntaxError):
        return pd.Series({"tp": 0, "fn": 0, "fp": 0})

def seek(df):

    with open("go_cache.pkl", "rb") as f:
        cache = pickle.load(f)
    
    all_terms = set(cache.keys())
    
    precision_vec = []
    recall_vec = []
    quick_count_vec = []
    pred_count_vec = []
    f1_vec = []

    tp_counts = Counter(term for s in df['tp'] if isinstance(s, set) for term in s)
    fp_counts = Counter(term for s in df['fp'] if isinstance(s, set) for term in s)
    fn_counts = Counter(term for s in df['fn'] if isinstance(s, set) for term in s)

    for term in all_terms:
        tp = tp_counts[term]
        fp = fp_counts[term]
        fn = fn_counts[term]

        if tp + fp + fn == 0:
            continue

        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = (2 * precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

        quick_count = tp + fn
        pred_count = tp +fp

        cache[term]["stats"] = {"tp": tp,
                                "fn": fn,
                                "fp": fp,
                                "precision": precision,
                                "recall": recall,
                                "f1": f1,
                                "quick_count": quick_count,
                                "pred_count": pred_count}
        
        precision_vec.append(precision)
        recall_vec.append(recall)
        quick_count_vec.append(quick_count)
        pred_count_vec.append(pred_count)
        f1_vec.append(f1)
        
    with open("go_cache.pkl", "wb") as f:
        pickle.dump(cache, f)
    
    # return 
    return precision_vec, recall_vec, quick_count_vec, pred_count_vec, f1_vec

def pr_main():

    df = pd.read_csv("cleaned_anew.csv")
    stats = df.progress_apply(retrieve, axis=1)
    together = pd.concat([df, stats], axis=1)
    # together.to_csv("together.csv")

    precision_vec, recall_vec, quick_count_vec, pred_count_vec, f1_vec = seek(together)

    fig, ax = plt.subplots(figsize=(10,9))
    hb = ax.hexbin(recall_vec, precision_vec, gridsize=20, mincnt=1, bins="log")
    plt.colorbar(hb, ax=ax, label="count")
    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title('Precision-Recall Density')
    plt.tight_layout()
    plt.savefig("PRDensity.png")
    plt.show()

def freq_main():
    # so here's the thing,, what do we plot frequency against?
    # gotta be something like accuracy (is that what f1 would do?)

    df = pd.read_csv("cleaned_anew.csv")
    stats = df.progress_apply(retrieve, axis=1)
    together = pd.concat([df, stats], axis=1)
    # together.to_csv("together.csv")

    precision_vec, recall_vec, quick_vec, pred_vec, f1_vec = seek(together)

    fig, axes = plt.subplots(1, 2, figsize=(20, 9))

    # ground truth frequency vs f1
    qhb = axes[0].hexbin(quick_vec, f1_vec, gridsize=50, mincnt=1)
    plt.colorbar(qhb, ax=axes[0], label="count")
    axes[0].set_xlabel('UniProt Annotation Frequency')
    axes[0].set_ylabel('F1 Value')
    axes[0].set_title('UniProt Annotation Frequency vs F1')

    # prediction frequency vs f1
    phb = axes[1].hexbin(pred_vec, f1_vec, gridsize=50, mincnt=1)
    plt.colorbar(phb, ax=axes[1], label="count")
    axes[1].set_xlabel('AlphaFunctor Prediction Frequency')
    axes[1].set_ylabel('F1 Value')
    axes[1].set_title('AlphaFunctor Prediction Frequency vs F1')

    plt.tight_layout()
    plt.savefig("FreqDensity.png")
    plt.show()

freq_main()
