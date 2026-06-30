# Basic Analysis

To run basic analysis, must have the `final_output.csv` (aka all the predictions collected from Alphafunctor) and the annotations collected from UniProt `annotated_swissprotkb.csv`.

somewhere here,, explain the kinds of things that this thing can do


## Import requirements

```python
import pickle
from goatools.obo_parser import GODag
import pandas as pd
import ast
import sklearn as sk
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, recall_score,precision_score, f1_score
import matplotlib.pyplot as plt
from tqdm import tqdm
import matplotlib.gridspec as gridspec
```

## Merge and clean the data

```python
def merge():
    '''Merges data from AlphaFunctor outputs and UniProt annotations'''

    main_data = pd.read_csv("annotated_swissprotkb.csv", sep=',')
    predictions = pd.read_csv("final_output.csv", sep=',')
    merged_data = pd.merge(main_data, predictions[['uniprot_id', 'predictions']], on='uniprot_id', how='left')
    merged_data.to_csv('merged_data.csv', index=False)
```

```python
def clean_csv():
    '''Remove all rows where prediction skipped, do this after merging'''

    df = pd.read_csv('merged_data.csv')
    df_clean = df[~df['predictions'].str.contains('SKIPPED:', case=False, na=False)]
    df_waste = df[df['predictions'].str.contains('SKIPPED:', case=False, na=False)]
    df_clean.to_csv('cleaned_data.csv', index=False)
    df_waste.to_csv('waste_data.csv', index=False)       # waste.csv to gather all proteins that got skipped
```

```python
def basic_stats(row):
    '''Return basic stats per row'''

    try:
        quick_set = set(ast.literal_eval(row["go_terms"]))
        pred_set = set(ast.literal_eval(row["predictions"]))

        return pd.Series({
            "tp": len(quick_set & pred_set),
            "fn": len(quick_set - pred_set),
            "fp": len(pred_set - quick_set),
            "len_quick": len(quick_set),
            "len_pred": len(pred_set)
        })
    except (ValueError, SyntaxError):
        return pd.Series({"tp": 0, "fn": 0, "fp": 0, "len_quick": 0, "len_pred": 0})
```

```python
def basic_stats_namesplit(row):
    '''Return basic stats seperated by namespace '''

    try:

        quick_set = set(ast.literal_eval(row["go_terms"]))
        pred_set = set(ast.literal_eval(row["predictions"]))

        namespaces = {"biological_process": {"tp":0,"fn":0,"fp":0},
                            "molecular_function": {"tp":0,"fn":0,"fp":0},
                            "cellular_component": {"tp":0,"fn":0,"fp":0}}
        
        for term in quick_set & pred_set:
            if term in cache:
                namespaces[cache[term]["namespace"]]["tp"] += 1
        for term in quick_set - pred_set:
            if term in cache:
                namespaces[cache[term]["namespace"]]["fn"] += 1
        for term in pred_set - quick_set:
            if term in cache:
                namespaces[cache[term]["namespace"]]["fp"] += 1

        return pd.Series({f"{ns}_{k}": v
                            for ns, counts in namespaces.items()
                            for k, v in counts.items()})
    except (ValueError, SyntaxError):
        # there might be a better option than pass here
        pass
```

```python
def basic_analysis(data):
    '''Calculate tp, fp, fn, and additional statistics for merged data.
        Shows single raw confusion matrix.'''
    
    merged_data = pd.read_csv(data)
    stats = merged_data.progress_apply(basic_stats, axis=1)
    merged_data = pd.concat([merged_data, stats], axis=1)
    
    merged_data.to_csv('basic_merged_anew.csv', index=False)


    tp = int(merged_data["tp"].sum())
    fp = int(merged_data["fp"].sum())
    fn = int(merged_data["fn"].sum())
    
    y_true = [1]*tp + [1]*fn + [0]*fp
    y_pred = [1]*tp + [0]*fn + [1]*fp

    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    print(f"Precision: {precision:.4f}")
    print(f"Recall: {recall:.4f}")
    print(f"f1: {f1:.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(15, 8))
    fig.suptitle("AlphaFunctor Predictions", fontsize=14)

    con_raw = confusion_matrix(y_true, y_pred)
    con_disp_raw = ConfusionMatrixDisplay(con_raw, display_labels=["Negative", "Positive"])
    con_disp_raw.plot(ax=axes[0], values_format='d')
    axes[0].set_title("Raw Counts")
    axes[0].set_xlabel("AlphaFunctor Prediction")
    axes[0].set_ylabel("UniProt Annotation")

    con_norm = confusion_matrix(y_true, y_pred, normalize='all')
    con_disp_norm = ConfusionMatrixDisplay(con_norm, display_labels=["Negative", "Positive"])
    con_disp_norm.plot(ax=axes[1], values_format='.2f')
    axes[1].set_title("Normalized")
    axes[1].set_xlabel("AlphaFunctor Prediction")
    axes[1].set_ylabel("UniProt Annotation")

    plt.figtext(0.5, 0.03, f"Precision: {precision:.4f} | Recall: {recall:.4f} | F1: {f1:.4f}", 
         ha='center', fontsize=11)
    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.subplots_adjust(bottom=0.15)
    plt.savefig("basic_stats_anew.png", bbox_inches="tight")
    plt.show()
```

```python
def basic_name_analysis(data):
    '''Calculate tp, fp, fn, and additional statistics for merged data seperated by namespace.
        Show 3 confusion matrices.'''

    merged_data = pd.read_csv(data)
    stats = merged_data.progress_apply(basic_stats_namesplit, axis=1)
    merged_data = pd.concat([merged_data, stats], axis=1)
    
    merged_data.to_csv('basic_namesplit_anew.csv', index=False)

    namespaces = ["biological_process", "molecular_function", "cellular_component"]
    labels = ["BP", "MF", "CC"]

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("AlphaFunction Predictions by Namespace")

    for i, (namespace, label) in enumerate(zip(namespaces, labels)):
        tp = int(merged_data[f"{namespace}_tp"].sum())
        fp = int(merged_data[f"{namespace}_fp"].sum())
        fn = int(merged_data[f"{namespace}_fn"].sum())

        y_true = [1]*tp + [1]*fn + [0]*fp
        y_pred = [1]*tp + [0]*fn + [1]*fp

        precision = precision_score(y_true, y_pred)
        recall = recall_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred)

        print(f"{namespace} precision: {precision:.4f}")
        print(f"{namespace} recall: {recall:.4f}")
        print(f"{namespace} f1: {f1:.4f}")
        print()

        # raw counts
        con_raw = confusion_matrix(y_true, y_pred)
        con_disp_raw = ConfusionMatrixDisplay(con_raw)
        con_disp_raw.plot(ax=axes[0][i], values_format='d')
        axes[0][i].set_title(f"{label} - Raw Counts")
        axes[0][i].set_xlabel(f"AlphaFunctor Prediction \nP: {precision:.4f} | R: {recall:.4f} | F1: {f1:.4f}")
        axes[0][i].set_ylabel("Uniprot Annotation")
        
        # normalized
        con_norm = confusion_matrix(y_true, y_pred, normalize='all')
        con_disp_norm = ConfusionMatrixDisplay(con_norm)
        con_disp_norm.plot(ax=axes[1][i], values_format='.2f')
        axes[1][i].set_title(f"{label} - Normalized")
        axes[1][i].set_xlabel("AlphaFunctor Prediction")
        axes[1][i].set_ylabel("UniProt Annotation")

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.savefig("basic_namesplit_stats2.png", bbox_inches="tight")
    plt.show()
```