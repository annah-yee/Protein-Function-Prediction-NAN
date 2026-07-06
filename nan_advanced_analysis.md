# Advanced Analysis

To run basic analysis, must have the `final_output.csv` (aka all the predictions collected from Alphafunctor) and the annotations collected from UniProt `annotated_swissprotkb.csv`.

more explanation here too


```python
def get_relatives(go_id):
    '''Get sets of parents and children for a go_id'''

    if go_id not in cache:
        return set(), set()
    term = cache[go_id]
    parents = set(term['parents'])
    children = set(term['children'])
    return parents, children


def classify_missing(go_term, comparison_set):
    '''Intake an unmatched go term compare it to the set it is mising from to find possible relatives'''

    if go_term not in cache:
        return "unknown"
    
    parents, children = get_relatives(go_term)

    if (parents & comparison_set) and (children & comparison_set):
        return "both_found"

    elif children & comparison_set:
        return "child_found"
    
    elif parents & comparison_set:
        return "parent_found"

    else:
        return "unrelated"
```

```python
def advanced_stats(row):
    '''Intake row, compare quick_go and predictions, return series of stats'''

    quick_set = set(ast.literal_eval(row["go_terms"]))
    pred_set = set(ast.literal_eval(row["predictions"]))

    matches = list(quick_set.intersection(pred_set))
    missing_in_quick = list(pred_set.difference(quick_set))
    missing_in_pred = list(quick_set.difference(pred_set))

    stats = {
        "true_positives": len(matches),
        "missing_in_quick": len(missing_in_quick),
        "missing_in_pred": len(missing_in_pred),

        # in quick but missed by pred
        "missed_sandwich": 0,   # prediction missed a sandwich term (how is this possible w/ true path?)
        "missed_parent": 0,     # prediction missed a parent term (this is ok)
        "missed_child": 0,      # prediction missed a child term (not going deep enough)
        "missed_unrelated": 0,  # prediction missed an unrelated term (wiff)
        "missed_unknown": 0,    # term isn't in the cache (obselete?)

        # in pred but not in quick
        "pred_sandwich": 0,     # prediction guessed a sandwich term (this is ok)
        "pred_parent": 0,       # prediction guessed a parent term (common)
        "pred_child": 0,        # prediction guessed a child term (digging deeper!)
        "pred_unrelated": 0,    # prediction guessed an unrelated term (EXCITING!)
        "pred_unknown": 0       # term isn't in the cache (obselete?)
        }

    for term in missing_in_pred:
        result = classify_missing(term, pred_set)
        if result == "both_found":
            stats["missed_sandwich"]+=1

        elif result == "child_found":
            stats["missed_parent"]+=1
            
        elif result == "parent_found":
            stats["missed_child"]+=1

        elif result == "unrelated":
            stats["missed_unrelated"]+=1

        elif result == "unknown":
            stats["missed_unknown"]+=1

    for term in missing_in_quick:
        result = classify_missing(term, quick_set)
        if result == "both_found":
            stats["pred_sandwich"]+=1

        elif result == "child_found":
            stats["pred_parent"]+=1

        elif result == "parent_found":
            stats["pred_child"]+=1

        elif result == "unrelated":
            stats["pred_unrelated"]+=1

        elif result == "unknown":
            stats["pred_unknown"]+=1

    return pd.Series(stats)
```

```python
def advanced_analysis(data, output_png="advanced.png"):
    '''Complete analysis with hierarchical context and produce graps'''

    # input merged data, extract advanced stats and add to df
    merged_data = pd.read_csv(data).dropna(subset=["go_terms", "predictions"])
    stats = merged_data.progress_apply(advanced_stats, axis=1)
    merged_data = pd.concat([merged_data, stats], axis=1)

    merged_data.to_csv('anew_advanced_stats.csv', index=False)

    fig = plt.figure(figsize=(20, 5))
    gs = gridspec.GridSpec(1, 4) 

    ax_cmr = fig.add_subplot(gs[0, 0])   # raw confusion matrix
    ax_cmn = fig.add_subplot(gs[0,1])    # normalized confusion matrix
    ax_bar1 = fig.add_subplot(gs[0, 2])  # false positive bar chart
    ax_bar2 = fig.add_subplot(gs[0, 3])  # false negative bar chart

    # making the general matrix 

    tp = int(merged_data["true_positives"].sum())
    fp = int(merged_data["missing_in_quick"].sum())
    fn = int(merged_data["missing_in_pred"].sum())

    y_true = [1]*tp + [1]*fn + [0]*fp
    y_pred = [1]*tp + [0]*fn + [1]*fp

    # raw confusion matrix
    gen_raw = confusion_matrix(y_true, y_pred)
    gen_disp_raw = ConfusionMatrixDisplay(gen_raw, display_labels=["Negative", "Positive"])
    gen_disp_raw.plot(ax=ax_cmr, values_format='d')
    ax_cmr.set_title("Raw Counts")
    ax_cmr.set_xlabel("AlphaFunctor Prediction")
    ax_cmr.set_ylabel("UniProt Annotation")

    # normalized confusion matrix
    gen_norm = confusion_matrix(y_true, y_pred, normalize='all')
    gen_disp_norm = ConfusionMatrixDisplay(gen_norm, display_labels=["Negative", "Positive"])
    gen_disp_norm.plot(ax=ax_cmn, values_format='.2f')
    ax_cmn.set_title("Normalized")
    ax_cmn.set_xlabel("AlphaFunctor Prediction")
    ax_cmn.set_ylabel("UniProt Annotation")

    # making the FP bar chart
    # FP = these are terms that AlphaFunctor predicted that were not included in quick go
    fp_categories = ["Sandwich", "Parent", "Child", "Unrelated"]
    fp_values = [
        int(merged_data["pred_sandwich"].sum()),    # quick didn't include a sandwich term
        int(merged_data["pred_parent"].sum()),      # quick didn't include all of the parent terms
        int(merged_data["pred_child"].sum()),       # quick didn't include a child term (interesting!)
        int(merged_data["pred_unrelated"].sum())    # quick didn't include any related terms
    ]
    bars1 = ax_bar1.bar(fp_categories, fp_values)
    ax_bar1.set_title("False Positive Breakdown")
    ax_bar1.set_xlabel("Category")
    ax_bar1.set_ylabel("Count")
    for bar in bars1:
        ax_bar1.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f'{int(bar.get_height())}', ha='center', va='bottom')

    # making the FN bar chart
    # FN = these are terms that quick go included but AlphaFunctor missed
    fn_categories = ["Sandwich", "Parent", "Child", "Unrelated"]
    fn_values = [
        int(merged_data["missed_sandwich"].sum()),  # prediction had parent and child of term
        int(merged_data["missed_parent"].sum()),    # prediction had child of term
        int(merged_data["missed_child"].sum()),     # prediction had parent of term
        int(merged_data["missed_unrelated"].sum())  # prediction didn't have any relatives of term
    ]
    bars2 = ax_bar2.bar(fn_categories, fn_values)
    ax_bar2.set_title("False Negative Breakdown")
    ax_bar2.set_xlabel("Category")
    ax_bar2.set_ylabel("Count")
    for bar in bars2:
        ax_bar2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    f'{int(bar.get_height())}', ha='center', va='bottom')

    # plot and save
    plt.tight_layout()
    plt.savefig(output_png)
```

```python
def create_namesplits(data):
    '''Create new dataframes with go_terms and predictions filtered by namespace'''

    df = pd.read_csv(data).dropna(subset=["go_terms", "predictions"])

    bp_df = df.copy()
    mf_df = df.copy()
    cc_df = df.copy()
    
    for idx, row in df.iterrows():
        quick_set = set(ast.literal_eval(row["go_terms"]))
        pred_set = set(ast.literal_eval(row["predictions"]))
        
        # filter quick terms by namespace
        bp_quick = [t for t in quick_set if t in cache and cache[t]["namespace"] == "biological_process"]
        mf_quick = [t for t in quick_set if t in cache and cache[t]["namespace"] == "molecular_function"]
        cc_quick = [t for t in quick_set if t in cache and cache[t]["namespace"] == "cellular_component"]
        
        # filter pred terms by namespace
        bp_pred = [t for t in pred_set if t in cache and cache[t]["namespace"] == "biological_process"]
        mf_pred = [t for t in pred_set if t in cache and cache[t]["namespace"] == "molecular_function"]
        cc_pred = [t for t in pred_set if t in cache and cache[t]["namespace"] == "cellular_component"]
        
        bp_df.at[idx, "go_terms"] = str(bp_quick)
        bp_df.at[idx, "predictions"] = str(bp_pred)
        
        mf_df.at[idx, "go_terms"] = str(mf_quick)
        mf_df.at[idx, "predictions"] = str(mf_pred)
        
        cc_df.at[idx, "go_terms"] = str(cc_quick)
        cc_df.at[idx, "predictions"] = str(cc_pred)
    
    bp_df.to_csv("anew_bp_split.csv", index=False)
    mf_df.to_csv("anew_mf_split.csv", index=False)
    cc_df.to_csv("anew_cc_split.csv", index=False)
```

```python
def advanced_analysis_namesplit():
    '''Complete advanced analysis on each namespace split data and generate graphs'''

    # create_namesplits("merged_anew.csv")
    advanced_analysis('anew_bp_split.csv', 'anew_bp_split.png')
    advanced_analysis('anew_mf_split.csv', 'anew_mf_split.png')
    advanced_analysis('anew_cc_split.csv', 'anew_cc_split.png')
```
