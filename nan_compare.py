import pickle
from goatools.obo_parser import GODag
import pandas as pd
import ast
import sklearn as sk
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay, recall_score,precision_score, f1_score
import matplotlib.pyplot as plt
import os
import re
from tqdm import tqdm
import matplotlib.gridspec as gridspec
from nan_tools import plot_confusion


'''Global things'''
tqdm.pandas()       # allow for progress tracking

with open("go_cache.pkl", "rb") as f:       # holds information about each go_term
        cache = pickle.load(f)

'''Some tools to work with the data'''

def merge(quick_data, prediction_data):
    '''Merges data from AlphaFunctor outputs and UniProt annotations'''

    main_data = pd.read_csv(quick_data, sep=',')
    predictions = pd.read_csv(prediction_data, sep=',')
    merged_data = pd.merge(main_data, predictions[['uniprot_id', 'predictions']], on='uniprot_id', how='left')
    merged_data.to_csv('merged_data_quick.csv', index=False)


def clean_csv(file):
    '''Remove all rows where prediction skipped, do this after merging'''

    df = pd.read_csv(file)
    df_clean = df[~df['predictions'].str.contains('SKIPPED:', case=False, na=False)]
    df_waste = df[df['predictions'].str.contains('SKIPPED:', case=False, na=False)]
    df_clean.to_csv('cleaned_data_quick.csv', index=False)
    df_waste.to_csv('waste_data_quick.csv', index=False)       # waste.csv to gather all proteins that got skipped

# merge the annoated kb and the alphafunctor outputs, then clean the data
# merge("quickgo.csv", "alphafunctor_output.csv")
# clean_csv('merged_data_quick.csv')

'''Basic analysis things'''


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
    
# def basic_evidence_stats(row, evidence_codes):
#     '''Return basic stats per row, only considering select evidence codes'''

#     try:
#         quick_set = set()
#         for code in evidence_codes:
#             add = ast.literal_eval(row[code])
#             quick_set.update(add)

#         pred_set = set(ast.literal_eval(row["predictions"]))

#         return pd.Series({
#             "tp": len(quick_set & pred_set),
#             "fn": len(quick_set - pred_set),
#             "fp": len(pred_set - quick_set),
#             "len_quick": len(quick_set),
#             "len_pred": len(pred_set)
#         })
#     except (ValueError, SyntaxError):
#         return pd.Series({"tp": 0, "fn": 0, "fp": 0, "len_quick": 0, "len_pred": 0})
    

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

def basic_analysis(data, plot):
    '''Calculate tp, fp, fn, and additional statistics for merged data.
        Shows single raw confusion matrix.'''
    
    merged_data = pd.read_csv(data)
    stats = merged_data.progress_apply(basic_stats, axis=1)
    merged_data = pd.concat([merged_data, stats], axis=1)
    
    merged_data.to_csv('basic_stats_data_quick.csv', index=False)

    tp = int(merged_data["tp"].sum())
    fp = int(merged_data["fp"].sum())
    fn = int(merged_data["fn"].sum())
    
    y_true = [1]*tp + [1]*fn + [0]*fp
    y_pred = [1]*tp + [0]*fn + [1]*fp

    precision = precision_score(y_true, y_pred)
    recall = recall_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred)

    metrics = {"tp": tp, 
                "fp": fp,
                "fn": fn,
                "precision": precision,
                "recall": recall,
                "f1": f1}

    if plot:
        plot_confusion(metrics, comparison = "QuickGo", ax=None) # may change to QuickGo

# basic_analysis('cleaned_data_quick.csv', plot=True)

def basic_name_analysis(data, plot):
    '''Calculate tp, fp, fn, and additional statistics for merged data seperated by namespace.
        Show 3 confusion matrices.'''

    merged_data = pd.read_csv(data)
    stats = merged_data.progress_apply(basic_stats_namesplit, axis=1)
    merged_data = pd.concat([merged_data, stats], axis=1)
    
    # merged_data.to_csv('basic_namesplit_data.csv', index=False)

    namespaces = ["biological_process", "molecular_function", "cellular_component"]
    labels = ["BP", "MF", "CC"]

    metrics = {}

    for (namespace, label) in zip(namespaces, labels):
        tp = int(merged_data[f"{namespace}_tp"].sum())
        fp = int(merged_data[f"{namespace}_fp"].sum())
        fn = int(merged_data[f"{namespace}_fn"].sum())

        y_true = [1]*tp + [1]*fn + [0]*fp
        y_pred = [1]*tp + [0]*fn + [1]*fp

        precision = precision_score(y_true, y_pred)
        recall = recall_score(y_true, y_pred)
        f1 = f1_score(y_true, y_pred)

        metrics[label] = {"tp": tp, 
                "fp": fp,
                "fn": fn,
                "precision": precision,
                "recall": recall,
                "f1": f1}
        
    if plot: 
        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        plot_confusion(metrics["BP"], "BP UniProt", ax=axes[0])
        plot_confusion(metrics["MF"], "MF UniProt", ax=axes[1])
        plot_confusion(metrics["CC"], "CC UniProt", ax=axes[2])
        plt.tight_layout()
        plt.savefig("namesplit_confusion.png", bbox_inches="tight")
        plt.savefig("namesplit_confusion.pdf", bbox_inches="tight")
        plt.show()


# basic_name_analysis('cleaned_data.csv', plot=True)

def basic_evidence_analysis(data):
    # would be deciding certain select set of evidence codes to focus on and creating different plots for each

    # first change up the df depending on what we want our go_terms column to include
    # then run basic stats, sum tp tn fp, calculate the metrics, plot the plots
    # do we super care about namesplits here? maybe only for an evidence split that we loVE yk

    # so would want to do some 2x3 tight layout of raw and normalized confusion matrices for these


    # want to go one by one for the plots and then plot at the end

    merged_data = pd.read

    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    fig.suptitle("AlphaFunction Predictions with Evidence Splits")

    # Everything first
    evidence_codes = [
        "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
        "HTP", "HDA", "HMP", "HGI", "HEP",
        "IBA", "IBD", "IKR", "IRD",
        "ISS", "ISO", "ISA", "ISM", "IGC",
        "RCA", "TAS", "NAS", "IC", "ND", "IEA"
    ]

    # need to do basic stats and then add metrics/raw/normalized plots to the fig

    # Removing electronic annotation
    manual_evidence = evidence_codes.copy()
    manual_evidence.remove("IEA")

    # need to do basic stats and then add metrics/raw/normalized plots to the fig

    # Removing author statement evidence, curatorial and inferred evidence
    EHC_evidence = manual_evidence.copy()
    EHC_evidence.remove("TAS", "NAS", "IC", "ND")

    # need to do basic stats and then add metrics/raw/normalized plots to the fig

    plt.tight_layout(rect=[0, 0.08, 1, 0.95])
    plt.savefig("basic_namesplit_stats2.png", bbox_inches="tight")
    plt.show()


'''Advanced analysis things (considering hierarchical relationships)'''

def create_DAG():
    ''' creates pickle object of a dict of dictionaries with all the important information (no goNames tho)''' 

    godag = GODag("go-basic.obo")
    cache = {}
    for go_id, term in godag.items():
        cache[go_id] = {
            "name": term.name,
            "namespace": term.namespace,
            "parents": term.get_all_parents(),
            "children": term.get_all_children(),
            "depth": term.depth,
            "level": term.level
        }
    with open("go_cache.pkl", "wb") as f:
        pickle.dump(cache, f)

# create_DAG()


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

# def selective_tp():
#     # return some re-imagined tp

# def selective_fp():
#     # ^^ 

# def selective_fn():
#     # ^^   

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

# create_namesplits("cleaned_anew.csv")


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

plt.rcParams.update({
    'font.size': 20,
    'axes.titlesize': 20,
    'axes.labelsize': 20,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
})

def advanced_analysis(data, output_png="advanced_stats_data.png"):
    '''Complete analysis with hierarchical context and produce graps'''

    # input merged data, extract advanced stats and add to df
    merged_data = pd.read_csv(data).dropna(subset=["go_terms", "predictions"])
    stats = merged_data.progress_apply(advanced_stats, axis=1)
    merged_data = pd.concat([merged_data, stats], axis=1)

    merged_data.to_csv('advanced_stats_quick.csv', index=False)

    # merged_data = pd.read_csv("advanced_stats_data.csv")

    fig = plt.figure(figsize=(10, 5))
    gs = gridspec.GridSpec(1, 2) 

    # ax_cmr = fig.add_subplot(gs[0, 0])   # raw confusion matrix
    # ax_cmn = fig.add_subplot(gs[0,1])    # normalized confusion matrix
    # ax_bar1 = fig.add_subplot(gs[0, 2])  # false positive bar chart
    # ax_bar2 = fig.add_subplot(gs[0, 3])
    
    ax_bar1 = fig.add_subplot(gs[0, 0]) 
    ax_bar2 = fig.add_subplot(gs[0, 1])  

    # making the general matrix 

    tp = int(merged_data["true_positives"].sum())
    fp = int(merged_data["missing_in_quick"].sum())
    fn = int(merged_data["missing_in_pred"].sum())

    y_true = [1]*tp + [1]*fn + [0]*fp
    y_pred = [1]*tp + [0]*fn + [1]*fp

    # raw confusion matrix
    # gen_raw = confusion_matrix(y_true, y_pred)
    # gen_disp_raw = ConfusionMatrixDisplay(gen_raw, display_labels=["Negative", "Positive"])
    # gen_disp_raw.plot(ax=ax_cmr, values_format='d')
    # ax_cmr.set_title("Raw Counts")
    # ax_cmr.set_xlabel("AlphaFunctor Prediction")
    # ax_cmr.set_ylabel("UniProt Annotation")

    # # normalized confusion matrix
    # gen_norm = confusion_matrix(y_true, y_pred, normalize='all')
    # gen_disp_norm = ConfusionMatrixDisplay(gen_norm, display_labels=["Negative", "Positive"])
    # gen_disp_norm.plot(ax=ax_cmn, values_format='.2f')
    # ax_cmn.set_title("Normalized")
    # ax_cmn.set_xlabel("AlphaFunctor Prediction")
    # ax_cmn.set_ylabel("UniProt Annotation")

    # # making the FP bar chart
    # # FP = these are terms that AlphaFunctor predicted that were not included in quick go
    fp_categories = ["Sandwich", "Parent", "Child", "Unrelated"]
    fp_values = [
        int(merged_data["pred_sandwich"].sum()),    # quick didn't include a sandwich term
        int(merged_data["pred_parent"].sum()),      # quick didn't include all of the parent terms
        int(merged_data["pred_child"].sum()),       # quick didn't include a child term (interesting!)
        int(merged_data["pred_unrelated"].sum())    # quick didn't include any related terms
    ]
    # bars1 = ax_bar1.bar(fp_categories, fp_values)
    # ax_bar1.set_title("False Positive Breakdown")
    # ax_bar1.set_xlabel("Category")
    # ax_bar1.set_ylabel("Count")
    # for bar in bars1:
    #     ax_bar1.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
    #                 f'{int(bar.get_height())}', ha='center', va='bottom')

    ax_bar1.pie(fp_values, labels=fp_categories, autopct='%1.1f%%', colors=["#c9e6b0","#f7d39d", "#b0d4f2", "#f2b6b0"])
    ax_bar1.set_title("False Positive Breakdown")


    # # making the FN bar chart
    # # FN = these are terms that quick go included but AlphaFunctor missed
    fn_categories = ["Sandwich", "Parent", "Child", "Unrelated"]
    fn_values = [
        int(merged_data["missed_sandwich"].sum()),  # prediction had parent and child of term
        int(merged_data["missed_parent"].sum()),    # prediction had child of term
        int(merged_data["missed_child"].sum()),     # prediction had parent of term
        int(merged_data["missed_unrelated"].sum())  # prediction didn't have any relatives of term
    ]
    # bars2 = ax_bar2.bar(fn_categories, fn_values)
    # ax_bar2.set_title("False Negative Breakdown")
    # ax_bar2.set_xlabel("Category")
    # ax_bar2.set_ylabel("Count")
    # for bar in bars2:
    #     ax_bar2.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
    #                 f'{int(bar.get_height())}', ha='center', va='bottom')

    ax_bar2.pie(fn_values, labels=fn_categories, autopct='%1.1f%%', colors=["#c9e6b0","#f7d39d", "#b0d4f2", "#f2b6b0"])
    ax_bar2.set_title("False Negative Breakdown")

    # plot and save
    plt.tight_layout()
    plt.savefig(output_png)

advanced_analysis("cleaned_data_quick.csv", output_png="advanced_stats_data_quick_big.png")

def advanced_analysis_namesplit():
    '''Complete advanced analysis on each namespace split data and generate graphs'''

    # create_namesplits("merged_anew.csv")
    advanced_analysis('anew_bp_split.csv', 'anew_bp_split.png')
    advanced_analysis('anew_mf_split.csv', 'anew_mf_split.png')
    advanced_analysis('anew_cc_split.csv', 'anew_cc_split.png')


# advanced_analysis_namesplit()