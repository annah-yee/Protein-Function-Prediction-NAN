# In order to retrieve the protein sequences and their uniprot annotations,
# must download all the reviewed SwissProt entries in tsv format

''' This file contains the tools to convert the SwissProt tsv into xxx
    and clean it up a bit'''

import pandas as pd
import ast


def deduplicate_go_terms(val):
    try:
        terms = ast.literal_eval(val)
        seen = set()
        unique = []
        for term in terms:
            key = tuple(term)
            if key not in seen:
                seen.add(key)
                unique.append(term)
        return str(unique)
    except Exception:
        return val 
    
def return_annotated(df):
    mask = (
        df['go_terms'].notna() &
        (df['go_terms'].str.strip() != '[]')
    )
    return df[mask]

def create_chonk():

    # normal tsv parsing and column renaming
    # replace with the name of your file if you want to redo this
    df = pd.read_csv('uniprotkb_AND_reviewed_true_2026_06_29.tsv', sep='\t')
    df["go_terms"] = df["Gene Ontology IDs"].str.split("; ")
    df = df.drop(columns=['Gene Ontology IDs'])
    df = df.rename(columns={'Entry':'uniprot_id', 
                            'Entry Name':'entry_name',
                            'Date of last modification':'date_mod',
                            'Reviewed':'reviewed',
                            'Sequence':'sequence',
                            })
    
    df["go_terms"] = df["go_terms"].apply(deduplicate_go_terms)

    # annotated is just the cleaned up version of the big data, will act as ground truth later
    annotated_df = return_annotated(df)
    annotated_df.to_csv("annotated_swissprotkb.csv", index=False)
    
    # condensed just includes the uniprot id and sequence. this is the model input
    condensed_df = annotated_df[['uniprot_id','sequence']]
    condensed_df.to_csv("condensed_swissprotkb.csv", index=False)


if __name__ == "__main__":
    create_chonk()


# like how streamlined do we want this to be?
# want to be able to do this whole process starting at a new uniprotkb file and getting out
# all the options? oh man
