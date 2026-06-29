import pandas as pd

def create_chonk():

    # normal tsv parsing and column renaming
    df = pd.read_csv('uniprotkb_anew.tsv', sep='\t')
    df["go_terms"] = df["Gene Ontology IDs"].str.split("; ")
    df = df.drop(columns=['Gene Ontology IDs'])
    df = df.rename(columns={'Entry':'uniprot_id', 
                            'Entry Name':'entry_name',
                            'Date of last modification':'date_mod',
                            'Reviewed':'reviewed',
                            'Sequence':'sequence',
                            })
    df.to_csv('anew_chonker.csv', index=False)

create_chonk()
