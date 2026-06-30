# Getting SwissProt Data
## Download tsv file
Before all of this, you have to make sure you have the swiss prot data downloaded
1. go to https://www.uniprot.org/uniprotkb?facets=reviewed%3Atrue&query=* 
2. customize columns to include entry name, entry, sequence, gene ontology IDs, and data of last modification
3. download all as a TSV file, not compressed

## Import requirements

```python
import pandas as pd
import ast
```
## Generate csv files

Using the code below, insert the name of your tsv file and run `create_chonk()`. Two csv files are generated: annotated and codensed versions of the data. 

```python
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
    df = pd.read_csv(<YOUR_TSV>, sep='\t')    # replace with name of your tsv file
    df["go_terms"] = df["Gene Ontology IDs"].str.split("; ")
    df = df.drop(columns=['Gene Ontology IDs'])
    df = df.rename(columns={'Entry':'uniprot_id', 
                            'Entry Name':'entry_name',
                            'Date of last modification':'date_mod',
                            'Sequence':'sequence',
                            })
    
    df["go_terms"] = df["go_terms"].apply(deduplicate_go_terms)

    # annotated is just the cleaned up version of the big data, will act as ground truth later
    annotated_df = return_annotated(df)
    annotated_df.to_csv("annotated_swissprotkb.csv", index=False)
    
    # condensed just includes the uniprot id and sequence. this is the model input
    condensed_df = annotated_df[['uniprot_id','sequence']]
    condensed_df.to_csv("condensed_swissprotkb.csv", index=False)
```

## Where to go from here? 
Using the condensed csv file, start running those bad boys through the hpcc #insertlink?