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

Wanna grab the quick go stuff too? 

first import these requirements
``` python
import time
import os
import requests
from datetime import datetime
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
```
```python
def GO_fetch(ids, batch_size):
    session = requests.Session()
    retry = Retry(
            total=5,
            backoff_factor=1,
            status_forcelist=[429, 500, 502, 503, 504],
            allowed_methods=["GET"]
        )
    adapter = HTTPAdapter(max_retries=retry)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
        
    session.headers.update({"User-Agent": "GOFetchScript/1.0"})

    results = {}

    for pid in ids:
        results[pid] = []
    batches = [ids[i:i+batch_size] for i in range(0, len(ids), batch_size)]
    url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    headers = {"Accept": "application/json"}

    for i, batch in enumerate(batches, 1):
        batch_ids = ",".join([f'UniProtKB:{id}' for id in batch])
        params = {
            "geneProductId": batch_ids,
            "limit": 100,
            "page": 1,
        }

        while True:
            try:
                response = session.get(
                                        url,
                                        params=params,
                                        headers=headers,
                                        timeout=30)
            except requests.exceptions.ConnectionError as e:
                print(f"Connection error on batch {i}, page {params['page']}: {e}. Retrying in 10s...")
                time.sleep(10)
                continue
            data = response.json()
            for result in data.get("results", []):
                pid = result["geneProductId"].split(":")[1]
                goId = result.get("goId")
                results.setdefault(pid, []).append(goId)
            page_info = data.get("pageInfo", {})
            if params["page"] >= page_info.get("total", 1):
                break
            params["page"] += 1

        time.sleep(0.5)
        print(f'Completed batch {i}/{len(batches)} <3')

    session.close()

    return results           

def create_quick():
    print("And so it begins . . .")
    df = pd.read_csv("condensed_swissprotkb.csv")

    ids = list(df["uniprot_id"])

    ball = GO_fetch(ids, batch_size=240)
    
    df["go_terms"] = df["uniprot_id"].map(ball)
    
    df.to_csv("quickgo.csv", index=False)
```
## Where to go from here? 
Using the condensed csv file, start running those bad boys through the hpcc #insertlink?