from Bio import SeqIO
import json
import pandas as pd
import time
import os
import requests
from datetime import datetime
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


df = pd.DataFrame([
    {
        "id": r.id,
        "description": r.description,
        "sequence": str(r.seq),
        "length": len(r.seq)
    }
    for r in SeqIO.parse("uniprot_sprot.fasta", "fasta")
]) 

df["uniprot_id"] = df["id"].str.split("|").str[1]

ids = list(df["uniprot_id"])

batch_size = 240

def parse_obo(filepath):
    go_names = {}
    with open(filepath) as f:
        current_id = None
        for line in f:
            line = line.strip()
            if line.startswith("id: GO:"):
                current_id = line.split("id: ")[1]
            elif line.startswith("name: ") and current_id:
                go_names[current_id] = line.split("name: ")[1]
                current_id = None
    return go_names

if os.path.exists("go_names.json"):
    with open("go_names.json") as f:
        go_names = json.load(f)
else:
    go_names = parse_obo("go.obo")
    with open("go_names.json", "w") as f:
        json.dump(go_names, f)
    print("GoNames dictionary generated!")

def add_terms(go_id):
    if go_id not in go_names:
        print(f"Warning: {go_id} not in OBO file")
        return None
    return go_names[go_id]

def GO_fetch(ids, checkpoint_every=10):
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
        
    # Add a User-Agent (EBI sometimes drops headless requests)
    session.headers.update({"User-Agent": "GOFetchScript/1.0"})

    results = {}
    batch_rows = []
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
                goName = add_terms(goId)
                goEvidence = result.get("goEvidence")
                results.setdefault(pid, []).append([goId, goEvidence]) # only adding the goID and the goEvidence
                batch_rows.append({"uniprot_id": pid, 
                                   "go_id": goId,
                                   "go_name": goName,
                                   "goEvidence": goEvidence
                                   })
            page_info = data.get("pageInfo", {})
            if params["page"] >= page_info.get("total", 1):
                break
            params["page"] += 1

        time.sleep(0.5)
        print(f'Completed batch {i}/{len(batches)} <3')

        if i % checkpoint_every == 0:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            checkpoint_file = (
                f"/home/annaeyee/nan_evaulation/checkpoints/"
                f"checkpoint_{timestamp}_{i}.csv"
            )

            pd.DataFrame(batch_rows).to_csv(checkpoint_file, index=False)
            print(f"Checkpoint saved at batch {i}")

            batch_rows = []

    if batch_rows:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        final_batch_file = (f"/home/annaeyee/nan_evaulation/checkpoints/checkpoint_{timestamp}_final_batch.csv")
        pd.DataFrame(batch_rows).to_csv(final_batch_file, index=False)

    session.close()

    return results           

def main():
    print("And so it begins . . .")
    ball = GO_fetch(ids)
    
    df["go_terms"] = df["uniprot_id"].map(ball)
    
    df.to_csv("chonker_evidence.csv", index=False)
    
if __name__ == "__main__":
    main()