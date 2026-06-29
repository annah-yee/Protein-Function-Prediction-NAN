import pandas as pd
import ast
import os
import re

# df = pd.read_csv('/home/annaeyee/nan_evaulation/chonker_data2.csv')

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
    
def remove_gonames(val):
    try:
        terms = ast.literal_eval(val)
        return [t[0] for t in terms]
    except Exception:
        return []
    
def evidence_split(filepath):
    df = pd.read_csv(filepath)
    
    # get all unique evidence codes
    evidence_codes = [
        "EXP", "IDA", "IPI", "IMP", "IGI", "IEP",
        "HTP", "HDA", "HMP", "HGI", "HEP",
        "IBA", "IBD", "IKR", "IRD",
        "ISS", "ISO", "ISA", "ISM", "IGC",
        "RCA", "TAS", "NAS", "IC", "ND", "IEA"
    ]
    
    # initialize all evidence code columns as empty lists
    for code in evidence_codes:
        df[code] = None
    
    for idx, row in df.iterrows():
        try:
            data = ast.literal_eval(row["go_terms"])
            splits = {code: [] for code in evidence_codes}
            
            for go_id, e_code in data:
                if e_code in splits:
                    splits[e_code].append(go_id)
            
            for code in evidence_codes:
                splits[code] = list(set(splits[code]))
                df.at[idx, code] = str(splits[code])
                
        except (ValueError, SyntaxError):
            for code in evidence_codes:
                df.at[idx, code] = str([])
    
    df.to_csv("evidence_split.csv", index=False)
    return df

# evidence_split('/home/annaeyee/nan_evaulation/annotated_evidence.csv')

def return_annotated(filepath):
    df = pd.read_csv(filepath)
    mask = (
        df['go_terms'].notna() &
        (df['go_terms'].str.strip() != '[]')
    )
    return df[mask]

# annotated = return_annotated('/home/annaeyee/nan_evaulation/chonker_evidence.csv')
# annotated.to_csv('data3_testagain.csv')

def return_condensed(filepath):
    df = pd.read_csv(filepath)
    new = df[['uniprot_id','sequence']]
    return new

def find_missing():
    # this will hopefully find what outputs we're missing
    files = os.listdir('/tank/anna/Protein-Function-Prediction/outputs')
    numbers = sorted(int(re.search(r'_(\d+)\.csv$', f).group(1)) for f in files if re.search(r'_(\d+)\.csv$', f))

    full_range = set(range(numbers[0], numbers[-1] + 1))
    missing = sorted(full_range - set(numbers))

    print(missing)

# find_missing()

def clean_csv(file):
    # go thru the chonker csv and remove all rows where "skipped is a thing"
    # this should filter out everything that has obsolete goTerms or are a strange length
    # do this after merging 
    df = pd.read_csv(file)
    df_clean = df[~df['predictions'].str.contains('SKIPPED:', case=False)]
    df_waste = df[df['predictions'].str.contains('SKIPPED:', case=False)]
    df_clean.to_csv('cleaned_evidence.csv', index=False)
    df_waste.to_csv('waste_evidence.csv', index=False)
    
def compare(one, two):
    uno = pd.read_csv(one)
    dos = pd.read_csv(two)

    uno_set = set(uno['sequence'].dropna())
    dos_set = set(dos['sequence'].dropna())

    diff = dos_set - uno_set

    print(f"Found {len(diff)} IDs in {two} that are not in {one}:")
    print(diff)

# compare('annotated_evidence.csv', 'final_output.csv')

# evidence history:
# chonker_evidence straight from quickgo
# annotated_evidence
# evidence_split


# df['go_terms'] = df['go_terms'].apply(deduplicate_go_terms)
# output_path = '/home/annaeyee/nan_evaulation/chonker_data2.csv'

# df['go_terms'] = df['go_terms'].apply(remove_gonames)
# output_path = '/home/annaeyee/nan_evaulation/chonker_data3.csv'

# annotated = return_annotated('/home/annaeyee/nan_evaulation/chonker_data3.csv')
# annotated.to_csv('/home/annaeyee/nan_evaulation/chonker_data_annotated.csv', index=False)

# condensed = return_condensed('/home/annaeyee/nan_evaulation/chonker_data_annotated.csv')
# condensed.to_csv('/home/annaeyee/nan_evaulation/chonker_data_condensed.csv', index=False)


# chonker_data.csv is the original big file
# chonker_data2.csv is deduplicated in the last column
# chonker_data3.csv is deduplicated and only goIds in the last column
# chonker_data_annotated.csv is ^^ but only includes the proteins that have go annotations
# chonker_data_condensed.csv is ^^ but stripped to only uniprot id and sequence

# trying to create new chonker data with all the evidence code splits. 
# first grab all the data (with the go_terms column having [[GO0123453, 'IEA'], [GO091309, 'ISS'], ...])
# then maybe return annotated? then evidence split?