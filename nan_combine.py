import glob
import pandas as pd
import re

'''Combine all the output files from the HPCC
    (optional: run job_combine.sb for this on HPCC then download just the final_output file after)'''

files = sorted(glob.glob("outputs/output_*.csv"), key=lambda f: int(re.search(r'_(\d+)\.csv$', f).group(1)))
df = pd.concat([pd.read_csv(f) for f in files])
df.to_csv("final_output.csv", index=False)