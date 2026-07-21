import torch
import os
import warnings
warnings.filterwarnings("ignore", category=UserWarning)
from src import get_data,get_prediction_5fold,check_input
import pandas as pd
import time
from tqdm import tqdm

'''Meant to be ran by big_job.sb eventually running all proteins through in chunks of 90'''
'''Essentially a tweaked version of the original prediction.py to work for the HPCC'''

tqdm.pandas()

task_id = int(os.environ.get('SLURM_ARRAY_TASK_ID', 0))
job_id = int(os.environ.get('SLURM_ARRAY_JOB_ID',0))

def add_prediction(sequence):
    ''' Return string of GoName predictions.'''
    try:

        seq = check_input(sequence)
        folder = f'./temp-data/task_{task_id}/'
        os.makedirs(folder, exist_ok=True)
        dataloader = get_data(seq, folder)
        model_folder = './src/checkpoint/'
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        prediction = get_prediction_5fold(model_folder,dataloader,device)
        return prediction
    except (ValueError, KeyError, TypeError) as e:
        print(f"Skipping sequence {sequence[:20]}...: {e}")
        return f"SKIPPED: {e}"


def nan_main():
    ''' Generate predictions for a selection of protein sequences '''
    begin = time.perf_counter()
    chunk_size = 90

    df = pd.read_csv('condensed_swissprotkb.csv')
    start = task_id * chunk_size
    end = start + chunk_size
    chunk = df.iloc[start:end].copy()

    chunk['predictions'] = chunk['sequence'].progress_apply(add_prediction)

    os.makedirs("outputs", exist_ok=True)
    chunk.to_csv(f"outputs/output_{job_id}_{task_id}.csv", index=False)

    duration = time.perf_counter() - begin
    print(f"Task {task_id} completed in {duration:.2f}s")

if __name__ == "__main__":
    nan_main()

# just wanted to add somethin somethin