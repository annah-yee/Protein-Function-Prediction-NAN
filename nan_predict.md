# Collecting Predictions from AlphaFunctor
To effeciently push over 500k proteins through AlphaFunctor and collect the outputs, we seperate the data into chunks of 90 proteins each and submitted chunks to the HPCC incrementally. At around 120s/it of AlphaFunctor, each chunk can be submitted with 4 hours wall time. 

Examples of my slurm scripts are included under `prediction_job.sb` and `combine_job`. 

Make sure the `condensed_swissprotkb.csv`, `prediction_job.sb`, `combine_job.sb`, `nan_prediction.py`, and `nan_combine.py` are in the HPCC working directory. 

## Submitting jobs
When I ran this, I had around 6150 tasks with an HPCC limit of 1000 tasks per user. This means my initial submission was the task array 0 to 999. About twice a day, I ran these these commands to check on the progress.

```bash
squeue -u <YOUR_USERNAME> -h -t r -r | wc -l 
```

```bash
squeue -u <YOUR_USERNAME> -h -t pd -r | wc -l 
```

Using running and pending counts, we can calculate how many more tasks to submit. Make sure to keep track of which task ID you're on!

```bash
sacct | grep TIMEOUT | awk '{print $1}' | cut -d'_' -f2 | paste -sd, - 
```

Some tasks may take longer than 4 hours, this will output a comma seperated list of the task numbers that have timed out. 
Copy/paste the list and resubmit the tasks with:
Note: If a task times out again, you'll need to retry with more time

```bash
sbatch --array=<PASTE_TASK_NUMBERS> --time=05:00:00 prediction_job.sb
```

Some tasks may run out of memory. Find them by doing:

```bash
i hope i have somehtnig here
```
and resubmit with

```bash
sbatch --array=<PASTE_TASK_NUMBERS> --mem=32G prediction_job.sb
```

Even with all of these tools, a task or two may slip through the cracks. Make sure `check_job.sb` and `nan_check.py` are in your directory and run 

```bash
sbatch check_job.sb
```

then check the output text file for the task numbers that may be missing. Either find the problem in the slurm_log of the task or just resubmit all of the problem tasks with both more time and memory. 

In the event where tasks were accidentally ran twice, manually compare the outputs or delete both files and submit the task again. 

