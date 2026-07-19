# Collecting Predictions from AlphaFunctor
To effeciently push over 500k proteins through AlphaFunctor and collect the outputs, we seperate the data into chunks of 90 proteins each and submitted chunks to the HPCC incrementally. At around 120s/it of AlphaFunctor, each chunk can be submitted with 4 hours wall time. 

Examples of my slurm scripts are included under `prediction_job.sb` and `combine_job`. 

Make sure the `condensed_swissprotkb.csv`, `prediction_job.sb`, `combine_job.sb`, `nan_prediction.py`, and `nan_combine.py` are in the HPCC working directory. 

## Submitting jobs
Submit with
```bash
sbatch prediction_job.sb
```

When I ran this, I had around 6150 tasks with an HPCC limit of 1000 tasks per user. This means my initial submission was the task array 0 to 999. About twice a day, I ran these these commands to check on the progress.

```bash
squeue -u <YOUR_USERNAME> -h -t r -r | wc -l 
```

```bash
squeue -u <YOUR_USERNAME> -h -t pd -r | wc -l 
```

Using running and pending counts, we can calculate how many more tasks to submit. Make sure to keep track of which task ID you're on!

Here's a snippet of my notes that kept me organized.

```text

    Submitted tasks 2101-2800 under job 11409821 

    Submitted tsks 2801-3200 under job 11448236 

Submitted tasks 888,891,900,908,922,923,962,986,1032,1034,1050,1072,1115,1128,1131,1135,1138,1141,1142,1160,1177 under job 11489469 with more time 

    Submitted tasks 3201-3700 under job 11489583 

Submitted tasks 2303,2314,2321,2330,2331,2337,2456,2477,2549,2593,2725 under job 11489635 with more time  

Submitted tasks 379,462,526,703,707,1981 under job 11490830 with more memory
```

```bash
sacct | grep TIMEOUT | awk '{print $1}' | cut -d'_' -f2 | paste -sd, - 
```

Some tasks may take longer than 4 hours, this will output a comma seperated list of the task numbers that have timed out. 
Copy/paste the list and resubmit the tasks with:
Note: If a task times out again, you'll need to retry with more time

```bash
sbatch --array=<PASTE_TASK_NUMBERS> --time=05:00:00 prediction_job.sb
```

Some tasks may run out of memory, find them with this commant (ex. start date 06/01)

```bash
sacct -S <START_DATE> --noheader --format=JobID%30,State%15 -X | grep OUT_OF_MEM
```
and resubmit with

```bash
sbatch --array=<TASK_NUMBERS> --mem=32G prediction_job.sb
```

Even with all of this monitoring, a task or two may slip through the cracks. There's a tool called `find_missing()` in `nan_tools.py` that can help identify missing and/or duplicate outputs. 

In the event where tasks were accidentally ran twice, manually compare the outputs or delete both files and submit the task again. 

## Combining outputs

Once you've confirmed all of the outputs have been generated, combine with 

```bash
sbatch combine_job.sb
```
