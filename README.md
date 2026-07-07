This repository contains helper functions to perform [UNRAVEL](https://github.com/DelinteNicolas/UNRAVEL) analysis on a [Elikopy](https://github.com/Hyedryn/elikopy) study.

![GitHub repo size](https://img.shields.io/github/repo-size/DelinteNicolas/Tract_Study)

## Setting up

1. Clone the repository
   ```
   git clone https://github.com/DelinteNicolas/Tract_Study.git
   ```
2. Install the required python packages specified in ```requirements.txt```. For instance,

   ```
   pip install unravel-python
   pip install pilab-regis
   pip install xlsxwriter
   ```
4. Modify the root, subject list and mail variables in ```submit_tract_analysis.sh```

## Launching the job

1. Move the working directory inside the folder ```cd Tract_Study```
2. Launch the job with slurm
   
   ```
   sbatch submit_tract_analysis.sh
   ```

## Modify the tracts

Upload the ROIs to the ```atlas_rois``` folder, and adapt the ```output_unravel.py``` and ```tracking_from_rois.py``` files.

## Skip the registration step

To gain time, if the regions are already registered but the tracking must be recomputed, the ```submit_tracking.sh``` file can be modified by commenting the registration step.

## Skip the tracking step

If the tracking is already done and you only want to run the UNRAVEL analysis, the ```submit_tract_analysis.sh``` file can be modified by replacing the ```submit_tracking.sh``` by the ```skip.sh``` step. This will bypass the tracking step.

>[!note]
>This readme is a work in progress, do not hesitate to leave suggestions on how to improve it.
