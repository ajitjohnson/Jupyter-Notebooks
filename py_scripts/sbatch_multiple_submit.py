#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 19:53:17 2021
@author: aj
"""

# %% Lib
import sys
import os
import pandas as pd
import subprocess

# %% WD
#cwd = os.getcwd()

# %% Import parameters
csv_file =  sys.argv[1] # prints python_script.py


# %% Data
# read the CSV file
#config = pd.read_csv(str(cwd) + '/' + str(csv_file))
config = pd.read_csv(csv_file)

try:
    # convert the description column to catergory
    config['description'] = config['description'].astype('category')
    
    for i in config['description'].cat.categories:
        read_1 = config[config['description'] == i]['samplename'].iloc[0]
        read_2 = config[config['description'] == i]['samplename'].iloc[1]
        name = i
        # sumbit command
        submit_command = ("sbatch " + "--job-name=" + str(name) +
                " --export=name=" + str(name) + ',' +          
                "read_1=" + str(read_1) + ',' +
                "read_2=" + str(read_2) + " sbatch_template.sh")
        #print(submit_command)
        # submit job
        exit_status = subprocess.call(submit_command, shell=True)
        if exit_status is 1:  # Check to make sure the job submitted
            print("Job" + str(name) + "failed to submit")
    
except Exception:
    pass   
    
