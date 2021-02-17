import os
import sys
import time

# This is the id of the parameters job which generated parameters for this job
PARAM_JOB_ID = sys.argv[1]
# This was set in driver.py
PARAMS_DIR_ROOT = 'parameter_output'
# This the the paramter directory for this set of jobs
CURRENT_PARAM_DIR = os.path.join(PARAMS_DIR_ROOT, PARAM_JOB_ID)
# This is the ID of the current job in the array (1-4 in my example code)
ARRAY_TASK_ID = os.environ["SLURM_ARRAY_TASK_ID"]

print("Optimization ARRAY Job ID: {}".format(os.environ["SLURM_ARRAY_JOB_ID"]))
print("\t Param Job ID {}".format(PARAM_JOB_ID))
print("\t Array ID: {}".format(ARRAY_TASK_ID))

# Read the params from the file written during the previous job
CURRENT_PARAM_FILE = os.path.join(CURRENT_PARAM_DIR, ARRAY_TASK_ID + '.txt')
with open(CURRENT_PARAM_FILE) as fh:
    print("\t PARAM FILE CONTENT: {}".format(fh.read()))

time.sleep(10)
