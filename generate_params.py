import os

# This is the directory created in the driver.py scripts
PARAMS_DIR_ROOT = 'parameter_output'
# This is the job id for this job - used by the next step to identify the parameters
PARAM_JOB_ID = os.environ["SLURM_JOB_ID"]
# Output for the current iterations params
CURRENT_PARAM_DIR = os.path.join(PARAMS_DIR_ROOT, PARAM_JOB_ID)

# This folder should not ever exist, but we still check
if not os.path.exists(CURRENT_PARAM_DIR):
    os.mkdir(CURRENT_PARAM_DIR)

print("generating params: {}".format(PARAM_JOB_ID))

# Write the parameters to a file for each of the jobs that are going to run
# I hard coded the job to consist of 4 subjobs that will run together
# I just write a file for each of the arrays jobs, this could also be a folder
#    What is important is the naming, I set the array to be from 1-4 so this must be consistant
xbar = numpy.add.reduce(sim[:-1], 0) / N
xr = (1 + rho) * xbar - rho * sim[-1]
fxr = func(xr)
doshrink = 0

if fxr < fsim[0]:
    xe = (1 + rho * chi) * xbar - rho * chi * sim[-1]
    fxe = func(xe)

    if fxe < fxr:
        sim[-1] = xe
        fsim[-1] = fxe
    else:
        sim[-1] = xr
        fsim[-1] = fxr
else:  # fsim[0] <= fxr
    if fxr < fsim[-2]:
        sim[-1] = xr
        fsim[-1] = fxr
    else:  # fxr >= fsim[-2]
        # Perform contraction
        if fxr < fsim[-1]:
            xc = (1 + psi * rho) * xbar - psi * rho * sim[-1]
            fxc = func(xc)

            if fxc <= fxr:
                sim[-1] = xc
                fsim[-1] = fxc
            else:
                doshrink = 1
        else:
            # Perform an inside contraction
            xcc = (1 - psi) * xbar + psi * sim[-1]
            fxcc = func(xcc)

            if fxcc < fsim[-1]:
                sim[-1] = xcc
                fsim[-1] = fxcc
            else:
                doshrink = 1

        if doshrink:
            for j in one2np1:
                sim[j] = sim[0] + sigma * (sim[j] - sim[0])
                fsim[j] = func(sim[j])
ind = numpy.argsort(fsim)
sim = numpy.take(sim, ind, 0)
fsim = numpy.take(fsim, ind, 0)
if callback is not None:
    callback(sim[0])
iterations += 1
if retall:
    allvecs.append(sim[0])

PARAM_SET = ["param1", "param2", "param3", "param4"]
for i in range(4):
    # Create a new file in the correct parameter directory for the array job
    with open(os.path.join(CURRENT_PARAM_DIR, '{}.txt'.format(i + 1)), 'w') as fh:
        fh.write("{}, {}\n".format(PARAM_SET[i], PARAM_JOB_ID))
