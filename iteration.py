import os
import numpy as np
import numpy
import re
import subprocess
import time
from decimal import Decimal

def writeCFR(ind,sim,LENGTH):
    stallL = sim[0:4]
    stallR = sim[4:8]
    stallLeftArray = np.zeros(LENGTH, dtype = np.double)
    stallRightARray = np.zeros(LENGTH, dtype = np.double)
    stallLeftArray[456]=stallL[0]
    stallLeftArray[808]=stallL[1]
    stallLeftArray[1178]=stallL[2]
    stallLeftArray[1592]=stallL[3]
    stallRightARray[456]=stallR[0]
    stallRightARray[808]=stallR[1]
    stallRightARray[1178]=stallR[2]
    stallRightARray[1592]=stallR[3]

    sticky = sim[8:13]
    domainInd=[456,808,1178,1592]
    domains=np.zeros(LENGTH)
    domains[0:domainInd[0]]=sticky[0]
    domains[domainInd[0]:domainInd[1]]=sticky[1]
    domains[domainInd[1]:domainInd[2]]=sticky[2]
    domains[domainInd[2]:domainInd[3]]=sticky[3]
    domains[domainInd[3]:LENGTH]=sticky[4]

    folder = "trajectory_{}".format(ind)

    fh = open('parameter_output/{}.txt'.format(ind), 'w')           # Open a new file in write mode
    fh.write('ind = {}\n'.format(ind))
    fh.write('LENGTH = {}\n'.format(LENGTH))
    fh.write('stallRightARray = {}\n'.format([round(Decimal(i),3) for i in stallRightARray])) # use Decimal to limit to 3 decimal points
    fh.write('stallLeftArray = {}\n'.format([round(Decimal(i),3) for i in stallLeftArray]))
    fh.write('stallLeft list = {}\n'.format([round(Decimal(i),3) for i in stallL]))
    fh.write('stallRight list = {}\n'.format([round(Decimal(i),3) for i in stallR]))
    fh.write('sticky = {}\n'.format([round(Decimal(i),3) for i in sticky]))
    fh.write('domains = {}\n'.format([round(Decimal(i),3) for i in domains]))
    fh.write('sim = {}\n'.format([round(Decimal(i),3) for i in sim]))
    
def checkparams(x):
    for i in range(8):
        if x[i] > 1:
            x[i] = 1
        if x[i] < 0:
            x[i] = 0
    for i in range(4):
        if x[i+8] > 2:
            x[i+8] = 2
        if x[i+8] < 0:
            x[i+8] = 0
    return x

start = time.time()
#read in results from first simulation as fsim, initial params as sim
#writeCFR(indx,sim,fsim,maxiter,maxfun,xatol,fatol,rho,chi,psi,sigma,one2np1)
with open('parameter_output/0.txt', 'r') as CFR:
    CFR.seek(0, 0)
    for params in CFR:
        if re.findall('params',params):
            N = int(re.findall('\d.*?\W', params)[0])
        if re.findall('LENGTH', params):
            LENGTH = int(re.findall('\d.*?\W', params)[0])
        if re.findall('maxiter', params):
            maxiter = int(re.findall('\d.*?\W', params)[0])
        if re.findall('maxfun', params):
            maxfun = float(re.findall('\d+\.\d*', params)[0])
        if re.findall('xatol', params):
            xatol = float(re.findall('\d+\.\d*', params)[0])
        if re.findall('fatol', params):
            fatol = float(re.findall('\d+\.\d*', params)[0])
        if re.findall('rho', params):
            rho = float(re.findall('\d+\.\d*', params)[0])
        if re.findall('chi', params):
            chi = float(re.findall('\d+\.\d*', params)[0])
        if re.findall('psi', params):
            psi = float(re.findall('\d+\.\d*', params)[0])
        if re.findall('sigma', params):
            sigma = float(re.findall('\d+\.\d*', params)[0])
        if re.findall('one2np1', params):
            one2np1 = np.array([int(i) for i in re.findall('\d', params)])

sim=np.zeros((N+1,N))
fsim=np.zeros(N+1)
Ctype='T cells' #optimizing parameters for T cells
fcalls=[0]

for indx in range(N+1):
    with open('parameter_output/{}.txt'.format(indx), 'r') as CFR:
        CFR.seek(0, 0)
        for params in CFR:
            if re.findall('sim', params):
                #for testing
                #sim[indx] = float(re.findall('\d+\.\d*', params)[0])
                sim[indx] = np.array([float(i) for i in re.findall('\d+\.\d*', params)],dtype=np.double)

for indx in range(N+1):
    with open('correlation_{}.txt'.format(indx), 'r') as CFR:
        CFR.seek(0, 0)
        for params in CFR:
            if re.findall(Ctype, params): 
                fsim[indx] = 1-float(re.findall('\d+\.\d*', params)[0])

#initial processing
ind = numpy.argsort(fsim)
fsim = numpy.take(fsim, ind, 0)
# sort so sim[0,:] has the lowest function value
sim = numpy.take(sim, ind, 0)

iterations = 1
indx = N+1

while (fcalls[0] < maxfun and iterations < maxiter):
    if (numpy.max(numpy.ravel(numpy.abs(sim[1:] - sim[0]))) <= xatol and
            numpy.max(numpy.abs(fsim[0] - fsim[1:])) <= fatol):
        break

    #This is the directory created in the driver.py scripts
    PARAMS_DIR_ROOT = 'parameter_output'
    # This is the job id for this job - used by the next step to identify the parameters
    #PARAM_JOB_ID = os.environ["SLURM_JOB_ID"]
    # Output for the current iterations params
    #CURRENT_PARAM_DIR = os.path.join(PARAMS_DIR_ROOT, PARAM_JOB_ID)

    # This folder should not ever exist, but we still check
    #if not os.path.exists(CURRENT_PARAM_DIR):
    #    os.mkdir(CURRENT_PARAM_DIR)

    print("generating params iteration {}".format(iterations))
#    print("Run time {}s".format(time.time()-start))

    # Write the parameters to a file for each of the jobs that are going to run
    xbar = numpy.add.reduce(sim[:-1], 0) / N
    xr = (1 + rho) * xbar - rho * sim[-1]
    #fxr = func(xr)
    xr = checkparams(xr)
    writeCFR(indx,xr,LENGTH)
    subprocess.run(["sh", "optimize.sh", "{}".format(indx)])
    #read result of the simulation
    if not os.path.exists('correlation_{}.txt'.format(indx)):
        print('Failed to wait for simulation to finish')
        break
#    time0=time.time()
#    time.sleep(8*60)#wait for 8min
#    while not os.path.exists('correlation_{}.txt'.format(indx)):
#        time.sleep(10)#try again in 10s
#        if (time.time()-time0 > 2*60*60):
#            print('fxr: Simulation failed to finish in 2h')
#            break
#    print('Simulation finished in {}s'.format(time.time()-time0))
    with open('correlation_{}.txt'.format(indx), 'r') as CFR:
        CFR.seek(0, 0)
        for params in CFR:
            if re.findall(Ctype, params): 
                fxr = 1-float(re.findall('\d+\.\d*', params)[0]) 
    print('Tested fxr')
    indx +=1
    doshrink = 0

    if fxr < fsim[0]:
        xe = (1 + rho * chi) * xbar - rho * chi * sim[-1]
        #fxe = func(xe)
        xe = checkparams(xe)
        writeCFR(indx,xe,LENGTH)
        subprocess.run(["sh", "optimize.sh", "{}".format(indx)])
#        time0=time.time()
#        time.sleep(8*60)#wait 8min
#        while not os.path.exists('correlation_{}.txt'.format(indx)):
#            time.sleep(10)#try again in 10s
#            if (time.time()-time0 > 2*60*60):
#                print('fxe: Simulation failed to finish in 2h')
#                break
        #read result of the simulation
        with open('correlation_{}.txt'.format(indx), 'r') as CFR:
            CFR.seek(0, 0)
            for params in CFR:
                if re.findall(Ctype, params):
                    fxe = 1-float(re.findall('\d+\.\d*', params)[0])
        print('Tested fxe')
        ind += 1

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
                #fxc = func(xc)
                xc = checkparams(xc)
                writeCFR(indx,xc,LENGTH)
                subprocess.run(["sh", "optimize.sh", "{}".format(indx)])
#                time0=time.time()
#                time.sleep(8*60)#wait 8min
#                while not os.path.exists('correlation_{}.txt'.format(indx)):
#                    time.sleep(10)#try again in 10s
#                    if (time.time()-time0 > 2*60*60):
#                        print('fxc: Simulation failed to finish in 2h')
#                        break
                #read result of the simulation
                with open('correlation_{}.txt'.format(indx), 'r') as CFR:
                    CFR.seek(0, 0)
                    for params in CFR:
                        if re.findall(Ctype, params):
                            fxc = 1-float(re.findall('\d+\.\d*', params)[0])
                print('Tested fxc')
                indx += 1

                if fxc <= fxr:
                    sim[-1] = xc
                    fsim[-1] = fxc
                else:
                    doshrink = 1
            else:
                # Perform an inside contraction
                xcc = (1 - psi) * xbar + psi * sim[-1]
                #fxcc = func(xcc)
                xcc = checkparams(xcc)
                writeCFR(indx,xcc,LENGTH)
                subprocess.run(["sh", "optimize.sh", "{}".format(indx)])
#                time0=time.time()
#                time.sleep(8*60)#wait 8min
#                while not os.path.exists('correlation_{}.txt'.format(indx)):
#                    time.sleep(10)#try again in 10s
#                    if (time.time()-time0 > 2*60*60):
#                        print('fxcc: Simulation failed to finish in 2h')
#                        break
                #read result of the simulation
                with open('correlation_{}.txt'.format(indx), 'r') as CFR:
                    CFR.seek(0, 0)
                    for params in CFR:
                        if re.findall(Ctype, params):
                            fxcc = 1-float(re.findall('\d+\.\d*', params)[0])
                print('Tested fxcc')
                indx += 1

                if fxcc < fsim[-1]:
                    sim[-1] = xcc
                    fsim[-1] = fxcc
                else:
                    doshrink = 1

            if doshrink:
                for j in one2np1:
                    sim[j] = sim[0] + sigma * (sim[j] - sim[0])
                    #fsim[j] = func(sim[j])
                    sim[j] = checkparams(sim[j])
                    writeCFR(indx,sim[j],LENGTH)
                    subprocess.run(["sh", "optimize.sh", "{}".format(indx)])
#                    time0=time.time()
#                    time.sleep(8*60)#wait 8min
#                    while not os.path.exists('correlation_{}.txt'.format(indx)):
#                        time.sleep(10)#try again in 10s
#                        if (time.time()-time0 > 2*60*60):
#                            print('doshrink: Simulation failed to finish in 2h')
#                            break
                    #read result of the simulation
                    with open('correlation_{}.txt'.format(indx), 'r') as CFR:
                        CFR.seek(0, 0)
                        for params in CFR:
                            if re.findall(Ctype, params):
                                fsim[j] = 1-float(re.findall('\d+\.\d*', params)[0])
                    indx += 1
                print('Tested doshrink')

    ind = numpy.argsort(fsim)
    sim = numpy.take(sim, ind, 0)
    fsim = numpy.take(fsim, ind, 0)
#    if callback is not None:
#        callback(sim[0])
    iterations += 1
#    if retall:
#        allvecs.append(sim[0])


#output final results
x = sim[0]
fval = numpy.min(fsim)

#writeCFR('best',x,fval,maxiter,maxfun,xatol,fatol,rho,chi,psi,sigma,one2np1)

stallL = x[0:4]
stallR = x[4:8]
stallLeftArray = np.zeros(LENGTH, dtype = np.double)
stallRightARray = np.zeros(LENGTH, dtype = np.double)
stallLeftArray[456]=stallL[0]
stallLeftArray[808]=stallL[1]
stallLeftArray[1178]=stallL[2]
stallLeftArray[1592]=stallL[3]
stallRightARray[456]=stallR[0]
stallRightARray[808]=stallR[1]
stallRightARray[1178]=stallR[2]
stallRightARray[1592]=stallR[3]

sticky = x[8:13]
domainInd=[456,808,1178,1592]
domains=np.zeros(LENGTH)
domains[0:domainInd[0]]=sticky[0]
domains[domainInd[0]:domainInd[1]]=sticky[1]
domains[domainInd[1]:domainInd[2]]=sticky[2]
domains[domainInd[2]:domainInd[3]]=sticky[3]
domains[domainInd[3]:LENGTH]=sticky[4]


fh = open('parameter_output/best.txt', 'w')           # Open a new file in write mode
fh.write('stallRightARray = {}\n'.format([float(i) for i in stallRightARray]))
fh.write('stallLeftArray = {}\n'.format([float(i) for i in stallLeftArray]))
fh.write('stallLeft list = {}\n'.format([float(i) for i in stallL]))
fh.write('stallRight list = {}\n'.format([float(i) for i in stallR]))
fh.write('sticky = {}\n'.format([float(i) for i in sticky]))
fh.write('domains = {}\n'.format([float(i) for i in domains]))
fh.write('sim = {}\n'.format([float(i) for i in x]))
fh.write('fsim = {}\n'.format(fval))
fh.write('maxiter = {}\n'.format(maxiter))
fh.write('maxfun = {}\n'.format(maxfun))
fh.write('xatol = {}\n'.format(xatol))
fh.write('fatol = {}\n'.format(fatol))
fh.write('rho = {}\n'.format(rho))
fh.write('chi = {}\n'.format(chi))
fh.write('psi = {}\n'.format(psi))
fh.write('sigma = {}\n'.format(sigma))
fh.write('iteration = {}\n'.format(iterations))
fh.write('simulations performed: {}'.format(indx))
