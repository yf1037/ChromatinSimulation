import joblib
#import pickle
import numpy as np
import os
import time
import pandas as pd
import scipy
import scipy.stats
import sys

#-------read configurations-----------
import argparse
import re

# initialize time
start_time = time.time();

# Create a parser:
parser = argparse.ArgumentParser(description='Simulation')

# Add positional arguments
parser.add_argument('file',                  # positional argument name;
                    nargs=1,           # num of argument values expected for this arg
                    type=str,          # for type checking at interface level
                    help='file name')

args = parser.parse_args()

with open('parameter_output/{}.txt'.format(args.file[0]), 'r') as CFR:
    CFR.seek(0, 0)
    for params in CFR:
        if re.findall('LENGTH',params):
            LENGTH = int(re.findall('\d', params)[0])

#plot Hi-C map

thr = 600
thr2 = thr*thr
#start_compute_time = time.time();

for i in range(10):#range(1) for testing, default range(10)
    for k in range(100):

        block=joblib.load('trajectory_{}_{}/block{}.dat'.format(args.file[0],i,k+1))
        data=block['data']

        N = len(data) # just shorter, could also use len(s)
        if k==0:
            HiC=np.zeros([N,N]);
            HiC0=np.zeros([N,N]);

        [x,y,z]=data.transpose()*50

        dx1 = np.repeat(x.reshape(N,1), N, axis=1)
        dx2 = np.repeat(x.reshape(1,N), N, axis=0)
        dx = dx1-dx2

        dy1 = np.repeat(y.reshape(N,1), N, axis=1)
        dy2 = np.repeat(y.reshape(1,N), N, axis=0)
        dy = dy1-dy2

        dz1 = np.repeat(z.reshape(N,1), N, axis=1)
        dz2 = np.repeat(z.reshape(1,N), N, axis=0)
        dz = dz1-dz2

        dis = dx*dx + dy*dy +dz*dz ;

        HiC0 = dis < thr2
        HiC=np.add(HiC,HiC0)
#print("--- performed computation in %s seconds ---" % (time.time() - start_compute_time))
HiC=np.divide(HiC,200)
#start_save_time = time.time();
np.savetxt('ContactFrequency_{}.txt'.format(args.file[0]),HiC)
print('Saved file as ContactFrequency_{}.txt'.format(args.file[0]))
#print("--- saved in %s seconds ---" % (time.time() - start_time))

HiC2=np.zeros([124,N])
j=0
while j < N :
    i = 0
    while (i < 16) & (j < N):
        HiC2[j//16,:]=np.add(HiC2[j//16,:],HiC[j,:])
        i+=1
        j+=1

#print("Shrinked one dimension in %s seconds" % (time.time() - start_time))

HiC3=np.zeros([124,124])
j=0
while j < N :
    i = 0
    while (i < 16) & (j<N) :
        HiC3[:,j//16]=np.add(HiC3[:,j//16],HiC2[:,j])
        i+=1
        j+=1

HiC=np.divide(HiC3,256)
np.savetxt('ContactFrequency_{}_40k.txt'.format(args.file[0]), HiC)
print('Saved file as ContactFrequency_{}_40k.txt'.format(args.file[0]))
#print("--- saved in %s seconds ---" % (time.time() - start_time))

data = pd.read_csv("../../../../HiCmatrix/CUTLL1_DMSO_repeat2_matrix.chr8.tsv",sep='\t')
HiC0=data.iloc[3168:3292,3168:3292]
del data
HiC0.astype('float64')

data1 = pd.read_csv("../../../../HiCmatrix/T_cell-Arima_matrix.chr8.tsv",sep='\t')
HiC1=data1.iloc[3168:3292,3168:3292]
del data1
HiC1.astype('float64')

with open('correlation_{}.txt'.format(args.file[0]), 'w') as fh:
    fh.write('\tRho\tPval\n')

    rho, pval = scipy.stats.spearmanr(HiC, HiC0, axis=None)
    fh.write('CUTLL\t{}\t{}\n'.format(rho, pval))

    rho, pval = scipy.stats.spearmanr(HiC, HiC1, axis=None)
    fh.write('T cells\t{}\t{}\n'.format(rho, pval))

print("Saved correlation_{}.txt".format(args.file[0]))
#print("--- exectued code in %s seconds ---" % (time.time() - start_time))

