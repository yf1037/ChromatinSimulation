import os
import subprocess
import time
import numpy 
import numpy as np
from decimal import Decimal
from numpy import (atleast_1d, eye, mgrid, argmin, zeros, shape, squeeze,
                   asarray, sqrt, Inf, asfarray, isinf)

def writeCFR(ind,sim,fsim,maxiter,maxfun,xatol,fatol,rho,chi,psi,sigma,one2np1,params):
    LENGTH = 1920
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
    fh.write('params = {}\n'.format(params))
    fh.write('LENGTH = {}\n'.format(LENGTH))
    fh.write('stallRightARray = {}\n'.format([round(Decimal(i),3) for i in stallRightARray])) #use Decimal to limit to 3 decimal points
    fh.write('stallLeftArray = {}\n'.format([round(Decimal(i),3) for i in stallLeftArray]))
    fh.write('stallLeft list = {}\n'.format([round(Decimal(i),3) for i in stallL]))
    fh.write('stallRight list = {}\n'.format([round(Decimal(i),3) for i in stallR]))
    fh.write('sticky = {}\n'.format([round(Decimal(i),3) for i in sticky]))
    fh.write('domains = {}\n'.format([round(Decimal(i),3) for i in domains]))
    fh.write('sim = {}\n'.format([round(Decimal(i),3) for i in sim]))
    fh.write('fsm = {}\n'.format([float(i) for i in fsim]))
    fh.write('maxiter = {}\n'.format(maxiter))
    fh.write('maxfun = {}\n'.format(float(maxfun)))
    fh.write('xatol = {}\n'.format(xatol))
    fh.write('fatol = {}\n'.format(fatol))
    fh.write('rho = {}\n'.format(float(rho)))
    fh.write('chi = {}\n'.format(float(chi)))
    fh.write('psi = {}\n'.format(float(psi)))
    fh.write('sigma = {}\n'.format(float(sigma)))
    fh.write('one2np1 = {}\n'.format(one2np1))

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

# Create a directory to hold the parameter sets if one does not exist already
if not os.path.exists('parameter_output'):
    os.mkdir("parameter_output")

start = time.time()
print('Initializing..')

#setting parameters for minimization
x0=[0.8, 0.4, 0.0, 1.0, 0.8, 0.9, 0.0, 1.0, 0.0, 0.1, 0.0, 0.1, 0.0]
callback=None
maxiter=500
maxfev=None
disp=False
return_all=False
initial_simplex=None
xatol=0.03
fatol=0.03
adaptive=False
maxfun=None

"""
    Minimization of scalar function of one or more variables using the
    Nelder-Mead algorithm.
    Options
    -------
    disp : bool
        Set to True to print convergence messages.
    maxiter, maxfev : int
        Maximum allowed number of iterations and function evaluations.
        Will default to ``N*200``, where ``N`` is the number of
        variables, if neither `maxiter` or `maxfev` is set. If both
        `maxiter` and `maxfev` are set, minimization will stop at the
        first reached.
    initial_simplex : array_like of shape (N + 1, N)
        Initial simplex. If given, overrides `x0`.
        ``initial_simplex[j,:]`` should contain the coordinates of
        the j-th vertex of the ``N+1`` vertices in the simplex, where
        ``N`` is the dimension.
    xatol : float, optional
        Absolute error in xopt between iterations that is acceptable for
        convergence.
    fatol : number, optional
        Absolute error in func(xopt) between iterations that is acceptable for
        convergence.
    adaptive : bool, optional
        Adapt algorithm parameters to dimensionality of problem. Useful for
        high-dimensional minimization [1]_.
    References
    ----------
    .. [1] Gao, F. and Han, L.
       Implementing the Nelder-Mead simplex algorithm with adaptive
       parameters. 2012. Computational Optimization and Applications.
       51:1, pp. 259-277
"""
fcalls = [0]

if adaptive:
    dim = float(len(x0))
    rho = 1
    chi = 1 + 2/dim
    psi = 0.75 - 1/(2*dim)
    sigma = 1 - 1/dim
else:
    rho = 1
    chi = 2
    psi = 0.5
    sigma = 0.5

nonzdelt = 0.3 #0.05
#zdelt = 0.1 #0.00025

x0 = asfarray(x0).flatten()

if initial_simplex is None:
    N = len(x0)

    sim = numpy.zeros((N + 1, N), dtype=x0.dtype)
    sim[0] = x0
    for k in range(N):
        y = numpy.array(x0, copy=True)
        if y[k] >0.5:
            y[k] = y[k]-nonzdelt 
        else:
            y[k] = y[k]+nonzdelt
        sim[k + 1] = y
else:
    sim = np.asfarray(initial_simplex).copy()
    if sim.ndim != 2 or sim.shape[0] != sim.shape[1] + 1:
        raise ValueError("`initial_simplex` should be an array of shape (N+1,N)")
    if len(x0) != sim.shape[1]:
        raise ValueError("Size of `initial_simplex` is not consistent with `x0`")
    N = sim.shape[1]

#if retall:
#    allvecs = [sim[0]]

# If neither are set, then set both to default
if maxiter is None and maxfun is None:
    maxiter = N * 200
    maxfun = N * 200
elif maxiter is None:
    maxiter = N * 200
elif maxfun is None:
    maxfun = N * 200

one2np1 = list(range(1, N + 1))
fsim = numpy.zeros((N + 1,), float)

#for testing
#print('initialized {} parameters'.format(sim.shape))
#maxiter=2
#N=1

#output all the settings and parmaters for subprocess
for i in range(N+1):
    sim[i] = checkparams(sim[i])
    writeCFR(i,sim[i],fsim,maxiter,maxfun,xatol,fatol,rho,chi,psi,sigma,one2np1,N)
print("initialized in {}".format(time.time() - start))

subprocess.run(["bash", "optimization_iteration.sh", "{}".format(N)])


