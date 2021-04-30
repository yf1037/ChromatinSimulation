import sys
import pickle
import os
import time
import numpy as np
from openmmlib import polymerutils
from openmmlib.polymerutils import scanBlocks
from openmmlib.openmmlib import Simulation
from openmmlib.polymerutils import grow_rw

sys.path.append('/gpfs/home/yf1037/simulation/codes/examples/loopExtrusion/')

import pyximport; pyximport.install()
from smcTranslocator import smcTranslocatorDirectional

import openmmlib.polymerutils as polymerutils
import scipy.io
import CompInt
from CompInt import InvertedNucleiSimulation #for compartment interactions

#-------Command line interface-------
import argparse
import re

# Create a parser:
parser = argparse.ArgumentParser(description='Simulation')

# Add positional arguments
parser.add_argument('CFR',                  # positional argument name;
                    nargs=1,           # num of argument values expected for this arg
                    type=str,          # for type checking at interface level
                    help='Configuration file')
parser.add_argument('ID',
                    nargs=1,
                    type=int,
                    help='Array ID')

args = parser.parse_args()
ID = int(args.ID[0])
with open('parameter_output/{}.txt'.format(args.CFR[0]), 'r') as CFR:
    CFR.seek(0, 0)
    for params in CFR:
        if re.findall('folder', params):
            folderID = re.findall('\w+', params)[-1]
        if re.findall('stallLeftArray', params):
            stallLeftArray = np.array([float(i) for i in re.findall('\d+\.\d*', params)],dtype=np.double)
        if re.findall('stallRightARray', params):
            stallRightARray = np.array([float(i) for i in re.findall('\d+\.\d*', params)],dtype=np.double)
        if re.findall('domains', params):
            domains = np.array([float(i) for i in re.findall('\d+\.\d*', params)],dtype=np.int)
        if re.findall('LENGTH', params):
            N = int(re.findall('\d.*?\W', params)[0])


# -------defining parameters----------
#  -- basic loop extrusion parameters
SEPARATION = 200
LIFETIME = 300
#N = 1920   # number of monomers
smcStepsPerBlock = 1  # now doing 1 SMC step per block
steps = 250   # steps per block (now extrusion advances by one step per block)
stiff = 2
dens = 0.2
box = (N / dens) ** 0.33  # density = 0.1.
data = polymerutils.grow_rw(N, int(box) - 2)  # creates a compact conformation

#create initial configurations for polymer - random walk in this case
#data = polymerutils.create_random_walk(1,N)

block = 0  # starting block

#folder 
folder = "trajectory_{}_{}".format(args.CFR[0],ID-1)

#  -- basic compartment interaction parameters
inv_chromFrac=6
superLength=N/inv_chromFrac

#define compartment identies (there are three compartments here, see inverted nucleus paper for details)
#domains = np.zeros(N)

# new parameters because some things changed 
saveEveryBlocks = 300   # save every 10 blocks (saving every block is now too much almost)
skipSavedBlocksBeginning = 4  # how many blocks (saved) to skip after you restart LEF positions
totalSavedBlocks = 100  # how many blocks to save (number of blocks done is totalSavedBlocks * saveEveryBlocks)
restartMilkerEveryBlocks = 600

# parameters for smc bonds 
smcBondWiggleDist = 0.2
smcBondDist = 0.5

# assertions for easy managing code below 

assert restartMilkerEveryBlocks % saveEveryBlocks == 0 
assert (skipSavedBlocksBeginning * saveEveryBlocks) % restartMilkerEveryBlocks == 0 
assert (totalSavedBlocks * saveEveryBlocks) % restartMilkerEveryBlocks == 0 

savesPerMilker = restartMilkerEveryBlocks // saveEveryBlocks
milkerInitsSkip = saveEveryBlocks * skipSavedBlocksBeginning  // restartMilkerEveryBlocks
milkerInitsTotal  = (totalSavedBlocks + skipSavedBlocksBeginning) * saveEveryBlocks // restartMilkerEveryBlocks
print("Milker will be initialized {0} times, first {1} will be skipped".format(milkerInitsTotal, milkerInitsSkip))

class smcTranslocatorMilker(object):

    def __init__(self, smcTransObject):
        """
        :param smcTransObject: smc translocator object to work with
        """
        self.smcObject = smcTransObject
        self.allBonds = []

    def setParams(self, activeParamDict, inactiveParamDict):
        """
        A method to set parameters for bonds.
        It is a separate method because you may want to have a Simulation object already existing

        :param activeParamDict: a dict (argument:value) of addBond arguments for active bonds
        :param inactiveParamDict:  a dict (argument:value) of addBond arguments for inactive bonds

        """
        self.activeParamDict = activeParamDict
        self.inactiveParamDict = inactiveParamDict


    def setup(self, bondForce,  blocks = 100, smcStepsPerBlock = 1):
        """
        A method that milks smcTranslocator object
        and creates a set of unique bonds, etc.

        :param bondForce: a bondforce object (new after simulation restart!)
        :param blocks: number of blocks to precalculate
        :param smcStepsPerBlock: number of smcTranslocator steps per block
        :return:
        """


        if len(self.allBonds) != 0:
            raise ValueError("Not all bonds were used; {0} sets left".format(len(self.allBonds)))

        self.bondForce = bondForce

        #precalculating all bonds
        allBonds = []
        for dummy in range(blocks):
            self.smcObject.steps(smcStepsPerBlock)
            left, right = self.smcObject.getSMCs()
            bonds = [(int(i), int(j)) for i,j in zip(left, right)]
            allBonds.append(bonds)

        self.allBonds = allBonds
        self.uniqueBonds = list(set(sum(allBonds, [])))

        #adding forces and getting bond indices
        self.bondInds = []
        self.curBonds = allBonds.pop(0)

        for bond in self.uniqueBonds:
            paramset = self.activeParamDict if (bond in self.curBonds) else self.inactiveParamDict
            ind = bondForce.addBond(bond[0], bond[1], **paramset)
            self.bondInds.append(ind)
        self.bondToInd = {i:j for i,j in zip(self.uniqueBonds, self.bondInds)}
        return self.curBonds,[]


    def step(self, context, verbose=False):
        """
        Update the bonds to the next step.
        It sets bonds for you automatically!
        :param context:  context
        :return: (current bonds, previous step bonds); just for reference
        """
        if len(self.allBonds) == 0:
            raise ValueError("No bonds left to run; you should restart simulation and run setup  again")

        pastBonds = self.curBonds
        self.curBonds = self.allBonds.pop(0)  # getting current bonds
        bondsRemove = [i for i in pastBonds if i not in self.curBonds]
        bondsAdd = [i for i in self.curBonds if i not in pastBonds]
        bondsStay = [i for i in pastBonds if i in self.curBonds]
        if verbose:
            print("{0} bonds stay, {1} new bonds, {2} bonds removed".format(len(bondsStay),
                                                                            len(bondsAdd), len(bondsRemove)))
        bondsToChange = bondsAdd + bondsRemove
        bondsIsAdd = [True] * len(bondsAdd) + [False] * len(bondsRemove)
        for bond, isAdd in zip(bondsToChange, bondsIsAdd):
            ind = self.bondToInd[bond]
            paramset = self.activeParamDict if isAdd else self.inactiveParamDict
            self.bondForce.setBondParameters(ind, bond[0], bond[1], **paramset)  # actually updating bonds
        self.bondForce.updateParametersInContext(context)  # now run this to update things in the context
        return self.curBonds, pastBonds

def initModel():
    # this jsut inits the simulation model. Put your previous init code here 
    birthArray = np.zeros(N, dtype=np.double) + 0.1
    deathArray = np.zeros(N, dtype=np.double) + 1. / LIFETIME
    stallDeathArray = np.zeros(N, dtype=np.double) + 1 / LIFETIME
    pauseArray = np.zeros(N, dtype=np.double)

    #stallList = [448, 784, 1536]#disrupted boundary 1130
    #stallLeftArray = np.zeros(N, dtype = np.double)
    #stallRightARray = np.zeros(N, dtype = np.double)
    #for i in stallList:
    #    stallLeftArray[i] = 0.8
    #    stallRightARray[i] = 0.8

    smcNum = N // SEPARATION
    SMCTran = smcTranslocatorDirectional(birthArray, deathArray, stallLeftArray, stallRightARray, pauseArray,
                                         stallDeathArray, smcNum)
    return SMCTran


SMCTran = initModel()  # defining actual smc translocator object 


# now polymer simulation code starts

# ------------feed smcTran to the milker---
SMCTran.steps(1000000)  # first steps to "equilibrate" SMC dynamics. If desired of course. Default 1000000
milker = smcTranslocatorMilker(SMCTran)   # now feed this thing to milker (do it once!)
#--------- end new code ------------

for milkerCount in range(milkerInitsTotal):
    doSave = milkerCount >= milkerInitsSkip

    # create object
    a = InvertedNucleiSimulation(timestep=80, thermostat=0.01)
    # specifiy simulation meta-options
    a.setup(platform="cuda", integrator="variableLangevin", errorTol=0.01, GPU=0, PBC=True, PBCbox=[box, box, box], precision="mixed")
    a.saveFolder(folder)

    ## simulation parameters are defined below
    #a = Simulation(timestep=80, thermostat=0.01)
    #a.setup(platform="CUDA", PBC=True, PBCbox=[box, box, box], GPU=0, precision="mixed")  # set up GPU here

    # assign stickinesses based on compartment assignment
#    chrm_stickyValues = [sticky if i else 0 for i in domains]
        #[1.5 if (s % chromLength < superLength) else .45 if (
         #       signal[s - (s // chromLength + 1) * superLength] == 0) else 0.04 for s in np.arange(chromLength * 2)]
    stickyValues = domains

    # set domains (important for some forces)
    #a.useDomains(domains)
    ## set breaks between chromosomes so it's not just one big one
    #a.setChains([(i * Nm // chromNum, (i + 1) * Nm // chromNum, 0) for i in np.arange(chromNum)])

    # add your actual force terms
    #a.addSphericalConfinement(k=5., density=0.35)
    a.addSuperDomainSSWForce(stickiness=stickyValues, repulsionEnergy=3.0, attractionRadius=1.8, N=N) #length scale??
    #repulsionE default 5.0
    #a.addHarmonicPolymerBonds(wiggleDist=0.05)
    #a.addStiffness(k=0.5)

    a.load(data)
    a.addHarmonicPolymerBonds(wiggleDist=0.1)
    if stiff > 0:
        a.addGrosbergStiffness(stiff)
    #a.addPolynomialRepulsiveForce(trunc=1.5, radiusMult=1.05)
    a.step = block
    #a.save(filename=folder + "/block%d.dat" % (j * sliceBlock))

    # ------------ initializing milker; adding bonds ---------
    # copied from addBond
    kbond = a.kbondScalingFactor / (smcBondWiggleDist ** 2)
    bondDist = smcBondDist * a.length_scale

    activeParams = {"length":bondDist,"k":kbond}
    inactiveParams = {"length":bondDist, "k":0}
    milker.setParams(activeParams, inactiveParams)
     
    # this step actually puts all bonds in and sets first bonds to be what they should be
    milker.setup(bondForce=a.forceDict["HarmonicBondForce"],
                blocks=restartMilkerEveryBlocks,   # default value; milk for 100 blocks
                 smcStepsPerBlock=smcStepsPerBlock)  # now only one step of SMC per step
    print("Restarting milker")

    a.doBlock(steps=steps, increment=False)  # do block for the first time with first set of bonds in

    for i in range(restartMilkerEveryBlocks - 1):
        curBonds, pastBonds = milker.step(a.context)  # this updates bonds. You can do something with bonds here
        #curBonds.append([(798,1427)])#E-P interaction
        if i % saveEveryBlocks == (saveEveryBlocks - 2):
            a.doBlock(steps=steps, increment = doSave)    
            if doSave: 
                a.save()
                pickle.dump(curBonds, open(os.path.join(a.folder, "SMC{0}.dat".format(a.step)),'wb'))
        else:
            a.integrator.step(steps)  # do steps without getting the positions from the GPU (faster)

    data = a.getData()  # save data and step, and delete the simulation
    block = a.step
    del a

    time.sleep(0.2)  # wait 200ms for sanity (to let garbage collector do its magic)

