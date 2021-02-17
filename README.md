# ChromatinSimulation  
These scripts use openmmlib (https://github.com/mirnylab/openmm-polymer-legacy) and simplex algorithm to simulate HiC map of a given section on a chromatin (e.g. hg19 chr8:129M-131M) on High Performance Clusters. It is setup to exam changes in TAD boudary strength and compartment strength between different cell lines/ conditions.  
  
How it works:  
1. Use sbatch to submit driver.sh  
2. driver.sh submits a job on a cpu node to run driver.py  
3. driver.py sets parameters for simplex including initial guesses, then subprocesses optimization_iteration.sh  
4. optimization_iteration.sh submit initial simulation jobs then submit iteration.sh  
5. iteration.sh submits iteration.py on a cpu node  
6. In each iteration, iteration.py generates next guess parameter set based on current correlations.  
  Ending criteria: Reach max iterations/ correlation stop improving and parameters stop changing  
  
Each simulation job includes:  
Submit simulation as 10 array jobs, each on a GPU node  
Calculate contact frequency maps from simulations (CPU node)  
Calculate Spearman correlation of HiC map and simulated contact frequency (CPU node)  
  
Simulation principle:  
Use loop extrusion to model TADs (dynamic simulation)  
Use sticky to model compartment interactions (steady-state simulation)  
Chromatin is modeled as self avoiding chain  
