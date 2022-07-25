# Random-Walk-with-Cardiac-Diffusion-Tensor-Imaging
Random walk algorithms are extensively used in diffusion problems, especially when analytical solutions are unavailable from 2 dimensions onwards. Throughout history, Monte carlo simulations have been the core of countless breakthroughs in this field, from nuclear physics to thermodynamics.
<br />
<br />
This code provides an easy-to-use Monte Carlo simulator to simulate an imaging voxel. Given a sufficient number of walkers, one can study the effects of extra-cellular volume fraction and the distribution of extra-cellular space<sup>1</sup>. It imports a myocyte datafile (in .mat format), and scans the voxel with one of the three pulse sequences: Stejskal–Tanner pulsed-gradient spin echo (PGSE), second-order motion-compensated spin echo (MCSE), and monopolar simulated echo acquisition mode (STEAM).
<br />
<br />
Results will be processed and outputted in .mat file format for visualisation.
# Compilation
1. To build the monte carlo simulator
``` 
make
``` 
2. To execute code
``` 
./run_sim.exe
```
To change the diffusion configuration, please see the three input files in the input folder (montecarlofile, sequencefile, and substratefile)
# Sample

# Original paper
 https://doi.org/10.1002/mrm.27561

