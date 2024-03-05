This analysis is adapted and modified from aDDM simulations introduced by Krajbich et al. (2011).
To save time, follow these steps (in order) to install HDDM and run the file JMP_SimulationHDDM_Vero.ipynb if you're using Windows:

conda create --name yourenv python=3.6
conda activate yourenv
conda install conda-build
conda install pymc3
conda install pandas patsy matplotlib scipy tqdm
pip install cython
pip install hddm


The prerequisites are Anaconda Navigator - also, always check if the stuff is actually loaded in your environment
Taken from this website: https://groups.google.com/g/hddm-users/c/RpS9O9O7Fwc/m/P8ZyCpiqAwAJ?pli=1

To run the other files in order:

1. sim_addm_fun_ctrlValueDiff_v18.R
2. Simulation_ADDM_Vero.R
3. Create folders EXP1 - EXP8 - change number depending on how many differently-conditioned aDDMs you want to generate data for
4. CombineSims_Vero.R
5. JMP_SimulationHDDM_Vero_2.ipynb
6. prop_correct_resp_Model2.ipynb