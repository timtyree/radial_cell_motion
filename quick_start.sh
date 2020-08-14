#!/bin/bash
# quick start for pde simulation
# Programmer: Tim Tyree
# Rappel Group, Univ. of Calif., San Diego
# Date:  7.8.2020
# Description: downloads github repository, installs dependencies, 
# integrates a radial cAMP profile for 20 minutes, and then 
# prints the final cAMP radial field. 

#initialize (can be redundant as written)
git clone https://github.com/timtyree/radial_cell_motion.git
brew link --overwrite python@3.8
pip3 install numpy scipy pandas
cd radial_cell_motion/pde-sim/pde-sim-transfer

python3 pde_sim_returns_camp.py 0.02 100.0 4.0 10.0 2.0 0.005 10.0
# where the parameters are the following
# - kPDE = 0.02/sec, 
# - LPDE = 100µm^2/sec, 
# - c0 = 4.0nM, 
# - T=10.0min, 
# - iter_no = 2.0 # the number of pulses to integrate (tmax = iter_no*T)
# - dt = 0.005sec
# - time_res = 10.0min # the frequency data is printed.