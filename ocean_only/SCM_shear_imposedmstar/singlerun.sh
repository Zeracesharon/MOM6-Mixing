#!/bin/bash
module load gcp
module load nco

echo "Model started:  " `date`
current_folder=$(basename "$PWD")
mkdir -p "INPUT"
mkdir -p "RESTART"
# Read values from files into arrays
y_loc='50'
mkdir -p "results_mpre"
mkdir -p "LES_results"
mkdir -p "LES_noLT_results"
mkdir -p "LES_results/${y_loc}_results"
mkdir -p "LES_noLT_results/${y_loc}_results"

python run_LES.py ${y_loc}
rm ./results_mpre/prog_ep*.nc
rm ./results_mpre/surffluxes*.nc
rm ./results_mpre/visc_ep*.nc

python run_LES_noLT.py ${y_loc}

rm ./results_mpre_noLT/prog_ep*.nc
rm ./results_mpre_noLT/surffluxes*.nc
rm ./results_mpre_noLT/visc_ep*.nc
echo "All timesteps completed successfully at $(date)"
