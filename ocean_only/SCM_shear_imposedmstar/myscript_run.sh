#!/bin/bash
module load gcp
module load nco

echo "Model started:  " `date`
current_folder=$(basename "$PWD")
mkdir -p "INPUT"
mkdir -p "RESTART"
# Read values from files into arrays
readarray -t mstar_time < time_50_LES.txt 
readarray -t mstar_prescribed < mstar_50_LES.txt
run_executable="/gpfs/f5/gfdl_o/scratch/Qian.Xiao/MOM6-examples/build/ocean_only/MOM6"
epbl_dir="./results_mpre/"
y_loc='50'
mkdir -p "$epbl_dir"
# ======= Loop over timesteps =======
ntime=$(( ${#mstar_time[@]} - 1 ))

for (( ii=0; ii<ntime; ii++ )); do
    echo "==== Running timestep $ii ===="

    # === Modify input.nml ===
    if [[ $ii -eq 0 ]]; then
	sed -i "s|input_filename = '.*'|input_filename = 'n'|g" input.nml
    else
    	sed -i "s|input_filename = '.*'|input_filename = 'r'|g" input.nml
    fi

    # === Modify MOM_override ===
    sed -i "s|#override MSTAR =.*|#override MSTAR = ${mstar_prescribed[ii]} |" MOM_override
    sed -i "s|#override DAYMAX =.*|#override DAYMAX = ${mstar_time[ii+1]} |" MOM_override
    
    sync && sleep 1
    # === Run the model ===
    $run_executable || { echo " MOM6 run failed at timestep $ii" >&2; exit 1; }

    # === Move output files ===
    cp prog.nc "$epbl_dir/"
    mv "$epbl_dir/prog.nc" "$epbl_dir/prog_ep$(printf '%03d' $ii).nc"

    cp surffluxes.nc "$epbl_dir/"
    mv "$epbl_dir/surffluxes.nc" "$epbl_dir/surffluxes$(printf '%03d' $ii).nc"

    cp visc.nc "$epbl_dir/"
    mv "$epbl_dir/visc.nc" "$epbl_dir/visc_ep$(printf '%03d' $ii).nc"
    
    rm prog.nc surffluxes.nc visc.nc
    # === Copy restart file for next run ===
    cp RESTART/MOM.res.nc INPUT/
done

mkdir -p ./${y_loc}_results

ncrcat "$epbl_dir"/prog_ep*.nc "./${y_loc}_results/prog.nc"
ncrcat "$epbl_dir"/surffluxes*.nc "./${y_loc}_results/surffluxes.nc"
ncrcat "$epbl_dir"/visc_ep*.nc "./${y_loc}_results/visc.nc"

#rm "$epbl_dir"/prog_ep*.nc
#rm "$epbl_dir"/surffluxes*.nc
#rm "$epbl_dir"/visc_ep*.nc
echo "All timesteps completed successfully at $(date)"
