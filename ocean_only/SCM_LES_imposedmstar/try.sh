#!/bin/bash

module load gcp
module load nco
echo "Model started:  " `date`
current_folder=$(basename "$PWD")
dt=300
DT=300

sed -i "/^\"prog\",/c\\\"prog\",                ${dt},\"seconds\",1,\"days\",\"Time\"" diag_table
sed -i "/^\"surffluxes\",/c\\\"surffluxes\",          ${dt},\"seconds\",1,\"days\",\"Time\"" diag_table
sed -i "/^\"visc\",/c\\\"visc\",                ${dt},\"seconds\",1,\"days\",\"Time\"," diag_table

new_line="#override DT = ${DT}"
sed -i '/^#override DT =/c\'"$new_line" MOM_override
mkdir -p "INPUT"
mkdir -p "RESTART"
mkdir -p "results_mpre"
mkdir -p "results_mpre_noLT"
mkdir -p "LES_results"
mkdir -p "LES_noLT_results"
###modify MOM_override file
#y_loc=(300 200 150 130 110 90 70 50 30 10 0 -20 -40 -60 -80 -100 -120 -140 -200 -300)
y_loc=(50 30 10 0 -20 -40 -60 -80 -200)
y_tag=(50 30 10 0 S20 S40 S60 S80 S200)
NUM_TOT=${#y_loc[@]}
for (( i_loop=0; i_loop<NUM_TOT; i_loop++ )); do
    rm -rf RESTART/*     2>/dev/null || true 
    rm -rf INPUT/*       2>/dev/null || true
    rm -r prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc
     rm -f ./results_mpre/prog_ep*.nc
   rm -f ./results_mpre/surffluxes*.nc
   rm -f ./results_mpre/visc_ep*.nc
    rm -f ./results_mpre_noLT/prog_ep*.nc
   rm -f ./results_mpre_noLT/surffluxes*.nc
   rm -f ./results_mpre_noLT/visc_ep*.nc    
# get the y‐value for this iteration (will be empty if i >= length)
    y_km="${y_loc[i_loop]}"
    y_m=$(( y_km * 1000 ))
    y_tag_i=${y_tag[i_loop]}
    y=$(printf "%.1E" "$y_m") 
    if [[ -z "$y" ]]; then
      echo "WARNING: y_loc[$i_loop] is empty, skipping"
      continue
    fi
    mkdir -p "LES_results/dt_${dt}_DT_${DT}/${y_tag_i}_results"
    mkdir -p "LES_noLT_results/dt_${dt}_DT_${DT}/${y_tag_i}_results"
    echo "  ---- Iteration $i_loop: setting IDL_HURR_SCM_LOCY = $y ----"

    # build the replacement line
    new_line="#override IDL_HURR_SCM_LOCY = ${y}"
    sed -i '/^#override IDL_HURR_SCM_LOCY =/c\'"$new_line" MOM_override
   python run_LES.py ${y_tag_i} ${dt} ${DT} || { echo "❌ Python LES failed at $y_tag_i"; exit 1; }
   rm -f ./results_mpre/prog_ep*.nc
   rm -f ./results_mpre/surffluxes*.nc
   rm -f ./results_mpre/visc_ep*.nc 

   python run_LES_noLT.py ${y_tag_i} ${dt} ${DT} || { echo "❌ Python noLT failed at $y_tag_i"; exit 1; }
   rm -f ./results_mpre_noLT/prog_ep*.nc
   rm -f ./results_mpre_noLT/surffluxes*.nc
   rm -f ./results_mpre_noLT/visc_ep*.nc
   echo "All timesteps completed successfully at $(date) for single location"

done
cp -rf MOM_parameter_doc.all MOM_override MOM_input diag_table ./LES_results/dt_${dt}_DT_${DT}/
cp -rf MOM_parameter_doc.all MOM_override MOM_input diag_table ./LES_noLT_results/dt_${dt}_DT_${DT}/
gcp -r ./LES_noLT_results ./LES_results gfdl:/archive/Qian.Xiao/Qian.Xiao/MOM-examples/${current_folder}/postPro_results/NK_1200/KV/
echo "All done: `date`"

