#!/bin/bash

module load gcp
echo "Model started:  " `date`
current_folder=$(basename "$PWD")

NUM_TOT=9  ###overall loop number=y location

###modify MOM_override file
y_loc=(50 30 10 0 -20 -40 -60 -80 -200)
file_name=('Hurr05_046.nc' 'Hurr05_048.nc' 'Hurr05_050.nc' 'Hurr05_051.nc' 'Hurr05_053.nc' 'Hurr05_055.nc' 'Hurr05_057.nc' 'Hurr05_059.nc' 'Hurr05_071.nc')
NUM_TOT=${#y_loc[@]}

for (( i_loop=0; i_loop<NUM_TOT; i_loop++ )); do
    rm -rf RESTART/*     2>/dev/null || true 
    rm -r 20110401.*.nc prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc
  fn="${file_name[i_loop]}"
  (
    cd INPUT || { echo "ERROR: missing INPUT dir"; exit 1; }
    rm -f Hurr.nc
    ln -sf "../Forcing/${fn}" Hurr.nc
  )
    # get the y‐value for this iteration (will be empty if i >= length)
    y_km="${y_loc[i_loop]}"
    y_m=$(( y_km * 1000 ))
    y=$(printf "%.1E" "$y_m") 
    if [[ -z "$y" ]]; then
      echo "WARNING: y_loc[$i_loop] is empty, skipping"
      continue
    fi

    echo "  ---- Iteration $i_loop: setting IDL_HURR_SCM_LOCY = $y ----"

    ../../build/ice_ocean_SIS2/MOM6   \
        &> MOM6_${y_fmt}.log &   # redirect stdout/stderr if you like
    PID=$!
    
    echo "  → Started PID $PID for y=$y_fmt at $(date)"
    
    #  ▶ poll every 5s until it exits
    while kill -0 "$PID" 2>/dev/null; do
      echo "  → PID $PID still running… $(date)"
      sleep 5
    done
    
    #  ▶ now it’s gone—collect its exit code
    wait "$PID"
    ret=$?
    if [[ $ret -eq 0 ]]; then
      echo "Local run PID $PID finished SUCCESS at $(date)"
    else
      echo "Local run PID $PID finished FAILURE (exit $ret) at $(date)"
      exit 1
    fi
    
    
    
    mkdir results
    mv 20110401.ave_prog.nc ave_prog.nc
    mv 20110401.prog.nc prog.nc
    mv 20110401.surffluxes.nc surffluxes.nc
    mv 20110401.visc.nc visc.nc
    cp -rf prog.nc ocean.stats* surffluxes.nc visc.nc Vertical_coordinate.nc ave_prog.nc results/
    cp -rf MOM_IC.nc MOM_parameter_doc.all MOM_override MOM_input results/
    mv results "./${y_km}_results"
    if [[ -e "./${y_km}_results" ]]; then
       echo "Task $i_loop completed successfully with results exist. Continuing..."
    else
       echo "Error in results store for Task $i_loop. Stopping."
       exit 1
    fi

done
shopt -s nullglob

for entry in -*; do
  # only directories
  if [[ -d $entry ]]; then
    new="S${entry#-}"
    echo "Renaming '$entry' ?~F~R '$new'"
    mv -- "$entry" "$new"
  fi
done

echo "All done: `date`"

