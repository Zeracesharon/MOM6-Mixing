#!/bin/bash


echo "Model started:  " `date`
export FI_VERBS_PREFER_XRC=0
export LD_LIBRARY_PATH=/opt/cray/pe/netcdf/4.9.0.9/intel/2023.2/lib:$LD_LIBRARY_PATH
../../../FMS_Wave_Coupling/build/ncrc6.intel23/wave_ice_ocean/REPRO/MOM6
#../../../EqPac_SCM/GLS_SCM_oldMOM/MOM6_experiments/build/ncrc5.intel22/ice_ocean_SIS2/DEBUG/MOM6
