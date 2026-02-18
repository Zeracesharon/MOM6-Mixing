#!/bin/csh
set platform="ncrc5.intel23"
set rootdir = `dirname $0`
source ${rootdir}/envs.${platform}
mkdir -p build/ice_ocean_SIS2/
(cd build/ice_ocean_SIS2/; rm -f path_names; \
../../src/mkmf/bin/list_paths -l ./ ../../src/MOM6/config_src/{infra/FMS1,memory/dynamic_symmetric,drivers/FMS_cap,external} ../../src/SIS2/config_src/dynamic_symmetric ../../src/MOM6/src/{*,*/*}/ ../../src/{atmos_null,coupler,land_null,ice_param,icebergs/src,SIS2,FMS/coupler,FMS/include}/)
(cd build/ice_ocean_SIS2/; \
../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/ncrc-intel.mk -o '-I../fms' -p MOM6 -l '-L../fms -lfms' -c '-Duse_AM3_physics -D_USE_LEGACY_LAND_' path_names )
(cd build/ice_ocean_SIS2/; make NETCDF=3 REPRO=1 MOM6 -j)
