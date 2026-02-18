#!/bin/csh
# Sample OceanOnlyComp.csh . intel18 repro

source build/intel/env

mkdir -p build/ocean_only/
(cd build/ocean_only/; rm -f path_names; \
../../src/mkmf/bin/list_paths -l ./ ../../src/MOM6/{config_src/infra/FMS1,config_src/memory/dynamic_symmetric,config_src/drivers/solo_driver,config_src/external,src/{*,*/*}}/ ; \
../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/ncrc-intel.mk -o "-I../fms/ -I../../src/MOM6/src/user" -p MOM6 -l "-L../fms -lfms"  path_names)
(cd build/ocean_only/; make NETCDF=3 REPRO=1 MOM6 -j)
