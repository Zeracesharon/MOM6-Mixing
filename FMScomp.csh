#!/bin/csh
#

source build/intel/env

mkdir -p build/fms/

(cd build/fms/; rm -f path_names; \
../../src/mkmf/bin/list_paths -l ../../src/FMS; \
../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/ncrc-intel.mk -p libfms.a -c "-Duse_libMPI -Duse_netCDF -DHAVE_GETTID" path_names)
(cd build/fms/; make NETCDF=3 REPRO=1 libfms.a -j)
