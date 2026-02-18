#!/bin/csh
#

set version=REPRO
set platform="ncrc5.intel23"
set rootdir = `dirname $0`
source ${rootdir}/envs.${platform}

mkdir -p build/fms/

(cd build/fms/; rm -f path_names; \
../../src/mkmf/bin/list_paths -l ../../src/FMS; \
../../src/mkmf/bin/mkmf -t ../../src/mkmf/templates/ncrc-intel.mk -p libfms.a -c "-DUSE_FMS2_IO -Duse_libMPI -Duse_netCDF -DHAVE_GETTID" path_names)
(cd build/fms/; make NETCDF=3 $version=1 libfms.a -j)
