Repository Overview

This repository adds two vertical mixing schemes on top of an ocean–ice coupling model. It is built on a more up-to-date MOM6 core than the repository FMS_Wave_Coupling_Mixing, but it does not couple with the WW3 wave model. Therefore, it does not support two-way ocean–wave coupling modes. To consider Langmuir turbulence (LT) effects, this repo supports one-way coupling and theoretical wave options only.

Key difference vs. FMS_Wave_Coupling_Mixing:

The main difference lies in the MOM6 core version. Here we incorporate the two vertical mixing schemes into a more updated MOM6.

All parameter settings, diag_table configurations for diagnosed quantities, and workflows are intended to be the same as FMS_Wave_Coupling_Mixing:
https://github.com/Zeracesharon/FMS_Wave_Coupling_Mixing/tree/main

Build & Compilation

You can compile the source (src/) with either:

an ocean-only driver, or

an ocean–ice (MOM6–SIS2) driver.

The compilation scripts are different for the two drivers.

Ocean-only driver

Includes built-in surface forcing modules (e.g., idealized hurricane surface forcing).

Recommended when you want to use internal surface-forcing modules directly.

Compile (we normally use FMS1):

./FMScomp.csh
./OceanOnly.csh

Ocean–ice (MOM6–SIS2) driver

The hurricane-forcing module is not compiled in the ocean–ice build.

You therefore need wind-forcing input files (external) to run hurricane SCM cases.

You can choose the compiling file to compile FMS1 or FMS2.
Important: if you choose FMS1/FMS2, go to src/ and link the corresponding version to FMS.

Build order (ocean–ice):

Compile FMS first

Compile MOM6–SIS2

If you choose FMS1
./FMScomp_iceFMS1.csh
./IceOceancomp_example_FMS1.csh

If you choose FMS2 (recommended for legacy I/O compatibility)

We normally use FMS2 because it supports some older data read types to accommodate old input files.

./FMScomp_ice.csh
./IceOceancomp_example.csh

Notes on FMS2 (Important for Successful Compilation)

The src/FMS2 source code in this repository may be slightly different from the version you obtain when you git --recursive download from NOAA-GFDL/MOM6-examples/. The reason is that the FMS2 link in their source code may not point to the specific version required to compile in this repo.

Related issue:

NOAA-GFDL/MOM6-examples/issues: REGISTER_UNLIMITED_COMPRESSED_AXIS does not exist or is not accessible with FMS2 compiler #565
https://github.com/NOAA-GFDL/MOM6-examples/issues/565

If that issue is closed/unavailable, the core idea and fix are summarized below.

FMS2 Compilation Failure and Fix (register_unlimited_compressed_axis)

When compiling with FMS2 and -DUSE_FMS2_IO, the file icebergs_fms2io.F90 requires the symbol:

register_unlimited_compressed_axis

However, in some FMS2 I/O source versions, this subroutine is missing (both declaration and implementation), causing compilation failure.

Proposed solution

You will need to:

1) Edit src/FMS2/fms2_io/fms2_io.F90

Add near ~line 55:

public :: register_unlimited_compressed_axis

2) Edit src/FMS2/fms2_io/netcdf_io.F90

Add:

public :: register_unlimited_compressed_axis


Then add the following implementation:

subroutine register_unlimited_compressed_axis(fileobj, dimension_name, dimension_length)
  class(FmsNetcdfFile_t), intent(inout) :: fileobj !< File object.
  character(len=*), intent(in) :: dimension_name   !< Dimension name.
  integer, intent(in) :: dimension_length          !< Dimension length for the current rank

  integer, dimension(:), allocatable :: npes_start !< Starting index for each PE
  integer, dimension(:), allocatable :: npes_count !< Local size for each PE
  integer :: i
  integer :: err
  integer :: dimid

  ! Gather all local dimension lengths on the I/O root PE.
  allocate(npes_start(size(fileobj%pelist)))
  allocate(npes_count(size(fileobj%pelist)))

  call mpp_gather((/dimension_length/), npes_count, pelist=fileobj%pelist)

  npes_start(1) = 1
  do i = 1, size(fileobj%pelist)-1
     npes_start(i+1) = npes_start(i) + npes_count(i)
  enddo

  call append_compressed_dimension(fileobj, dimension_name, npes_start, npes_count)

  if (fileobj%is_root .and. .not. fileobj%is_readonly) then
     call set_netcdf_mode(fileobj%ncid, define_mode)
     err = nf90_def_dim(fileobj%ncid, trim(dimension_name), unlimited, dimid)
     call check_netcdf_code(err, "Netcdf_add_dimension: file:"//trim(fileobj%path)//" dimension name:"// &
                           trim(dimension_name))
  endif

end subroutine register_unlimited_compressed_axis


After making these changes:

recompile FMS with FMS linked to FMS2 (instead of FMS1),

then compile the ocean–ice driver (MOM6–SIS2).

This enables compiling ocean_SIS2 with FMS2 and allows FMS2 I/O mode required for certain override/input data types.

Example Folders

After compilation, you will find two example folders corresponding to your solver/driver:

ocean_only/ → for ocean-only builds

ocean_ice/ → for ocean–ice (MOM6–SIS2) builds

Ocean-only Examples

The ocean_only/ folder mainly demonstrates how to:

impose/replace mstar in ePBL using input files

postprocess results

configure idealized hurricane surface forcing

configure constant wind forcing (taux = constant, tauy = 0 or vice versa)

(1) SCM_LES_imposedmstar
Goal

Impose mstar using input files. The input mstar can come from:

LES results, or

shear-diagnosed energy (offline hybrid approach coupling κ-shear with ePBL)

Core idea

For each location, in MOM_override after setting:

#override IDL_HURR_SCM_LOCY = 1.0E+04


set:

#override EPBL_MSTAR_SCHEME = "CONSTANT"


instead of "REICHL_H18".

Then at each time step you:

update MSTAR to the prescribed value (from your input file), and

update DAYMAX to the newest step time until the end.

This means your input file must contain (for SCM, per location):

time

mstar value

Outputs and file management

After each time step, you need to save and rename output files (so they won’t be overwritten):

prog.nc

surffluxes.nc

visc.nc

However, outputs may not be written at every time step if your diag_table output frequency is larger than your timestep.

If there is no new output file at that step, simply do not rename/save for that step.

Output filenames will depend on the model time when they were produced; they do not need to be consecutive.

Restart handling (why this workflow exists)

You must also save each time step’s restart file and put it under INPUT/ so the next time step restarts from the last ocean state.

We use this workflow because:

MOM6 only reads MOM_override at initialization (start of each run), and

we set restart file generation only when the simulation ends,
so to impose time-varying mstar we need to run step-by-step.

Practical steps (mstar imposed run)
Step 1: Prepare mstar input file

We provide examples generating mstar from LES data:

LES_mstar/mstar_LES_generate.ipynb for 5 m/s cases

LES_mstar/mstar_LES_generate_10mps.ipynb for 10 m/s cases

Main differences between notebooks:

LES file path

case number name

Definition of mstar in notebooks:
mstar is obtained by integrating all negative buoyancy fluxes within the mixed-layer depth, which is diagnosed by the depth of maximum temperature gradient (also stated in Xiao & Reichl, 2026).

The results are interpolated to a fixed output time frequency, normally set to match the MOM6 timestep.

Important note on time steps:
You can set the mstar file output dt larger than MOM6 dt, meaning several MOM6 steps use the same mstar value, but it must be an integer multiple of MOM6 dt to avoid time misalignment.

We provide LES data with:

no Langmuir turbulence: LES_noLT_results

with Langmuir turbulence: LES_results

Step 2: Place generated files

Put mstar files under:

mstar_files/

Inside mstar_files/, there are subfolders indicating different output frequencies. Example:

mstar_files/LES_noLT_results/dt_100/mstar_0_LES_noLT.txt

Meaning:

LES_noLT case

mstar output interval = 100 s

location/case label = 0 (0 km away from storm center on the right side of track)

S20 means 20 km away from center on the left side of track

Step 3: Run scripts

Scripts:

run_LES_noLT.py

run_LES.py

You can use them after changing:

directory settings (mstar_time, mstar_prescribed)

paths to executables (run_executable)

output combine directory, e.g.:

combine_dir = os.path.join(os.getcwd(), f'LES_noLT_results/dt_{dt}_DT_{DT}/{y_loc}_results')


to ensure the correct output directory is used.

Step 4: Run single location or multiple locations

Single location: use singlerun.sh (runs one location for all time steps)

Multiple locations: use try.sh
Example command:

nohup bash try.sh > submit_jobs.log 2>&1 &


In try.sh:

dt = timestep for the input mstar file (LES mstar file)

DT = timestep used by MOM6

Before running, create temporary output directories:

results_mpre/

results_mpre_noLT/

After runs:

individual per-step outputs are combined into one dataset per location

combined outputs are stored under:

LES_noLT_results/

LES_results/

Example output folder:

LES_noLT_results/dt_100_DT_50/

Inside, you will find subfolders corresponding to each location; these are final outputs for postprocessing.

Optional experiments

You can also play with vertical mixing schemes by changing MOM_override settings, e.g.:

turn κ-shear on/off

adjust MIX_LEN_EXPONENT

(2) SCM_shear_imposedmstar
Goal

Generate mstar input files from shear-only runs (pre-run κ-shear results), instead of LES.

For the pre-run shear case:

only κ-shear is on

other vertical mixing schemes are off

Required physical variables:

temperature profile

Tflx_dia_diff (temperature fluxes)

taux, tauy

Output frequency should be at least equal or smaller than desired mstar file output frequency.

Generate mstar from shear results

Use:

mstar_LES_generate_shear_5mps.ipynb

mstar_LES_generate_shear_10mps.ipynb

Since shear cases do not include LT effects:

only noLT results are simulated by default.

If you want LT effects:

set USE_WAVES=TRUE

choose theoretical waves or wave-forcing input files
This adds mstar_LT on top of the constant mstar value from the file.

Clarification of the three time steps

There are three different time steps:

time step used for mstar file generation

time step used for MOM6 running (DT)

time step for MOM6 ocean outputs (set in diag_table)

Recommendations:

use identical time step for mstar generation and MOM6 model steps, OR

run MOM6 with smaller dt while keeping mstar dt as an integer multiple

For MOM6 ocean output files:

you can set any output interval as long as it is ≥ MOM6 model dt

typically we set 1800 s

(3) SCM_wind folder

This folder shows how to set constant wind forcing using the surface wind forcing module:

#override WIND_CONFIG = "const"
#override CONST_WIND_TAUX = 0.1
#override CONST_WIND_TAUY = 0.0
#override VARIABLE_WINDS = False

LES Data Naming Rules (MLD = 32 m cases)

The naming rule applies similarly for MLD = 10 m cases.

Location indexing:

stations arranged as 1 (500 km north) to 101 (500 km south)

interval: 10 km per station

Storm translation:

storm moves east → west

< 51: storm-relative right of track

> 51: storm-relative left of track

5 m/s translation speed cases
With Langmuir

TC013 LOC30: 210 km / 200 km north

TC011 LOC26: 250 km

TC008 LOC020: 310 / 300 km

TC004 LOC012: 390 / 400 km

TC002 LOC008: 430 km

TC021 DATA055065/LOC046 LT # 500-(46-1)*10 = 50 km north (right side 50 km)

TC022 DATA055065/LOC048 LT # 30 km north

TC023 DATA055065/LOC050 LT # 10 km north

TC024 DATA055065/LOC051 LT # 0 km north

TC025 DATA055065/LOC053 LT # -20 km

TC026 DATA055065/LOC055 LT # -40 km

TC027 DATA055065/LOC057 LT # -60 km

TC028 DATA055065/LOC059 LT # -80 km

TC029 DATA055065/LOC061 LT # -100 km

TC030 DATA055065/LOC071 LT # -200 km

No Langmuir

TC031 DATA055065/LOC046 no LT

TC032 DATA055065/LOC048 no LT

TC033 DATA055065/LOC050 no LT

TC034 DATA055065/LOC051 no LT

TC035 DATA055065/LOC053 no LT

TC036 DATA055065/LOC055 no LT

TC037 DATA055065/LOC057 no LT

TC038 DATA055065/LOC059 no LT

TC039 DATA055065/LOC061 no LT

TC040 DATA055065/LOC071 no LT

10 m/s translation speed cases
No Langmuir

TC041 DATA105065/LOC046 no LT

TC042 DATA100565/LOC048 no LT

TC043 DATA105065/LOC050 no LT

TC044 DATA105065/LOC051 no LT

TC045 DATA105065/LOC053 no LT

TC046 DATA105065/LOC055 no LT

TC047 DATA105065/LOC057 no LT

TC048 DATA105065/LOC059 no LT

TC049 DATA105065/LOC061 no LT

TC050 DATA105065/LOC071 no LT

TC051 DATA105065/LOC081 no LT

With Langmuir

TC052 DATA105065/LOC046 LT

TC053 DATA105065/LOC048 LT

TC054 DATA105065/LOC050 LT

TC055 DATA105065/LOC051 LT

TC056 DATA105065/LOC053 LT

TC057 DATA105065/LOC055 LT

TC058 DATA105065/LOC057 LT

TC059 DATA105065/LOC059 LT

TC060 DATA105065/LOC061 LT

TC061 DATA105065/LOC071 LT

TC062 DATA105065/LOC081 LT

Ocean–ice Examples

The ocean_ice/ folder provides a case demonstrating how to set up the idealized hurricane SCM case with wind–wave forcing input files.

Run:

try.sh or

try_10mps.sh

This generates results at 9 locations for:

mixed-layer depth = 32 m

translation speed = 5 m/s and 10 m/s
These can be compared with LES results.

Reference

Xiao, Q. and Reichl, B.G., 2026.
A hybrid $\kappa$-shear–ePBL vertical mixing scheme with a Langmuir turbulence parameterization for tropical cyclones in a coupled ocean-wave model.
