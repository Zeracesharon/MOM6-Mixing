import os
import numpy as np
import netCDF4 as ncd
import time
import subprocess
import glob
import sys
import shutil
import math
os.system("module load nco")  # 有时必须在 shell 启动脚本中设置
if len(sys.argv) < 4:
    print("Usage: python script.py <y_loc> <dt> <DT>")
    sys.exit(1)
    
y_loc=sys.argv[1]
dt=sys.argv[2]
DT=sys.argv[3]
pwd=os.getcwd()  #'/Users/aakash/mom6_files/MOM6-examples-eqdisc/ocean_only/neural_network/EPBL/'
# epbl_dir=pwd+'/results_mpre_noLT/' ##you should modify this to your directory
# os.makedirs(epbl_dir, exist_ok=True)
epbl_dir = os.path.join(pwd, 'results_mpre_noLT')
if os.path.exists(epbl_dir):
    subprocess.run(['rm', '-rf', epbl_dir], check=True)
os.makedirs(epbl_dir, exist_ok=True)

combine_dir = os.path.join(os.getcwd(), f'LES_noLT_results/dt_{dt}_DT_{DT}/{y_loc}_results')
os.makedirs(combine_dir, exist_ok=True)

mstar_time = np.loadtxt(f'mstar_files/LES_noLT_results/dt_{dt}/time_{y_loc}_LES_noLT.txt')
mstar_prescribed = np.loadtxt(f'mstar_files/LES_noLT_results/dt_{dt}/mstar_{y_loc}_LES_noLT.txt')

# tx=tx4mtime(mstar_time,0.2,5)
ntime=len(mstar_time)-1
k=np.arange(0,ntime) #250)    #len(mstar_time)-1) overall timestep 
#k=np.arange(0,4) 
ii=0  # zeroth oteration

run_executable = '/gpfs/f5/gfdl_o/scratch/Qian.Xiao/MOM6-examples/build/ocean_only/MOM6'

### with ePBL:
for ii in k:
    
    # === Modify input.nml ===
    with open('input.nml', 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    for line in lines:
        if "input_filename =" in line:
            new_value = "'n'" if ii == 0 else "'r'"
            newline = f"         input_filename = {new_value}\n"
            new_lines.append(newline)
        else:
            new_lines.append(line)
    
    with open('input.nml', 'w') as f:
        f.writelines(new_lines)
    
    
    # === Modify MOM_override ===
    with open('MOM_override', 'r') as f:
        lines = f.readlines()
    
    new_lines = []
    for line in lines:
        if "#override MSTAR =" in line:
            newline = f"#override MSTAR = {mstar_prescribed[ii]} \n"
        elif "#override DAYMAX =" in line:
            newline = f"#override DAYMAX = {mstar_time[ii+1]} \n"
        else:
            newline = line
        new_lines.append(newline)
    
    with open('MOM_override', 'w') as f:
        f.writelines(new_lines)

    ### MOM_override modified!

    
    os.system(run_executable)
    # List of expected output files
    output_files = {
        'prog.nc': f'prog_ep{ii:04d}.nc',
        'surffluxes.nc': f'surffluxes{ii:04d}.nc',
        'visc.nc': f'visc_ep{ii:04d}.nc',
    }
    
    # Move and rename only if output file exists
    for src, renamed in output_files.items():
        if os.path.exists(src):
            shutil.move(src, os.path.join(epbl_dir, renamed))
            print(f"Saved and removed: {src} → {renamed}")
        else:
            print(f"Skipping: {src} not found for step {ii}")
    
    # Always copy restart for next step
    if os.path.exists('RESTART/MOM.res.nc'):
        shutil.copy('RESTART/MOM.res.nc', 'INPUT/')

    # prog_src = 'prog.nc'
    # surf_src = 'surffluxes.nc'
    # visc_src = 'visc.nc'
    
    # prog_dst = os.path.join(epbl_dir, f'prog_ep{ii:03d}.nc')
    # surf_dst = os.path.join(epbl_dir, f'surffluxes{ii:03d}.nc')
    # visc_dst = os.path.join(epbl_dir, f'visc_ep{ii:03d}.nc')
    
    # # Copy and rename
    # shutil.copy(prog_src, prog_dst)
    # shutil.copy(surf_src, surf_dst)
    # shutil.copy(visc_src, visc_dst)
    
    # # Restart file
    # shutil.copy('RESTART/MOM.res.nc', 'INPUT/')
def batch_ncrcat(input_dir, pattern, output_file, batch_size=500):
    tmp_files = []
    files = sorted(glob.glob(os.path.join(input_dir, pattern)))

    n_batches = math.ceil(len(files) / batch_size)

    for i in range(n_batches):
        batch_files = files[i*batch_size : (i+1)*batch_size]
        tmp_out = os.path.join(input_dir, f"tmp_{i:04d}.nc")
        subprocess.run(["ncrcat", "-O"] + batch_files + [tmp_out], check=True)
        tmp_files.append(tmp_out)
    
    subprocess.run(["ncrcat", "-O"] + tmp_files + [output_file], check=True)

    # Optional: clean up intermediate files
    for tmp in tmp_files:
        os.remove(tmp)

# Combine in batches
batch_ncrcat(epbl_dir, "prog_ep*.nc", os.path.join(combine_dir, "prog.nc"))
batch_ncrcat(epbl_dir, "surffluxes*.nc", os.path.join(combine_dir, "surffluxes.nc"))
batch_ncrcat(epbl_dir, "visc_ep*.nc", os.path.join(combine_dir, "visc.nc"))




