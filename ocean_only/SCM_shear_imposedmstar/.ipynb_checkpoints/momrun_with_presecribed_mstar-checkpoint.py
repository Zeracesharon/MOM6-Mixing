import os
import numpy as np
import netCDF4 as ncd
import time
import subprocess
import glob


os.system("module load nco")  # 有时必须在 shell 启动脚本中设置
y_loc='50'
pwd=os.getcwd()  #'/Users/aakash/mom6_files/MOM6-examples-eqdisc/ocean_only/neural_network/EPBL/'
epbl_dir=pwd+'/results_mpre/' ##you should modify this to your directory
combine_dir = os.path.join(os.getcwd(), y_loc + '_results/')
os.makedirs(combine_dir, exist_ok=True)

mstar_time=np.loadtxt('time_50_LES.txt')
mstar_prescribed=np.loadtxt('mstar_50_LES.txt')

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

    os.system('cp prog.nc '+epbl_dir)
    os.system('mv '+epbl_dir+'prog.nc '+epbl_dir+'prog_ep'+str(ii).zfill(3)+'.nc')
    
    os.system('cp surffluxes.nc '+epbl_dir)
    os.system('mv '+epbl_dir+'surffluxes.nc '+epbl_dir+'surffluxes'+str(ii).zfill(3)+'.nc')

    os.system('cp visc.nc '+epbl_dir)
    os.system('mv '+epbl_dir+'visc.nc '+epbl_dir+'visc_ep'+str(ii).zfill(3)+'.nc')

    os.system('cp RESTART/MOM.res.nc INPUT/')


subprocess.run(f"ncrcat {epbl_dir}/prog_ep*.nc {combine_dir}/prog.nc", shell=True, check=True)
subprocess.run(f"ncrcat {epbl_dir}/surffluxes*.nc {combine_dir}/surffluxes.nc", shell=True, check=True)
subprocess.run(f"ncrcat {epbl_dir}/visc_ep*.nc {combine_dir}/visc.nc", shell=True, check=True)



