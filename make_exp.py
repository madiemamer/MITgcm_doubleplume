import numpy as np
import matplotlib.pyplot as plt
import run_config_funcs as rcf # type: ignore # import helpter functions
import platform
import os
import shutil
import params_dp as p

dir = '/storage/home/hcoda1/1/mmamer3/scratch/MITgcm'
dir = '/Users/mmamer3/Documents/MITgcm'
email = 'mmamer3@gatech.edu'
OSX = platform.system()
makeDirs = True

print('run_config is %s' %p.run_config)

if makeDirs:
    os.makedirs(p.run_config['run_dir'], exist_ok = True)

#========================================================================================
# Generate new experiment directory and copy over defaults

default_dir = os.path.join(p.run_config['exps_dir'], 'DEFAULT_doubleplume')
# create experimentary directory on SCRATCH and copy over default configuration
# NOTE: this step does not overwrite existing directories. 
if(makeDirs):
    run_subdir_list = ['build', 'code', 'input', 'run']
    for subdir in run_subdir_list:
        p.run_config['%s_dir'% subdir] = os.path.join(p.run_config['run_dir'], subdir)
        src_fpath = os.path.join(p.run_config['exps_dir'], 'DEFAULT_doubleplume', subdir)
        dst_dir = os.path.join(p.run_config['run_dir'], subdir)
        print(f"src: {src_fpath}, dst: {dst_dir}")
        shutil.copytree(src_fpath, dst_dir)

print('run directory and subdirectories:')
print(p.run_config['run_dir'])
print(os.listdir(p.run_config['run_dir']))

#==================#
# Model Parameters
#==================#

p.run_config['run_dir'] = '/Users/mmamer3/Documents/MITgcm/experiments/' + p.run_config['run_name']
if(makeDirs):
    data_params = [p.params01, p.params02, p.params03, p.params04, p.params05]
    #data file
    rcf.write_data(p.run_config, data_params, group_name='data', lf=p.run_config['lf'])
    #SIZE file
    rcf.createSIZEh(p.run_config, p.grid_params)

#======================#
# Diagnostic Parameters
#======================#


if(makeDirs):
    rcf.write_data(p.run_config, p.diag_params, group_name='diagnostics')
    
    ## create DIAGNOSTICS_SIZE.h
    Nlevels = p.grid_params['Nr']
    rcf.createDIAGSIZEh(p.run_config, p.Ndiags, p.numperlist, Nlevels)

    # create eedata
    rcf.create_eedata(p.run_config, p.grid_params['nTx'], p.grid_params['nTy'])

    #create data.obcs
    rcf.write_data(p.run_config, p.obcs_params, group_name='obcs')

#========================================================================================
# Copy Relevant files over
src = 'experiments/DEFAULT_doubleplume/plotPlumeData.py'
dst = p.run_config['run_dir']
shutil.copy2(src, dst)

src = 'experiments/DEFAULT_doubleplume/run_mitgcm.sh'
shutil.copy2(src,dst)

src = 'params_dp.py'
shutil.copy2(src,dst)
