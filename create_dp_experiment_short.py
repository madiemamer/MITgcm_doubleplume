import numpy as np
import matplotlib.pyplot as plt
import run_config_funcs as rcf # type: ignore # import helpter functions
import platform
import os
import shutil

exps_dir = '/storage/home/hcoda1/1/mmamer3/scratch/MITgcm/experiments/'
exps_dir = '/Users/mmamer3/Documents/MITgcm/experiments/'
email = 'mmamer3@gatech.edu'
OSX = platform.system()
makeDirs = True
run_name = 'doublePlume_demo'

#========================================================================================
#main values to imput 

run_dir = exps_dir + run_name
print(f'run_dir is {run_dir}')

if makeDirs:
    os.makedirs(run_dir, exist_ok = True)

#========================================================================================
# Generate new experiment directory and copy over defaults

default_dir = exps_dir + 'DEFAULT_doubleplume'
# create experimentary directory on SCRATCH and copy over default configuration
# NOTE: this step does not overwrite existing directories. 
if(makeDirs):
    subDirs = ['build', 'code', 'input', 'run']
    for subDir in subDirs:

        shutil.copytree(default_dir + '/' + subDir, run_dir + '/' + subDir)

print('run directory and subdirectories:')
print(os.listdir(run_dir))



# #==================#
# # Domain Parameters
# #==================#
# domain_params = {}

# domain_params['Lx'] = 10000 # m
# domain_params['Ly'] = 40000 # m 
# domain_params['H'] = 500 # m

# #==================#
# # Grid Domain
# #==================#
# grid_params = {}
# grid_params['sNx'] = 50
# grid_params['sNy'] = 11
# grid_params['nTx'] = 1 # num of threads per processor in x-direction
# grid_params['nTy'] = 1 # num of threads per processor in y-direction
# grid_params['OLx'] = 3
# grid_params['OLy'] = 3
# grid_params['nSx'] = 1
# grid_params['nSy'] = 1
# grid_params['nPx'] = 1
# grid_params['nPy'] = 1
# grid_params['Nx'] = grid_params['sNx']*grid_params['nSx']*grid_params['nPx']
# grid_params['Ny'] = grid_params['sNy']*grid_params['nSy']*grid_params['nPy']
# grid_params['Nr'] = 50

# for key, param  in grid_params.items():
#     assert param%1==0, "grid parameter needs to be an integer"
#     grid_params[key] = int(param)

# # Grid Spacing
# grid_params['delX'] = (domain_params['Lx']/grid_params['Nx'])*np.ones(grid_params['Nx'])
# grid_params['delY'] = (domain_params['Ly']/grid_params['Ny'])*np.ones(grid_params['Ny'])
# grid_params['delZ'] = (domain_params['H']/grid_params['Nr'])*np.ones(grid_params['Nr'])

# #==================#
# # Model Parameters
# #==================#

# # Continous Equation Parameters - Parameters 1
# params01 = {}

# params01['viscAz'] = 0.00001
# params01['viscC2smag'] = 2.2
# params01['diffKzT'] = 0.00001
# params01['diffKhT'] = 20
# params01['diffKzS'] = 0.00001
# params01['diffKhS'] = 20

# params01['no_slip_sides'] = False
# params01['no_slip_bottom'] = True
# params01['rigidLid'] = False
# params01['implicitFreeSurface'] = True

# params01['readBinaryPrec'] = 64

# params01['staggerTimeStep'] = True

# params01['tempAdvScheme'] = 33
# params01['saltAdvScheme'] = 33

# params01['eosType'] = 'JMD95Z'

# params01['f0'] = 0.e-4
# params01['beta'] = 0.E-11
# params01['selectAddFluid'] = 1

# params01['tRef'] = np.ones(grid_params['Nr'])*0. #ref temp
# params01['sRef'] = np.ones(grid_params['Nr'])*34. #ref salt

# # Elliptic Solver Parameters - Parameters 2
# params02 = {}

# params02['cg2dMaxIters'] = 300
# params02['cg2dTargetResidual'] = 1.E-13
# params02['cg3dMaxIters'] = 20
# params02['cg3dTargetResidual'] = 1.E-8

# # Time Stepping Parameters - Parameters 3
# params03 = {}

# params03['nIter0'] = 0
# params03['endTime'] = 10000.0
# params03['deltaT'] = 25.0
# params03['abEps'] = 0.01
# params03['pChkptFreq'] = 0.0
# params03['chkptFreq'] = 0.0
# params03['dumpFreq'] = 1000
# params03['monitorSelect'] = 1
# params03['monitorFreq'] = 1000.0
# params03['periodicExternalForcing'] = False

# # Gridding Parameters - Parameters 4
# params04 = {}

# params04['usingCartesianGrid'] = True

# delx = np.asarray(grid_params['delX'], dtype = '>f8')
# dely = np.asarray(grid_params['delY'], dtype = '>f8')
# delz = np.asarray(grid_params['delZ'], dtype = '>f8')

# run_config['run_dir'] = run_config['exps_dir'] + '/' + run_config['run_name']
# delx.tofile(run_config['run_dir'] + '/input/' + 'delx.bin')
# dely.tofile(run_config['run_dir'] + '/input/' + 'dely.bin')
# delz.tofile(run_config['run_dir'] + '/input/' + 'delz.bin')

# params04['delXfile'] = 'delX.bin'
# params04['delYfile'] = 'delY.bin'
# params04['delZ'] = grid_params['delZ'][0]

# # Input datasets
# params05 = {}

# params05['bathyFile'] = 'bathymetry.bin'
# params05['hydrogThetaFile'] = 'tempini.bin'
# params05['hydrogSaltFile'] = 'saltini.bin'

# run_config['run_dir'] = '/Users/mmamer3/Documents/MITgcm/experiments/' + run_config['run_name']
# if(makeDirs):
#     data_params = [params01, params02, params03, params04, params05]
#     #data file
#     rcf.write_data(run_config, data_params, group_name='data', lf=run_config['lf'])
#     #SIZE file
#     rcf.createSIZEh(run_config, grid_params)

# #======================#
# # Diagnostic Parameters
# #======================#

# # Time averaged fields
# diag_fields_avg = [['THETA','SALT','UVEL','WVEL','VVEL'],
#                     ['icefrtQi','icefrtMi','icefrtHi','icefrtBi','icefrtCi', \
#                      'icefrtQo', 'icefrtMo', 'icefrtHo', 'icefrtBo'],
#                     ]
# diag_fields_max = 0
# diag_fields_avg_name = ['dynDiag','plumeDiag']

# numdiags_avg = len(diag_fields_avg)
# numdiags_avg_total = 0
# for i in range(len(diag_fields_avg)):
#     numdiags_avg_total += len(diag_fields_avg[i])
# diag_phase_avg = 0.0

# run_config['inst_freq'] = 12
# run_config['tavg_freq'] = 12
# diag_freq_inst = -run_config['inst_freq']*3600
# diag_freq_avg = run_config['tavg_freq']*3600

# diag_params01 = {}
# diag_params01['diag_mnc'] = False #<---you would need to modify this if you want netcdf output

# for ii in range(numdiags_avg):  
#     n = ii+1
#     if len(diag_fields_avg[ii]) > diag_fields_max:
#         diag_fields_max = len(diag_fields_avg[ii])
#     diag_params01['fields(1:%i,%s)'%(len(diag_fields_avg[ii]),n)] ="','".join(diag_fields_avg[ii])
#     diag_params01['fileName(%s)'%n] = diag_fields_avg_name[ii]
#     diag_params01['frequency(%s)'%n] = diag_freq_avg
#     diag_params01['timePhase(%s)'%n] = diag_phase_avg

# # Instatneous fields
# diag_fields_inst = [['THETA','SALT','UVEL','WVEL','VVEL'],
#                     ['icefrtQi','icefrtMi','icefrtHi','icefrtBi','icefrtCi', \
#                      'icefrtQo', 'icefrtMo', 'icefrtHo', 'icefrtBo']
#                     ]
# diag_fields_names = ['dynDiag','plumeDiag']


# numdiags_inst = len(diag_fields_inst)
# numdiags_inst_total = 0
# for i in range(len(diag_fields_inst)):
#     numdiags_inst_total += len(diag_fields_inst[i])
# diag_phase_inst = 0.0

# for ii in range(numdiags_inst):
#     n = numdiags_avg+ii+1
#     if len(diag_fields_inst[ii]) > diag_fields_max:
#         diag_fields_max = len(diag_fields_inst[ii])
#     diag_params01['fields(1:%i,%s)'%(len(diag_fields_inst[ii]),n)] = "','".join(diag_fields_inst[ii])
#     diag_params01['fileName(%s)'%n] = diag_fields_names[ii] + '_inst'
#     diag_params01['frequency(%s)'%n] = diag_freq_inst
#     diag_params01['timePhase(%s)'%n] = diag_phase_inst

# Ndiags = n

# diag_params02 = {}
# diag_params = [diag_params01, diag_params02]

# #======================#
# # Boundary Conditions
# #======================#
# obcs_params01 = {}

# obcs_params01['useOBCSprescribe'] = True
# obcs_params01['useOBCSsponge'] = True
# obcs_params01['OB_Ieast'] = 11*-1
# obcs_params01['OBEtFile'] = 'EBCt.bin'
# obcs_params01['OBEsFile'] = 'EBCs.bin'
# obcs_params01['OBEuFile'] = 'EBCu.bin'

# obcs_params03 = {}

# obcs_params03['spongeThickness'] = 10
# obcs_params03['Urelaxobcsbound'] = 10.0
# obcs_params03['Urelaxobcsinner'] = 100000.0
# obcs_params03['Vrelaxobcsbound'] = 10.0
# obcs_params03['Vrelaxobcsinner'] = 100000.0

# obcs_params = [obcs_params01, obcs_params03]

# if(makeDirs):
#     rcf.write_data(run_config, diag_params, group_name='diagnostics')
    
#     ## create DIAGNOSTICS_SIZE.h
#     Nlevels = grid_params['Nr']
#     rcf.createDIAGSIZEh(run_config, Ndiags, Nlevels)

#     # create eedata
#     rcf.create_eedata(run_config, grid_params['nTx'], grid_params['nTy'])

#     #create data.obcs
#     rcf.write_data(run_config, obcs_params, group_name='obcs')

# #========================================================================================