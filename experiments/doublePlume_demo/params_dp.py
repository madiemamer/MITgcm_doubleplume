import numpy as np
import platform
import os

dir = '/storage/home/hcoda1/1/mmamer3/scratch/MITgcm'
dir = '/Users/mmamer3/Documents/MITgcm'
email = 'mmamer3@gatech.edu'
OSX = platform.system()
makeDirs = True
#========================================================================================
#main values to imput 

run_config = {}
run_config['ncpus_xy'] = [1,1] # cpu distribution in the x and y directions
run_config['run_name'] = 'doublePlume_demo'
run_config['ndays'] = 30 # simulation time (days)
run_config['test'] = False # if True, run_config['nyrs'] will be shortened to a few time steps

wallWidthInd = 2 #width of walls in units of dy
run_config['horiz_res_m'] = 500 # horizontal grid spacing (m)
run_config['Lx_m'] = 40000 # domain size in x (m)
run_config['Ly_m'] = 10000 + (2 * wallWidthInd * run_config['horiz_res_m']) # domain size in y (m) with walls (1 wall each side)
# NOTE: the number of grid points in x and y should be multiples of the number of cpus.
run_config['terminus_m'] = run_config['horiz_res_m']
terminus_index = int(run_config['terminus_m']/run_config['horiz_res_m'])

MITgcm_release = 'MITgcm-checkpoint68z' #Sept 2024 release
run_config['exps_dir'] = '/Users/mmamer3/Documents/MITgcm/experiments'
run_config['run_dir'] = os.path.join(run_config['exps_dir'], run_config['run_name'])
run_config['lf'] = '\r\n' # linebreak characters 

#==================#
# Domain Parameters
#==================#
domain_params = {}

domain_params['Lx'] = 10000 # m
domain_params['Ly'] = 40000 # m 
domain_params['H'] = 500 # m

#==================#
# Grid Domain
#==================#
grid_params = {}
grid_params['sNx'] = 50
grid_params['sNy'] = 11
grid_params['nTx'] = 1 # num of threads per processor in x-direction
grid_params['nTy'] = 1 # num of threads per processor in y-direction
grid_params['OLx'] = 3
grid_params['OLy'] = 3
grid_params['nSx'] = 1
grid_params['nSy'] = 1
grid_params['nPx'] = 1
grid_params['nPy'] = 1
grid_params['Nx'] = grid_params['sNx']*grid_params['nSx']*grid_params['nPx']
grid_params['Ny'] = grid_params['sNy']*grid_params['nSy']*grid_params['nPy']
grid_params['Nr'] = 50

for key, param  in grid_params.items():
    assert param%1==0, "grid parameter needs to be an integer"
    grid_params[key] = int(param)

# Grid Spacing
grid_params['delX'] = (domain_params['Lx']/grid_params['Nx'])*np.ones(grid_params['Nx'])
grid_params['delY'] = (domain_params['Ly']/grid_params['Ny'])*np.ones(grid_params['Ny'])
grid_params['delZ'] = (domain_params['H']/grid_params['Nr'])*np.ones(grid_params['Nr'])

#==================#
# Model Parameters
#==================#

# Continous Equation Parameters - Parameters 1
params01 = {}

params01['viscAz'] = 0.00001
params01['viscC2smag'] = 2.2
params01['diffKzT'] = 0.00001
params01['diffKhT'] = 20
params01['diffKzS'] = 0.00001
params01['diffKhS'] = 20

params01['no_slip_sides'] = False
params01['no_slip_bottom'] = True
params01['rigidLid'] = False
params01['implicitFreeSurface'] = True

params01['readBinaryPrec'] = 64

params01['staggerTimeStep'] = True

params01['tempAdvScheme'] = 33
params01['saltAdvScheme'] = 33

params01['eosType'] = 'JMD95Z'

params01['f0'] = 0.e-4
params01['beta'] = 0.E-11
params01['selectAddFluid'] = 1

params01['tRef'] = np.ones(grid_params['Nr'])*0. #ref temp
params01['sRef'] = np.ones(grid_params['Nr'])*34. #ref salt

# Elliptic Solver Parameters - Parameters 2
params02 = {}

params02['cg2dMaxIters'] = 300
params02['cg2dTargetResidual'] = 1.E-13
params02['cg3dMaxIters'] = 20
params02['cg3dTargetResidual'] = 1.E-8

# Time Stepping Parameters - Parameters 3
params03 = {}

params03['nIter0'] = 0
params03['endTime'] = 10000.0
params03['deltaT'] = 25.0
params03['abEps'] = 0.01
params03['pChkptFreq'] = 0.0
params03['chkptFreq'] = 0.0
params03['dumpFreq'] = 1000
params03['monitorSelect'] = 1
params03['monitorFreq'] = 1000.0
params03['periodicExternalForcing'] = False

# Gridding Parameters - Parameters 4
params04 = {}

params04['usingCartesianGrid'] = True

params04['delX'] = grid_params['delX']
params04['delY'] = grid_params['delY']
params04['delZ'] = grid_params['delZ']

# Input datasets
params05 = {}

params05['bathyFile'] = 'bathymetry.bin'
params05['hydrogThetaFile'] = 'tempini.bin'
params05['hydrogSaltFile'] = 'saltini.bin'

run_config['run_dir'] = '/Users/mmamer3/Documents/MITgcm/experiments/' + run_config['run_name']

#======================#
# Diagnostic Parameters
#======================#

# Time averaged fields
diag_fields_avg = [['THETA','SALT','UVEL','WVEL','VVEL'],
                    ['icefrtQi','icefrtMi','icefrtHi','icefrtBi','icefrtCi']
                    ]
diag_fields_max = 0
diag_fields_avg_name = ['dynDiag','plumeDiag']

numdiags_avg = len(diag_fields_avg)
numdiags_avg_total = 0
for i in range(len(diag_fields_avg)):
    numdiags_avg_total += len(diag_fields_avg[i])
diag_phase_avg = 0.0

run_config['inst_freq'] = 1 # in hours
run_config['tavg_freq'] = 1
diag_freq_inst = -params03['dumpFreq'] #-run_config['inst_freq']*3600
diag_freq_avg = params03['dumpFreq']


diag_params01 = {}
diag_params01['diag_mnc'] = False #<---you would need to modify this if you want netcdf output

for ii in range(numdiags_avg):  
    n = ii+1
    if len(diag_fields_avg[ii]) > diag_fields_max:
        diag_fields_max = len(diag_fields_avg[ii])
    diag_params01['fields(1:%i,%s)'%(len(diag_fields_avg[ii]),n)] ="','".join(diag_fields_avg[ii])
    diag_params01['fileName(%s)'%n] = diag_fields_avg_name[ii]
    diag_params01['frequency(%s)'%n] = diag_freq_avg
    diag_params01['timePhase(%s)'%n] = diag_phase_avg

# Instatneous fields
diag_fields_inst = [['THETA','SALT','UVEL','WVEL','VVEL'],
                    ['icefrtQi','icefrtMi','icefrtHi','icefrtBi','icefrtCi']
                    ]
diag_fields_names = ['dynDiag','plumeDiag']


numdiags_inst = len(diag_fields_inst)
numdiags_inst_total = 0
for i in range(len(diag_fields_inst)):
    numdiags_inst_total += len(diag_fields_inst[i])
diag_phase_inst = 0.0

for ii in range(numdiags_inst):
    n = numdiags_avg+ii+1
    if len(diag_fields_inst[ii]) > diag_fields_max:
        diag_fields_max = len(diag_fields_inst[ii])
    diag_params01['fields(1:%i,%s)'%(len(diag_fields_inst[ii]),n)] = "','".join(diag_fields_inst[ii])
    diag_params01['fileName(%s)'%n] = diag_fields_names[ii] + '_inst'
    diag_params01['frequency(%s)'%n] = diag_freq_inst
    diag_params01['timePhase(%s)'%n] = diag_phase_inst

Ndiags = n
numperlist = len(diag_fields_inst[1])

diag_params02 = {}
diag_params = [diag_params01, diag_params02]

#======================#
# Boundary Conditions
#======================#
obcs_params01 = {}

obcs_params01['useOBCSprescribe'] = True
obcs_params01['useOBCSsponge'] = True
obcs_params01['OB_Ieast'] = 11*-1
obcs_params01['OBEtFile'] = 'EBCt.bin'
obcs_params01['OBEsFile'] = 'EBCs.bin'
obcs_params01['OBEuFile'] = 'EBCu.bin'

obcs_params02 = {}

obcs_params03 = {}

obcs_params03['spongeThickness'] = 10
obcs_params03['Urelaxobcsbound'] = 10.0
obcs_params03['Urelaxobcsinner'] = 100000.0
obcs_params03['Vrelaxobcsbound'] = 10.0
obcs_params03['Vrelaxobcsinner'] = 100000.0

obcs_params = [obcs_params01, obcs_params02, obcs_params03]

#========================================================================================