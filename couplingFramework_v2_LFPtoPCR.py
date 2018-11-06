# coding: utf-8
"""

GLOFRIM 2.0

This script executes the 2way coupling between large-scale hydrology (PCR-GLOBWB) and hydrodynamics (Delft3D FM or LISFLOODFP).
It requires two files with information about the coupling settings and the paths to the models in each specific environmentself.

As such, the script should be executed like:
    python couplingFramework_v2_LFPtoPCR.py settingsFile.set paths.env True/False discharge/runoff

If the last switch is True, then a preparatory run is initiated to perform the grid coupling onlyself.
If it is False, then the already prepared grid coupling py-objects are used to perform the remaining computations and actual exchange.

Further information can be found on GitHub or by contacting the author.

Author: Jannis Hoch, M.Sc.
Mail: j.m.hoch@uu.nl
Institute: Utrecht University / Deltares

(C) 2018

"""

# -------------------------------------------------------------------------------------------------
# LOAD REQUIRED LIBRARIES
# -------------------------------------------------------------------------------------------------

import netCDF4
from distutils.util import strtobool
import pylab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pcraster as pcr
import sys
import os
import platform
import numpy as np
import pyproj as pyproj
import datetime
import pdb
import time
import pickle
import bmi.wrapper
import model as pcrglobwb_bmi_v203
from model import pcrglobwb_bmi
from model import disclaimer
from coupling_PCR_FM_2way import coupling_functions
from coupling_PCR_FM_2way import model_functions
from coupling_PCR_FM_2way import utils
from coupling_PCR_FM_2way import configuration

# -------------------------------------------------------------------------------------------------
# IMPORT MODEL SETTINGS FROM INI-FILE/SET-FILE
# -------------------------------------------------------------------------------------------------

config = configuration.Configuration()
config.parse_configuration_file(sys.argv[1])

# -------------------------------------------------------------------------------------------------
# DEFINE WHETHER TO PERFORM PREPARATORY OR ACTUAL RUN
# -------------------------------------------------------------------------------------------------

prepRun = strtobool(sys.argv[3])

if prepRun == True:
    print '\n###############################'
    print '### PREPARATORY RUN ACTIVE! ###'
    print '###############################\n'
elif prepRun == False:
    print '\n########################'
    print '### ACTUAL RUN ACTIVE! ###'
    print '########################\n'
else:
    sys.exit('\nplease indicate properly whether to perform preparatory or actual run!\n')

# -------------------------------------------------------------------------------------------------
# SPECIFY MODEL SETTINGS
# -------------------------------------------------------------------------------------------------

model_type = config.model_type['model_type']

use_2way = strtobool(config.general_settings['2way coupling'])

latlon = strtobool(config.general_settings['latlon'])
if latlon == False:
    inProj  = pyproj.Proj(init=config.model_settings['model_projection'])

use_Fluxes = strtobool(config.general_settings['use_Fluxes'])

use_RFS = strtobool(config.general_settings['use_RFS'])
if use_2way == True:
    use_RFS = True

verbose = strtobool(config.general_settings['verbose'])

couple_channelStorage = True

try:
    sys.argv[4]
except NameError:
    forcing = 'discharge'
else:
    forcing = str(sys.argv[4])


# -------------------------------------------------------------------------------------------------
# SPECIFY NUMERICAL SETTINGS
# -------------------------------------------------------------------------------------------------

nr_model_timesteps                    = int(config.numerical_settings['number_of_timesteps'])
update_step                           = int(config.numerical_settings['update_step'])

secPerDay                             = 86400.
end_time 							  = nr_model_timesteps * secPerDay
fraction_timestep 					  = secPerDay / update_step

threshold_inundated_depth             = float(config.numerical_settings['threshold_inundated_depth'])

# other
missing_value_landmask                = 255
missing_value_pcr                     = -999

# -------------------------------------------------------------------------------------------------
# SET PATHS TO MODELS
# -------------------------------------------------------------------------------------------------

# hydrodynamics
model_dir       	= config.model_settings['model_dir']
model_file      	= config.model_settings['model_file']
model_proj			= config.model_settings['model_projection']

# hydrology
config_pcr       	=  config.PCR_settings['config_pcr']
configPCR           = configuration.Configuration()
configPCR.parse_configuration_file(config_pcr )
inputDIR 			= configPCR.globalOptions['inputDir']
clone_pcr 			= os.path.join(inputDIR, configPCR.globalOptions['cloneMap'])
landmask_pcr 		= os.path.join(inputDIR, configPCR.globalOptions['landmask'])
LDD                 = os.path.join(inputDIR, configPCR.routingOptions['lddMap'])

# -------------------------------------------------------------------------------------------------
# SET PATHS TO .SO / .DLL FILES
# -------------------------------------------------------------------------------------------------

configPaths = configuration.Configuration()
configPaths.parse_configuration_file(sys.argv[2])

# these may be changed according to personal file and folder structure
if model_type == 'DFM':
    model_path = configPaths.paths['DFM_path']

elif model_type == 'LFP':
    model_path = configPaths.paths['LFP_path']

else:
    sys.exit('\nno adequate model defined in configuration file - define either DFM or LFP!\n')

# -------------------------------------------------------------------------------------------------
# INITIALIZE AND SPIN-UP PCR-GLOBWB
# -------------------------------------------------------------------------------------------------

# initiate logging and define folder for verbose-output
verbose_folder = model_functions.write2log(model_dir, model_file, latlon, use_2way, use_Fluxes, use_RFS, prepRun, use_floodplain_infiltration_factor=False, adjust_initial_groundwater=False, verbose=False, moment='start')
startTime = datetime.datetime.now()
print 'Model Start-Time: ', str(startTime) + os.linesep

# initiate PCR-GLOBWB
hydrologicModel = pcrglobwb_bmi_v203.pcrglobwb_bmi.pcrglobwbBMI()
hydrologicModel.initialize(config_pcr)
print '\n>>> PCR Initialized <<<\n'

# spin-up PCR-GLOBWB (only for actual run)
if prepRun == False:
    hydrologicModel.spinup()

# -------------------------------------------------------------------------------------------------
# INITIALIZING HYDRODYNAMIC MODEL
# -------------------------------------------------------------------------------------------------

# initiate hydraulic model
hydrodynamicModel = bmi.wrapper.BMIWrapper(engine = model_path, configfile = (os.path.join(model_dir, model_file)))
hydrodynamicModel.initialize()
print '\n>>> ',model_type,' Initialized <<<\n'

# -------------------------------------------------------------------------------------------------
# EXCTRACTING RELEVANT DATA FROM MODELS
# -------------------------------------------------------------------------------------------------

#- from Delft3D FM
if model_type == 'DFM':

    #- retrieving data from Delft3D FM
    bottom_lvl_1D, bottom_lvl_2D, cellAreaSpherical_1D, cellAreaSpherical_2D,         hydrodynamic_coords_1D, hydrodynamic_coords_2D, separator_1D  = model_functions.extractModelData_DFM(hydrodynamicModel,
                                                                                                             use_RFS,
                                                                                                             use_2way)

    print '\n>>> DFM data retrieved <<<\n'

#- from LISFLOODFP
elif model_type == 'LFP':

    #- retrieving data from LISFLOOD-FP
    DEM, cellArea, SGCQin, SGCwidth, SGCz, \
        list_x_coords, list_y_coords, \
        list_x_coords_2way, list_y_coords_2way, \
        i_1d, j_1d, \
        i_2d, j_2d, \
        dx, dy  = model_functions.extractModelData_LFP(hydrodynamicModel, use_RFS, use_2way)

    #- construct coords of vertices from midpoint coords and grid resolution
    hydrodynamic_coords_1D = coupling_functions.getVerticesFromMidPoints(list_x_coords, list_y_coords, dx, dy, verbose)
    if use_2way == True:
        hydrodynamic_coords_2D = coupling_functions.getVerticesFromMidPoints(list_x_coords_2way, list_y_coords_2way, dx, dy, verbose)

    #- construcut arrays with indices pointing to 1way/2way coupled LFP cells in a 2D grid
    coupledFPindices_1way = zip(i_1d, j_1d)
    if use_2way == True:
        coupledFPindices_2way = zip(i_2d, j_2d)

    #- align LFP with DFM
    bottom_lvl_1D = SGCz[i_1d,j_1d]
    if use_2way == True:
        bottom_lvl_2D = SGCz[i_2d,j_2d]
    cellAreaSpherical_1D = cellArea[i_1d,j_1d]
    if use_2way == True:
        cellAreaSpherical_2D = cellArea[i_2d,j_2d]

    separator_1D = 0.

    print '\n>>> LFP data retrieved <<<\n'

#- from PCR-GLOBWB
cellarea_data_pcr, landmask_data_pcr, clone_data_pcr = model_functions.extractModelData_PCR(hydrologicModel,
                                                                                            landmask_pcr,
                                                                                            clone_pcr)
PCRcoords = coupling_functions.getPCRcoords(landmask_data_pcr)
print '\n>>> PCR data retrieved <<<\n'

# -------------------------------------------------------------------------------------------------
# PLOT LFP DEM, CHANNEL DEPTH, CHANNEL WIDTH
# -------------------------------------------------------------------------------------------------

if (model_type == 'LFP') and (verbose==True):
    plt.figure(figsize=(16,8))
    plt.subplot(131)
    plt.imshow(DEM, vmax=200)
    plt.colorbar()
    plt.subplot(132)
    plt.imshow(SGCz, vmax=200)
    plt.colorbar()
    plt.subplot(133)
    plt.imshow(SGCwidth, vmax=500)
    plt.colorbar()
    plt.savefig(os.path.join(verbose_folder, 'LFP_properties.png'), dpi=300)

# -------------------------------------------------------------------------------------------------
# ACTIVATING COUPLING FOR RELEVANT PCR SECTIONS
# -------------------------------------------------------------------------------------------------

hydrologicModel.set_var(('grassland','ActivateCoupling'), 'True') #2way
hydrologicModel.set_var(('forest','ActivateCoupling'), 'True') #2way
hydrologicModel.set_var(('routing','ActivateCoupling'), 'True')
hydrologicModel.set_var(('WaterBodies', 'ActivateCoupling'), 'True')

# -------------------------------------------------------------------------------------------------
# COUPLING THE GRIDS
# -------------------------------------------------------------------------------------------------

if prepRun == True:

    print '\ncoupling hydrology with hydrodynamics for 1D...'

    # linking PCR-cells with 1D hydrodynamic cells (since use_RFS=True)
    couple_HDYN_2_HLOG, couple_HLOG_2_HDYN, coupled_HLOG_indices = coupling_functions.assignPCR2cells(landmask_pcr,
                                                                                                      hydrodynamic_coords_1D,
                                                                                                      verbose)

    # print len(hydrodynamic_coords_1D)
    # print len(couple_HDYN_2_HLOG)

    print '...pickling data for 1D to files...'

    f1 = open(os.path.join(verbose_folder,'couple_HDYN_2_HLOG.pkl'), 'wb')
    f2 = open(os.path.join(verbose_folder,'couple_HLOG_2_HDYN.pkl'), 'wb')
    f3 = open(os.path.join(verbose_folder,'coupled_HLOG_indices.pkl'), 'wb')

    pickle.dump(couple_HDYN_2_HLOG, f1, -1)
    pickle.dump(couple_HLOG_2_HDYN, f2, -1)
    pickle.dump(coupled_HLOG_indices, f3, -1)

    f1.close()
    f2.close()
    f3.close()

    print '...succesfully pickled data for 1D!\n'

    if use_2way == True:

        print 'coupling hydrology with hydrodynamics for 2D...'

        # linking PCR-cells with 2D hydrodynamic cells
        couple_HDYN_2_HLOG_2way, couple_HLOG_2_HDYN_2way, coupled_HLOG_indices_2way = coupling_functions.assignPCR2cells(landmask_pcr,
                                                                                                                         hydrodynamic_coords_2D,
                                                                                                                         verbose)

        # print len(hydrodynamic_coords_2D)
        # print len(couple_HDYN_2_HLOG_2way)

        print '...pickling data for 2D to files...'

        f4 = open(os.path.join(verbose_folder,'couple_HDYN_2_HLOG_2way.pkl'), 'wb')
        f5 = open(os.path.join(verbose_folder,'couple_HLOG_2_HDYN_2way.pkl'), 'wb')
        f6 = open(os.path.join(verbose_folder,'coupled_HLOG_indices_2way.pkl'), 'wb')

        pickle.dump(couple_HDYN_2_HLOG_2way, f4, -1)
        pickle.dump(couple_HLOG_2_HDYN_2way, f5, -1)
        pickle.dump(coupled_HLOG_indices_2way, f6, -1)

        f4.close()
        f5.close()
        f6.close()

        print '...succesfully pickled data for 2D!\n'

    sys.exit('\nAll preparatory files have successfully been pickled\n')

elif prepRun == False:

    print '\nloading pickled files...'

    f1 = open(os.path.join(verbose_folder,'couple_HDYN_2_HLOG.pkl'), 'rb')
    f2 = open(os.path.join(verbose_folder,'couple_HLOG_2_HDYN.pkl'), 'rb')
    f3 = open(os.path.join(verbose_folder,'coupled_HLOG_indices.pkl'), 'rb')

    couple_HDYN_2_HLOG = pickle.load(f1)
    couple_HLOG_2_HDYN = pickle.load(f2)
    coupled_HLOG_indices = pickle.load(f3)

    if use_2way == True:

        f4 = open(os.path.join(verbose_folder,'couple_HDYN_2_HLOG_2way.pkl'), 'rb')
        f5 = open(os.path.join(verbose_folder,'couple_HLOG_2_HDYN_2way.pkl'), 'rb')
        f6 = open(os.path.join(verbose_folder,'coupled_HLOG_indices_2way.pkl'), 'rb')

        couple_HDYN_2_HLOG_2way = pickle.load(f4)
        couple_HLOG_2_HDYN_2way = pickle.load(f5)
        coupled_HLOG_indices_2way = pickle.load(f6)

    print '...pickled files succesfully loaded!\n'

# -------------------------------------------------------------------------------------------------
# EXTRACT CURRENT WATER DEPTH AND ASSIGNT TO 1D OR 2D CELLS
# -------------------------------------------------------------------------------------------------

if model_type == 'DFM':
    current_water_depth = hydrodynamicModel.get_var('s1') - hydrodynamicModel.get_var('bl')
    current_water_depth_1D = current_water_depth[separator_1D:]
    if use_2way == True:
        current_water_depth_2D = current_water_depth[:separator_1D]
elif model_type == 'LFP':
    current_water_depth = np.copy(hydrodynamicModel.get_var('H'))
    current_water_depth_1D = current_water_depth[i_1d,j_1d]
    if use_2way == True:
        current_water_depth_2D = current_water_depth[i_2d,j_2d]

# -------------------------------------------------------------------------------------------------
# DETERMINE INUNDATION AREA OF 1D CELLS
# -------------------------------------------------------------------------------------------------

inundatedArea_HDYN1D_2_HLOG,\
    inundatedArea_HDYN1D_2_HLOG_BMI = model_functions.determine_InundationArea1D_Hydrodynamics(couple_HLOG_2_HDYN,
                                                                                               coupled_HLOG_indices,
                                                                                               current_water_depth_1D,
                                                                                               threshold_inundated_depth,
                                                                                               cellAreaSpherical_1D,
                                                                                               cellarea_data_pcr,
                                                                                               landmask_pcr)

# -------------------------------------------------------------------------------------------------
# DETERMINE INUNDATION AREA AND INUNDATED AREA FRACTION OF 2D CELLS
# -------------------------------------------------------------------------------------------------

if use_2way == True:
    inundatedArea_HDYN2D_2_HLOG,\
        inundatedArea_HDYN2D_2_HLOG_BMI,\
            inundatedFraction_HDYN2D_2_HLOG,\
                inundatedFraction_HDYN2D_2_HLOG_BMI = model_functions.determine_InundationArea2D_Hydrodynamics(couple_HLOG_2_HDYN_2way,
                                                                                                               coupled_HLOG_indices_2way,
                                                                                                               current_water_depth_2D,
                                                                                                               threshold_inundated_depth,
                                                                                                               cellAreaSpherical_2D,
                                                                                                               cellarea_data_pcr,
                                                                                                               landmask_pcr)

# -------------------------------------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------------------------------------

if (use_2way == True) and (verbose == True):
    plt.figure(figsize=(16,8))
    plt.subplot(131)
    plt.title('inundatedArea (1D) [m2]')
    plt.imshow(inundatedArea_HDYN1D_2_HLOG_BMI, vmax=np.max(inundatedArea_HDYN1D_2_HLOG_BMI))
    plt.colorbar(orientation='horizontal')
    plt.subplot(132)
    plt.title('inundatedArea (2D) [m2]')
    plt.imshow(inundatedArea_HDYN2D_2_HLOG_BMI, vmax=np.max(inundatedArea_HDYN1D_2_HLOG_BMI))
    plt.colorbar(orientation='horizontal')
    plt.subplot(133)
    plt.title('inundatedFraction (2D) [m2/m2]')
    plt.imshow(inundatedFraction_HDYN2D_2_HLOG_BMI, vmin=0, vmax=1)
    plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(verbose_folder,'inundatedAreaAndFraction_HDYN_2_HLOG_BMI.png'), dpi=300)

# -------------------------------------------------------------------------------------------------
# EXTRACT CURRENT WATER VOLUME AND ASSIGNT TO 1D OR 2D CELLS
# -------------------------------------------------------------------------------------------------

if model_type == 'DFM':
    current_water_volume = np.copy(hydrodynamicModel.get_var('vol1')) # all hydrodyancmi cells
    current_water_volume_1D = current_water_volume[separator_1D:] # 1D channel cells
    if use_2way == True:
        current_water_volume_2D = current_water_volume[:separator_1D] # 2D floodplain cells
elif model_type == 'LFP':
    current_water_volume = cellArea * hydrodynamicModel.get_var('H')
    current_water_volume_1D = current_water_volume[i_1d, j_1d]
    if use_2way == True:
        current_water_volume_2D = current_water_volume[i_2d, j_2d]

# -------------------------------------------------------------------------------------------------
# DETERMINE INUNDATION VOLUME AND DEPTH OF 1D CELLS
# -------------------------------------------------------------------------------------------------

waterVolume_HDYN1D_2_HLOG,\
    waterVolume_HDYN1D_2_HLOG_BMI,\
        waterDepth_HDYN1D_2_HLOG_BMI = model_functions.determine_inundationVolume_HDYN_1D(hydrologicModel,
                                                                                          inundatedArea_HDYN1D_2_HLOG,
                                                                                          current_water_volume_1D,
                                                                                          couple_HLOG_2_HDYN,
                                                                                          coupled_HLOG_indices,
                                                                                          landmask_pcr)

# -------------------------------------------------------------------------------------------------
# DETERMINE INUNDATION VOLUME AND DEPTH OF 2D CELLS
# -------------------------------------------------------------------------------------------------

if use_2way == True:
    waterVolume_HDYN2D_2_HLOG_BMI,\
        waterDepth_HDYN2D_2_HLOG_BMI = model_functions.determine_inundationVolume_HDYN_2D(current_water_volume_2D,
                                                                                          inundatedArea_HDYN2D_2_HLOG,
                                                                                          waterVolume_HDYN1D_2_HLOG_BMI,
                                                                                          couple_HLOG_2_HDYN_2way,
                                                                                          coupled_HLOG_indices_2way,                                                                                          landmask_pcr)

if use_2way == False:
    waterVolume_HDYN2D_2_HLOG_BMI = np.zeros_like(waterVolume_HDYN1D_2_HLOG_BMI)

# -------------------------------------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------------------------------------

if (use_2way == True) and (verbose == True):
    plt.figure(figsize=(16,8))
    plt.subplot(141)
    plt.title('waterVolume (1D) [m3]')
    plt.imshow(waterVolume_HDYN1D_2_HLOG_BMI)
    plt.colorbar(orientation='horizontal')
    plt.subplot(142)
    plt.title('waterDepth (1D) [m3]')
    plt.imshow(waterDepth_HDYN1D_2_HLOG_BMI)
    plt.colorbar(orientation='horizontal')
    plt.subplot(143)
    plt.title('waterVolume (2D) [m3]')
    plt.imshow(waterVolume_HDYN2D_2_HLOG_BMI)
    plt.colorbar(orientation='horizontal')
    plt.subplot(144)
    plt.title('waterDepth (2D) [m]')
    plt.imshow(waterDepth_HDYN2D_2_HLOG_BMI)
    plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(verbose_folder,'inundatedVolumeAndDepth_HDYN_2_HLOG_BMI.png'), dpi=300)

# -------------------------------------------------------------------------------------------------
# UPDATE STORAGES IN PCR: SETTING WATERBODYSTORAGE TO 0 AND CHANNELSTORAGE TO 1D VOLUME FROM HYD
# -------------------------------------------------------------------------------------------------

if use_2way == True:
    model_functions.updateStorage(hydrologicModel,
                                  missing_value_pcr,
                                  waterVolume_HDYN1D_2_HLOG_BMI)

# -------------------------------------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------------------------------------

if (use_2way == True) and (verbose == True):
    plt.figure(figsize=(16,8))
    plt.subplot(131)
    plt.title('initial channelStorage [m3]')
    plt.imshow(hydrologicModel.get_var('channelStorage'), vmax=np.max(waterVolume_HDYN1D_2_HLOG_BMI))
    plt.colorbar(orientation='horizontal')
    plt.subplot(132)
    plt.title('waterVolume (1D) [m3]')
    plt.imshow(waterVolume_HDYN1D_2_HLOG_BMI, vmax=np.max(waterVolume_HDYN1D_2_HLOG_BMI))
    plt.colorbar(orientation='horizontal')
    plt.subplot(133)
    plt.title('updated channelStorage [m3]')
    plt.imshow(hydrologicModel.get_var('channelStorage'), vmax=np.max(waterVolume_HDYN1D_2_HLOG_BMI))
    plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(verbose_folder,'iniAndUpdatedChannelStoragePCR.png'), dpi=300)

# -------------------------------------------------------------------------------------------------
# UPDATING PCR WITH FLOODPLAIN INUNDATION FRACTION AND VOLUME ON SURFACE
# -------------------------------------------------------------------------------------------------

if use_2way == True:
    model_functions.updateHydrologicVariables(hydrologicModel,
                                              inundatedFraction_HDYN2D_2_HLOG_BMI,
                                              waterVolume_HDYN2D_2_HLOG_BMI)

# -------------------------------------------------------------------------------------------------
# UPDATING PCR TO NEXT TIME STEP
# -------------------------------------------------------------------------------------------------

hydrologicModel.update(1)

# -------------------------------------------------------------------------------------------------
# COMPUTING DELTA VOLUME TO BE COUPLED FROM PCR TO HYDRODYNAMIC MODEL
# -------------------------------------------------------------------------------------------------

delta_volume_PCR, delta_volume_PCR_1way = model_functions.calculateDeltaVolumes(hydrologicModel,
                                                                                missing_value_pcr,
                                                                                secPerDay,
                                                                                coupled_HLOG_indices,
                                                                                cellarea_data_pcr,
                                                                                waterVolume_HDYN2D_2_HLOG_BMI,
                                                                                forcing)

# -------------------------------------------------------------------------------------------------
# PLOT
# -------------------------------------------------------------------------------------------------

if verbose == True:
    r = np.copy(hydrologicModel.get_var('landSurfaceRunoff'))
    r = r * cellarea_data_pcr
    r[r<=0.] = -999
    Q = np.copy(hydrologicModel.get_var('discharge'))
    Q = Q*86400.
    Q[Q<=0.] = -999

    plt.figure(figsize=(16,8))
    plt.subplot(131)
    plt.title('r from PCR [m3]')
    plt.imshow(r)
    plt.colorbar(orientation='horizontal')
    plt.subplot(132)
    plt.title('Q from PCR [m3]')
    plt.imshow(Q)
    plt.colorbar(orientation='horizontal')
    plt.subplot(133)
    plt.title('delta volume PCR [m3/d]')
    plt.imshow(delta_volume_PCR)
    plt.colorbar(orientation='horizontal')
    plt.savefig(os.path.join(verbose_folder,'coupledPCRvolumes.png'), dpi=300)

# -------------------------------------------------------------------------------------------------
# ASSIGNING DELTA VOLUME TO HYDRODYNAMIC CELLS
# -------------------------------------------------------------------------------------------------

delta_water_DFM_1way, verbose_volume_DFM_1way = model_functions.calculateDeltaWater(hydrodynamicModel,
                                                                                    couple_HDYN_2_HLOG,
                                                                                    couple_HLOG_2_HDYN,
                                                                                    coupled_HLOG_indices,
                                                                                    delta_volume_PCR_1way,
                                                                                    cellAreaSpherical_1D,
                                                                                    fraction_timestep,
                                                                                    model_type,
                                                                                    use_Fluxes)

# -------------------------------------------------------------------------------------------------
# UPDATING HYDRODYNAMIC MODEL
# -------------------------------------------------------------------------------------------------

if model_type == 'LFP':

    delta_water_DFM_1way = model_functions.fillLFPgrid(hydrodynamicModel,
                                                       coupledFPindices_1way,
                                                       delta_water_DFM_1way,
                                                       DEM,
                                                       verbose_folder,
                                                       verbose=False)

    model_functions.updateModel(hydrodynamicModel,
                                delta_water_DFM_1way,
                                update_step,
                                separator_1D,
                                use_Fluxes,
                                use_RFS,
                                model_type,
                                verbose=False)

while hydrodynamicModel.get_current_time() < (hydrologicModel.get_time_step() * secPerDay):

    if model_type == 'LFP':
        hydrodynamicModel.update()

    elif model_type == 'DFM':
        model_functions.updateModel(hydrodynamicModel,
                                delta_water_DFM_1way,
                                update_step,
                                separator_1D,
                                use_Fluxes,
                                use_RFS,
                                model_type,
                                verbose)

# -------------------------------------------------------------------------------------------------
# REPEATING STEPS UNTIL END OF MODELLING PERIOD
# -------------------------------------------------------------------------------------------------

while hydrodynamicModel.get_time_step() < end_time:

    if model_type == 'DFM':
        current_water_depth = hydrodynamicModel.get_var('s1') - hydrodynamicModel.get_var('bl')
        current_water_depth_1D = current_water_depth[separator_1D:]
        if use_2way == True:
            current_water_depth_2D = current_water_depth[:separator_1D]
    elif model_type == 'LFP':
        current_water_depth = np.copy(hydrodynamicModel.get_var('H'))
        current_water_depth_1D = current_water_depth[i_1d,j_1d]
        if use_2way == True:
            current_water_depth_2D = current_water_depth[i_2d,j_2d]

    inundatedArea_HDYN1D_2_HLOG,\
        inundatedArea_HDYN1D_2_HLOG_BMI = model_functions.determine_InundationArea1D_Hydrodynamics(couple_HLOG_2_HDYN,
                                                                                                   coupled_HLOG_indices,
                                                                                                   current_water_depth_1D,
                                                                                                   threshold_inundated_depth,
                                                                                                   cellAreaSpherical_1D,
                                                                                                   cellarea_data_pcr,
                                                                                                   landmask_pcr)

    if use_2way == True:
        inundatedArea_HDYN2D_2_HLOG,\
            inundatedArea_HDYN2D_2_HLOG_BMI,\
                inundatedFraction_HDYN2D_2_HLOG,\
                    inundatedFraction_HDYN2D_2_HLOG_BMI = model_functions.determine_InundationArea2D_Hydrodynamics(couple_HLOG_2_HDYN_2way,
                                                                                                                   coupled_HLOG_indices_2way,
                                                                                                                   current_water_depth_2D,
                                                                                                                   threshold_inundated_depth,
                                                                                                                   cellAreaSpherical_2D,
                                                                                                                   cellarea_data_pcr,
                                                                                                                   landmask_pcr)

    if model_type == 'DFM':
        current_water_volume = np.copy(hydrodynamicModel.get_var('vol1')) # all hydrodyancmi cells
        current_water_volume_1D = current_water_volume[separator_1D:] # 1D channel cells
        if use_2way == True:
            current_water_volume_2D = current_water_volume[:separator_1D] # 2D floodplain cells
    elif model_type == 'LFP':
        current_water_volume = cellArea * hydrodynamicModel.get_var('H')
        current_water_volume_1D = current_water_volume[i_1d, j_1d]
        if use_2way == True:
            current_water_volume_2D = current_water_volume[i_2d, j_2d]

    waterVolume_HDYN1D_2_HLOG,\
        waterVolume_HDYN1D_2_HLOG_BMI,\
            waterDepth_HDYN1D_2_HLOG_BMI = model_functions.determine_inundationVolume_HDYN_1D(hydrologicModel,
                                                                                              inundatedArea_HDYN1D_2_HLOG,
                                                                                              current_water_volume_1D,
                                                                                              couple_HLOG_2_HDYN,
                                                                                              coupled_HLOG_indices,
                                                                                              landmask_pcr)

    if use_2way == True:
        waterVolume_HDYN2D_2_HLOG_BMI,\
            waterDepth_HDYN2D_2_HLOG_BMI = model_functions.determine_inundationVolume_HDYN_2D(current_water_volume_2D,
                                                                                              inundatedArea_HDYN2D_2_HLOG,
                                                                                              waterVolume_HDYN1D_2_HLOG_BMI,
                                                                                              couple_HLOG_2_HDYN_2way,
                                                                                              coupled_HLOG_indices_2way,
                                                                                              landmask_pcr)

    if use_2way == False:
        waterVolume_HDYN2D_2_HLOG_BMI = np.zeros_like(waterVolume_HDYN1D_2_HLOG_BMI)

    if use_2way == True:
        model_functions.updateStorage(hydrologicModel,
                                      missing_value_pcr,
                                      waterVolume_HDYN1D_2_HLOG_BMI)


    if use_2way == True:
        model_functions.updateHydrologicVariables(hydrologicModel,
                                                  inundatedFraction_HDYN2D_2_HLOG_BMI,
                                                  waterVolume_HDYN2D_2_HLOG_BMI)

    hydrologicModel.update(1)

    delta_volume_PCR, delta_volume_PCR_1way = model_functions.calculateDeltaVolumes(hydrologicModel,
                                                                            missing_value_pcr,
                                                                            secPerDay,
                                                                            coupled_HLOG_indices,
                                                                            cellarea_data_pcr,
                                                                            waterVolume_HDYN2D_2_HLOG_BMI,
                                                                            forcing)

    delta_water_DFM_1way, verbose_volume_DFM_1way = model_functions.calculateDeltaWater(hydrodynamicModel,
                                                                                couple_HDYN_2_HLOG,
                                                                                couple_HLOG_2_HDYN,
                                                                                coupled_HLOG_indices,
                                                                                delta_volume_PCR_1way,
                                                                                cellAreaSpherical_1D,
                                                                                fraction_timestep,
                                                                                model_type,
                                                                                use_Fluxes)
    if model_type == 'LFP':

        delta_water_DFM_1way = model_functions.fillLFPgrid(hydrodynamicModel,
                                                           coupledFPindices_1way,
                                                           delta_water_DFM_1way,
                                                           DEM,
                                                           verbose_folder,
                                                           verbose=False)

        model_functions.updateModel(hydrodynamicModel,
                                    delta_water_DFM_1way,
                                    update_step,
                                    separator_1D,
                                    use_Fluxes,
                                    use_RFS,
                                    model_type,
                                    verbose=False)

    while hydrodynamicModel.get_current_time() < (hydrologicModel.get_time_step() * secPerDay):

        if model_type == 'LFP':
            hydrodynamicModel.update()

        elif model_type == 'DFM':
            model_functions.updateModel(hydrodynamicModel,
                                    delta_water_DFM_1way,
                                    update_step,
                                    separator_1D,
                                    use_Fluxes,
                                    use_RFS,
                                    model_type,
                                    verbose)

# -------------------------------------------------------------------------------------------------
# FINALIZING MODELS
# -------------------------------------------------------------------------------------------------

hydrodynamicModel.finalize()
hydrologicModel.finalize()

endTime = datetime.datetime.now()
print 'Model Start-Time: ', str(startTime)
print 'Model End-Time: ', str(endTime) + os.linesep
