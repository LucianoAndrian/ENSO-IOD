"""
Intengo de imagenes conceptuales de la circulacion ENSO-IOD
a partir de snr observado y del modelo cfvs2 y la regresion observada
"""
################################################################################
data_dir = '/pikachu/datos/luciano.andrian/esquemas/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/'
dir_results = 'esquemas'
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
from ENSO_IOD_Funciones import CreateDirectory, DirAndFile
import warnings
from shapely.errors import ShapelyDeprecationWarning
from scipy.signal import convolve2d
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
################################################################################
save = True
CreateDirectory(out_dir, dir_results)
#------------------------------------------------------------------------------#
if save:
    dpi = 600
else:
    dpi = 100
################################################################################
# Fase positiva -------------------------------------------------------------- #
ninio_obs = xr.open_dataset(data_dir +  'HGT200SON_N34_un_pos_mer_d_w_DMI_standard.nc')

dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
dmi = dmi.rename({'var_polyfit_coefficients':'var'})

sim = xr.open_dataset(data_dir + 'HGT200SON_DMI_sim_pos_mer_d_w_DMI_standard.nc')
dmi = xr.open_dataset(data_dir +  'HGT200SON_DMI_un_pos_mer_d_w_DMI_standard.nc')


Plot3(sim, [-2, -.5, 0, .5, 2], ninio_obs, [-.5, .5], dmi*0, [-.6, .6], '',
      DirAndFile(out_dir, dir_results, 'esquema_z200', '1'), dpi, save, 1)

# Fase negativa -------------------------------------------------------------- #
ninia_obs = xr.open_dataset(data_dir +  'HGT200SON_N34_un_neg_mer_d_w_DMI_standard.nc')

dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
dmi = dmi.rename({'var_polyfit_coefficients':'var'})

sim = xr.open_dataset(data_dir + 'HGT200SON_DMI_sim_neg_mer_d_w_DMI_standard.nc')
dmi = xr.open_dataset(data_dir +  'HGT200SON_DMI_un_neg_mer_d_w_DMI_standard.nc')


Plot3(sim, [-2, -.5, 0, .5, 2], ninia_obs, [-.5, .5], dmi*0, [-.4, .4],
      'El Niño (lineas solidas), IOD pos. (lineas punteadas) ' +
      '\n' + 'El Niño + IOD positivo (sombreado)',
      DirAndFile(out_dir, dir_results, 'esquema', '1.5'), dpi, save, 1)


# PP ------------------------------------------------------------------------- #
# Fase positiva -------------------------------------------------------------- #
ninio_mod = xr.open_dataset(data_dir + 'SNR_prec_CFSv2n34_puros_pos_SON.nc')
ninio_obs = xr.open_dataset(data_dir + 'ppgpccSON_N34_un_pos_mer_d_w_DMI_standard.nc')
ninio_obs = ninio_obs.interp(lon=ninio_mod.lon.values, lat=ninio_mod.lat.values)
ninio_obs=xr.where(np.isnan(ninio_obs), 0, ninio_obs)
# dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
# dmi = dmi.rename({'var_polyfit_coefficients':'var'})

sim_mod = xr.open_dataset(data_dir + 'SNR_prec_CFSv2sim_pos_SON.nc')
sim_obs = xr.open_dataset(data_dir + 'ppgpccSON_DMI_sim_pos_mer_d_w_DMI_standard.nc')
sim_obs = sim_obs.interp(lon=sim_mod.lon.values, lat=sim_mod.lat.values)
sim_obs=xr.where(np.isnan(sim_obs), 0, sim_obs)

ninio = (ninio_mod + 1.1*ninio_obs)/2
sim = (sim_mod*2 + sim_obs)/2

dmi = xr.open_dataset(data_dir + 'ppgpccSON_DMI_un_pos_mer_d_w_DMI_standard.nc')*0

Plot3b(sim, [-5,-.5, -.2, 0, .2, .5, 5],
      ninio, [-5,-.5, -.2, .2, .5, 5],
      dmi, [-2, 2], '',
      DirAndFile(out_dir, dir_results, 'esquema_pp', '2'),
      dpi, save, 1, True,
      smoothing=True, win1=(5,5), win2=(6,6))


####

