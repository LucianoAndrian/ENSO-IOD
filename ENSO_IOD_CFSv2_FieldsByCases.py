"""
Anomalía de PP (y T cuando ande) en SA según la magnitud de los índices
(Similar a *2D_bins_DMI_N34_3.0.py, pero en lugar de regiones para todo SA)
"""
########################################################################################################################
import xarray as xr
import numpy as np
import os
from matplotlib import colors
from ENSO_IOD_Funciones import SelectNMMEFiles, fix_calendar, ComputeFieldsByCases
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/' # índices por estaciones
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/' # campos de las variables PP  ( y T cuando ande)
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/CuadByCases/'
out_data_dir = '/pikachu/datos/luciano.andrian/cases_fields/'

save = True
dpi = 500
########################################################################################################################
cases = ['dmi_puros_pos', 'dmi_puros_neg',
        'n34_puros_pos', 'n34_puros_neg',
        'sim_pos', 'sim_neg',
        'dmi_neg_n34_pos', 'dmi_pos_n34_neg',
        'neutros']

bin_limits = [[-4.5,-1], [-1, -0.5], #1
              [-0.5, 0.5], #2
              [0.5, 1], [1, 4.5]] #4

bins_by_cases_dmi = [[3, 4], [0, 1],
                     [2], [2],
                     [3, 4], [0, 1],
                     [0, 1],[3, 4],
                     [2]]

bins_by_cases_n34 = [[2], [2],
                     [3, 4], [0,1],
                     [3, 4], [0, 1],
                     [3, 4], [0, 1],
                     [2]]

cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cbar_snr = colors.ListedColormap(['#002A3D','#074D4F', '#1E6D5A' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#DCBC75', '#995D13','#6A3D07','#543005','#3F2404'][::-1])
cbar_snr.set_under('#3F2404')
cbar_snr.set_over('#002A3D')
cbar_snr.set_bad(color='white')

cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')

cbar_snr_t = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#FEB77E', '#FB8761','#CA3E72','#782281','#251255'])
cbar_snr_t.set_over('#251255')
cbar_snr_t.set_under('#070B4F')
cbar_snr_t.set_bad(color='white')

cbar_hgt = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar_hgt.set_over('#641B00')
cbar_hgt.set_under('#012A52')
cbar_hgt.set_bad(color='white')

########################################################################################################################
# Prec #################################################################################################################
# HINDCAST para climatología ------------------------------------------------------------------------------------------#
try:
    data = xr.open_dataset(out_data_dir + 'prec_data_full.nc')
except:
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable='prec',
                            dir='/pikachu/datos/osman/nmme/monthly/hindcast/', All=True)
    files = sorted(files, key=lambda x: x.split()[0])

    # abriendo todos los archivos
    data = xr.open_mfdataset(files, decode_times=False)  # xr no entiende la codificacion de Leads, r y las fechas
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
    data = data.sel(r=slice(1, 24))
    data = data.rolling(time=3, center=True).mean()
    #data.to_netcdf(out_data_dir + 'prec_data_full.nc')
    data = data.compute()
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# Composite
levels = np.linspace(-30, 30, 13)
levels_clim = np.linspace(0,150,11)

ComputeFieldsByCases(v='prec', v_name='prec', fix_factor=30, snr=False,
                     data=data,
                     cases=cases, bin_limits=bin_limits,
                     bins_by_cases_dmi=bins_by_cases_dmi,
                     bins_by_cases_n34=bins_by_cases_n34,
                     cases_dir=cases_dir, dates_dir=dates_dir,
                     levels_main=levels, cbar_main=cbar_pp,
                     levels_clim=levels_clim, cbar_clim='terrain_r',
                     title_var='Prec', name_fig='prec', dpi=dpi, save=save,
                     out_dir=out_dir)

# Signal-to-Noise ratio
levels = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]

ComputeFieldsByCases(v='prec', v_name='prec', fix_factor=30, snr=True,
                     data=data,
                     cases=cases, bin_limits=bin_limits,
                     bins_by_cases_dmi=bins_by_cases_dmi,
                     bins_by_cases_n34=bins_by_cases_n34,
                     cases_dir=cases_dir, dates_dir=dates_dir,
                     levels_main=levels, cbar_main=cbar_snr,
                     levels_clim=levels_clim, cbar_clim='terrain_r',
                     title_var='Prec', name_fig='prec', dpi=dpi, save=save,
                     out_dir=out_dir)

########################################################################################################################
# Temp #################################################################################################################
# HINDCAST para climatología ------------------------------------------------------------------------------------------#
try:
    data = xr.open_dataset(out_data_dir + 'temp_data_full.nc')
except:
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable='tref',
                            dir='/pikachu/datos/osman/nmme/monthly/hindcast/', All=True)
    files = sorted(files, key=lambda x: x.split()[0])

    # abriendo todos los archivos
    data = xr.open_mfdataset(files, decode_times=False)  # xr no entiende la codificacion de Leads, r y las fechas
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
    data = data.sel(r=slice(1, 24))

    # media movil de 3 meses para separar en estaciones
    data = data.rolling(time=3, center=True).mean()
    #data = data.to_netcdf(out_data_dir + 'temp_data_full.nc')
    data = data.compute()

#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# Composite
levels =  np.linspace(-1.2,1.2,13)
levels_clim = np.linspace(0,25,11)

ComputeFieldsByCases(v='tref', v_name='tref', fix_factor=1, snr=False,
                     data=data-273,
                     cases=cases, bin_limits=bin_limits,
                     bins_by_cases_dmi=bins_by_cases_dmi,
                     bins_by_cases_n34=bins_by_cases_n34,
                     cases_dir=cases_dir, dates_dir=dates_dir,
                     levels_main=levels, cbar_main=cbar_t,
                     levels_clim=levels_clim, cbar_clim='Spectral_r',
                     title_var='Temp', name_fig='temp', dpi=dpi, save=save,
                     out_dir=out_dir)

# Signal-to-Noise ratio
levels = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]

ComputeFieldsByCases(v='tref', v_name='tref', fix_factor=1, snr=True,
                     data=data-273,
                     cases=cases, bin_limits=bin_limits,
                     bins_by_cases_dmi=bins_by_cases_dmi,
                     bins_by_cases_n34=bins_by_cases_n34,
                     cases_dir=cases_dir, dates_dir=dates_dir,
                     levels_main=levels, cbar_main=cbar_snr_t,
                     levels_clim=levels_clim, cbar_clim='Spectral_r',
                     title_var='Tref', name_fig='tref', dpi=dpi, save=save,
                     out_dir=out_dir)

########################################################################################################################
# HGT #################################################################################################################
# HINDCAST para climatología ------------------------------------------------------------------------------------------#

try:
    data = xr.open_dataset(out_data_dir + 'hgt_data_full.nc')
except:
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable='hgt',
                            dir='/pikachu/datos/luciano.andrian/hindcast/', All=True)
    files = sorted(files, key=lambda x: x.split()[0])

    # abriendo todos los archivos
    data = xr.open_mfdataset(files, decode_times=False)  # xr no entiende la codificacion de Leads, r y las fechas
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    #data = data.sel(lon=slice(30, 340), lat=slice(-80, 20), P=200)
    data = data.sel(lat=slice(-90, 20), P=200)
    data = data.drop('P')
    data = data.sel(r=slice(1, 24))

    # media movil de 3 meses para separar en estaciones
    data = data.rolling(time=3, center=True).mean()
    data = data.to_netcdf(out_data_dir + 'hgt_data_full_0_360.nc')
    #data = data.compute()

#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

levels = [-300,-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250,300]
levels_clim = np.linspace(10000,15000,11)

ComputeFieldsByCases(v='hgt', v_name='hgt', fix_factor=9.8, snr=False,
                     data=data*9.8,
                     cases=cases, bin_limits=bin_limits,
                     bins_by_cases_dmi=bins_by_cases_dmi,
                     bins_by_cases_n34=bins_by_cases_n34,
                     cases_dir=cases_dir, dates_dir=dates_dir,
                     levels_main=levels, cbar_main=cbar_hgt,
                     levels_clim=levels_clim, cbar_clim='Spectral',
                     title_var='HGT200_NC', name_fig='hgt200_neutro_clim', dpi=dpi,
                     x_lon=np.arange(30, 340, 25), x_lat=np.arange(-80, 20, 10),
                     figsize=[20,10], usemask=False, hcolorbar=True, save=save,
                     out_dir=out_dir)

# Signal-to-Noise ratio
levels = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]

ComputeFieldsByCases(v='hgt', v_name='hgt', fix_factor=9.8, snr=True,
                     data=data,
                     cases=cases, bin_limits=bin_limits,
                     bins_by_cases_dmi=bins_by_cases_dmi,
                     bins_by_cases_n34=bins_by_cases_n34,
                     cases_dir=cases_dir, dates_dir=dates_dir,
                     levels_main=levels, cbar_main=cbar_snr_t,
                     levels_clim=levels_clim, cbar_clim='YlGnBu',
                     title_var='HGT200_NC', name_fig='hgt200_neutro_clim', dpi=dpi,
                     x_lon=np.arange(30, 340, 25), x_lat=np.arange(-80, 20, 10),
                     figsize=[20, 10], usemask=False, hcolorbar=True, save=save,
                     out_dir=out_dir)
########################################################################################################################