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
    dpi = 300
else:
    dpi = 100
################################################################################
def Plot3(data = None, levels = None, data2 = None, levels2 = None,
          data3 = None, levels3 = None, title = 'title', name = 'name',
          dpi = 70, save = False, step = 1, high=3.5):
    crs_latlon = ccrs.PlateCarree()

    fig = plt.figure(figsize=(7.08661, high), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_extent([0, 359, -80, 20], crs=crs_latlon)

    try:
        data_var = data['var']
        ax.contourf(data.lon[::step], data.lat[::step], data_var[::step, ::step],
                    linewidths=1, alpha=.4, levels=levels,
                    transform=crs_latlon,
                    colors=['dodgerblue', 'white', 'white', 'firebrick'])
    except:
        pass

    try:
        data2_var = data2['var']
        ax.contour(data2.lon[::step], data2.lat[::step], data2_var[::step, ::step],
                linewidths=1.5, alpha=.7, levels=levels2,
                transform=crs_latlon, colors=['Blue', 'Red'])
    except:
        pass

    try:
        data3_var = data3['var']
        ax.contour(data3.lon[::step], data3.lat[::step],
                   data3_var[::step, ::step],
                   linewidths=1.5, alpha=1, levels=levels3, linestyles='dashed',
                   transform=crs_latlon, colors=['#00B0CC', '#D00038'])
    except:
        pass

    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor='grey')
    # ax.set_xticks(np.arange(0, 360, 30), crs=crs_latlon)
    # ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color='grey', linestyle='-', alpha=1)
    # ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', zorder=20)
    # lon_formatter = LongitudeFormatter(zero_direction_label=True)
    # lat_formatter = LatitudeFormatter()
    # ax.xaxis.set_major_formatter(lon_formatter)
    # ax.yaxis.set_major_formatter(lat_formatter)
    # ax.tick_params(labelsize=7)
    # plt.title(title, fontsize=10)
    fig.text(0.5, 0.025, title, ha='center', fontsize=8)
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)

    fig.subplots_adjust(bottom=0.1, wspace=0, hspace=0.25, left=0.01,
                        right=0.99, top=0.975)
    #plt.tight_layout()
    if save:
        print('save: ' + name)
        plt.savefig(name, dpi=dpi)
        plt.close()

    else:
        plt.show()

# Fase positiva -------------------------------------------------------------- #
ninio_obs = xr.open_dataset(data_dir +  'HGT200SON_N34_un_pos_mer_d_w_DMI_standard.nc')

dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
dmi = dmi.rename({'var_polyfit_coefficients':'var'})

sim = xr.open_dataset(data_dir + 'HGT200SON_DMI_sim_pos_mer_d_w_DMI_standard.nc')
dmi = xr.open_dataset(data_dir +  'HGT200SON_DMI_un_pos_mer_d_w_DMI_standard.nc')

Plot3(sim, [-2, -.5, 0, .5, 2], ninio_obs, [-.5, .5], dmi*0, [-.6, .6],
      f"El Niño (solid lines), El Niño + positive IOD (shading)",
      DirAndFile(out_dir, dir_results, 'esquema_z200', '1', 'png'), dpi, save,
      1, 2.3)

Plot3(sim, [-2, -.5, 0, .5, 2], ninio_obs, [-.5, .5], dmi*0, [-.6, .6],
      f"El Niño (solid lines), El Niño + positive IOD (shading)",
      DirAndFile(out_dir, dir_results, 'esquema_z200', '1', 'pdf'), dpi, save,
      1, 2.3)


# Fase negativa -------------------------------------------------------------- #
ninia_obs = xr.open_dataset(data_dir +  'HGT200SON_N34_un_neg_mer_d_w_DMI_standard.nc')

dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
dmi = dmi.rename({'var_polyfit_coefficients':'var'})

sim = xr.open_dataset(data_dir + 'HGT200SON_DMI_sim_neg_mer_d_w_DMI_standard.nc')
dmi = xr.open_dataset(data_dir +  'HGT200SON_DMI_un_neg_mer_d_w_DMI_standard.nc')


Plot3(sim, [-2, -.5, 0, .5, 2], ninia_obs, [-.5, .5], dmi*0, [-.4, .4],
      'El Niño (lineas solidas), IOD pos. (lineas punteadas) ' +
      '\n' + 'El Niño + IOD positivo (sombreado)',
      DirAndFile(out_dir, dir_results, 'esquema', '1.5', 'pdf'), dpi, save, 1)


# # PP ------------------------------------------------------------------------- #
# # Fase positiva -------------------------------------------------------------- #
# ninio_mod = xr.open_dataset(data_dir + 'SNR_prec_CFSv2n34_puros_pos_SON.nc')
# ninio_obs = xr.open_dataset(data_dir + 'ppgpccSON_N34_un_pos_mer_d_w_DMI_standard.nc')
# ninio_obs = ninio_obs.interp(lon=ninio_mod.lon.values, lat=ninio_mod.lat.values)
# ninio_obs=xr.where(np.isnan(ninio_obs), 0, ninio_obs)
# # dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
# # dmi = dmi.rename({'var_polyfit_coefficients':'var'})
#
# sim_mod = xr.open_dataset(data_dir + 'SNR_prec_CFSv2sim_pos_SON.nc')
# sim_obs = xr.open_dataset(data_dir + 'ppgpccSON_DMI_sim_pos_mer_d_w_DMI_standard.nc')
# sim_obs = sim_obs.interp(lon=sim_mod.lon.values, lat=sim_mod.lat.values)
# sim_obs=xr.where(np.isnan(sim_obs), 0, sim_obs)
#
# ninio = (ninio_mod + 1.1*ninio_obs)/2
# sim = (sim_mod*2 + sim_obs)/2
#
# dmi = xr.open_dataset(data_dir + 'ppgpccSON_DMI_un_pos_mer_d_w_DMI_standard.nc')*0
#
# Plot3b(sim, [-5,-.5, -.2, 0, .2, .5, 5],
#       ninio, [-5,-.5, -.2, .2, .5, 5],
#       dmi, [-2, 2], '',
#       DirAndFile(out_dir, dir_results, 'esquema_pp', '2'),
#       dpi, save, 1, True,
#       smoothing=True, win1=(5,5), win2=(6,6))
#
#
# ####
#
