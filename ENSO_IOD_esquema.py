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
warnings.filterwarnings("ignore", category=ShapelyDeprecationWarning)
################################################################################
save = False
test = True
CreateDirectory(out_dir, dir_results)
#------------------------------------------------------------------------------#
if save:
    dpi = 300
else:
    dpi = 100
################################################################################
# Test: Que resutla más adecuado?
ninio_obs = xr.open_dataset(data_dir +  'HGT200SON_N34_un_pos_mer_d_w_DMI_standard.nc')

dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
dmi = dmi.rename({'var_polyfit_coefficients':'var'})

sim = xr.open_dataset(data_dir + 'HGT200SON_DMI_sim_pos_mer_d_w_DMI_standard.nc')

def Plot3(data = None, levels = None, data2 = None, levels2 = None,
          data3 = None, levels3 = None, title = 'title', name = 'name',
          dpi = 70, save = False, step = 1):
    crs_latlon = ccrs.PlateCarree()

    fig = plt.figure(figsize=(9, 3.5), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_extent([0, 359, -80, 20], crs=crs_latlon)

    try:
        data_var = data['var']
        ax.contourf(data.lon[::step], data.lat[::step], data_var[::step, ::step],
                    linewidths=1, alpha=.4, levels=levels,
                    transform=crs_latlon, colors=['dodgerblue', 'white', 'white', 'firebrick'])
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
    ax.set_xticks(np.arange(0, 360, 30), crs=crs_latlon)
    ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color='grey', linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', zorder=20)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)
    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        print('save: ' + name)
        plt.savefig(name)
        plt.close()

    else:
        plt.show()



def Plot(data = None, levels = None, data2 = None,
         levels2 = None, levels2_cont = None, data3 = None, levels3 = None, title = 'title', name = 'name',
          dpi = 70, save = False, step = 1):
    crs_latlon = ccrs.PlateCarree()

    fig = plt.figure(figsize=(9, 3.5), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.set_extent([0, 359, -80, 20], crs=crs_latlon)


    try:
        data3_var = data3['var']
        ax.contourf(data3.lon[::step], data3.lat[::step], data3_var[::step, ::step],
                    linewidths=1, alpha=.4, levels=levels3,
                    transform=crs_latlon,
                    colors=['dodgerblue', 'white', 'white', 'firebrick'])
    except:
        pass

    try:
        data_var = data['var']
        ax.contour(data.lon[::step], data.lat[::step], data_var[::step, ::step],
                   linewidths=1.5, alpha=.7, levels=levels,
                    transform=crs_latlon, colors=['Blue', 'Red'])
    except:
        pass

    try:
        data2_var = data2['var']
        intersection_var = np.where((data_var >= levels[0]) &
                                    (data_var <= levels[-1]),
                                    np.nan, data2_var)
        ax.contourf(data.lon[::step], data.lat[::step],
                   intersection_var[::step, ::step], alpha=.4, levels=levels2,
                   transform=crs_latlon, colors=['#6D8BAB', 'white', 'white', '#BD8A8A'])
#            ['#00B0CC', '#D00038'])

        ax.contour(data2.lon[::step], data2.lat[::step],
                    data2_var[::step, ::step],
                    linewidths=1.5, alpha=1, levels=levels2_cont,
                    linestyles='dashed',
                    transform=crs_latlon,
                    colors=['#00B0CC', '#D00038'])
    except:
        pass



    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor='grey')
    ax.set_xticks(np.arange(0, 360, 30), crs=crs_latlon)
    ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color='grey', linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', zorder=20)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)
    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        print('save: ' + name)
        plt.savefig(name)
        plt.close()

    else:
        plt.show()


ninio_obs = xr.open_dataset(data_dir +  'HGT200SON_N34_un_pos_mer_d_w_DMI_standard.nc')

dmi = xr.open_dataset(data_dir +  'HGT200_SON_mer_d_w_SON_1940_2020_DMI_woN34.nc')
dmi = dmi.rename({'var_polyfit_coefficients':'var'})

sim = xr.open_dataset(data_dir + 'HGT200SON_DMI_sim_pos_mer_d_w_DMI_standard.nc')
dmi = xr.open_dataset(data_dir +  'HGT200SON_DMI_un_pos_mer_d_w_DMI_standard.nc')


Plot3(sim, [-2, -.5, 0, .5, 2], ninio_obs, [-.5, .5], dmi, [-.6, .6],
      'El Niño (lineas solidas), IOD pos. (lineas punteadas) ' +
      '\n' + 'El Niño + IOD positivo (sombreado)',
      DirAndFile(out_dir, dir_results, 'esquema', '1.5'), dpi, save, 1)


Plot(data = sim, levels = [-.5,.5], data2 = dmi,
      levels2 = [-500, -40, 0, 40, 500], levels2_cont=[-40, 40],
     title='El Niño (lineas solidas), IOD pos. (lineas punteadas) ' + '\n' +
           'El Niño + IOD positivo (sombreado)',
      name=DirAndFile(out_dir, dir_results, 'esquema', '2'), dpi=dpi, save=save,
      step=1)


