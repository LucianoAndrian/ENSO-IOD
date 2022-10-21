"""
Composiciones de SST a partir de los outputs de ENSO_IOD_CFSv_fixSELECT_variables.py
"""
########################################################################################################################
import xarray as xr
import numpy as np
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings('ignore')
########################################################################################################################
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/SST/'
save = True
dpi = 200
# Funciones ############################################################################################################
def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title', color_map='#4B4B4B'):

    import matplotlib.pyplot as plt

    levels_contour = levels.copy()
    comp_var = comp['var']
    if isinstance(levels, np.ndarray):
        levels_contour = levels[levels != 0]
    else:
        levels_contour.remove(0)
    fig = plt.figure(figsize=(7, 2), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([50,270, -20,20], crs_latlon)
    # ax.contour(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
    #            linewidths=.8, levels=levels_contour, transform=crs_latlon, colors='black')

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey', edgecolor=color_map)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(50, 270, 60), crs=crs_latlon)
    ax.set_yticks(np.arange(-20, 40, 20), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)

    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

########################################################################################################################
seasons = ['JJA', 'JAS', 'ASO', 'SON']
cases = ['dmi_puros_pos', 'dmi_puros_neg', 'n34_puros_pos', 'n34_puros_neg', 'sim_pos', 'sim_neg']

title_case = ['DMI pure - positive',
              'DMI pure - negative',
              'El Niño pure', 'La Niña pure',
              'DMI positive - El Niño',
              'DMI negative - La Niña']

cbar_sst = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_sst.set_over('#9B1C00')
cbar_sst.set_under('#014A9B')
cbar_sst.set_bad(color='white')

scale = [-1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5]

for s in seasons:
    neutro = xr.open_dataset(cases_dir + 'neutros_' + s + '.nc').rename({'sst':'var'})
    c_count = 0
    for c in cases:
        case = xr.open_dataset(cases_dir + c + '_' + s + '.nc').rename({'sst':'var'})

        try:
            num_case = len(case.time)
            comp = case.mean('time') - neutro.mean('time')

            Plot(comp, levels=scale, cmap=cbar_sst, dpi=dpi, step=1,
                 name_fig='sst_' + c + '_' + s,
                 title='Mean Composite - CFSv2 - ' + s + '\n' + title_case[c_count] + '\n' + ' ' + 'SST'
                       + ' - ' + 'Cases: ' + str(num_case),
                 save=save)
        except:
            print('Error in ' + c + ' - ' + s)

        c_count += 1
########################################################################################################################