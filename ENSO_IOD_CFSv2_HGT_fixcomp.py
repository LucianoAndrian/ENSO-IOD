"""
Composiciones de HGT200 a partir de los outputs de ENSO_IOD_CFSv2_preSELECT_HGT.py
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
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/cfsv2/'
save = False
dpi = 100
# Funciones ############################################################################################################
def Plot(comp, comp_var, levels = np.linspace(-1,1,11),
         cmap='RdBu', dpi=100, save=True, step=1, name_fig='fig', title='title', color_map='grey'):

    import matplotlib.pyplot as plt
    levels_contour = levels.copy()
    levels_contour.remove(0)
    #comp_var = comp['var']
    fig = plt.figure(figsize=(9, 3.5), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([0,359, -80,20], crs_latlon)

    ax.contour(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step], linewidths=.8,
                     levels=levels_contour, transform=crs_latlon, colors='black')

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')

    ax.add_feature(cartopy.feature.LAND, facecolor='white', edgecolor=color_map)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(0, 360, 30), crs=crs_latlon)
    ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
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

def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w
########################################################################################################################
seasons = ['SON']
cases = ['dmi_puros_pos', 'dmi_puros_neg', 'n34_puros_pos', 'n34_puros_neg', 'sim_pos', 'sim_neg']

title_case = ['DMI pure - positive',
              'DMI pure - negative',
              'El Niño pure', 'La Niña pure',
              'DMI positive - El Niño',
              'DMI negative - La Niña']

cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')

cbar_snr = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#6FFE9B', '#FFFFFF',
                                  '#FFFFFF', '#FFFFFF',
                                  '#FEB77E','#CA3E72','#782281','#251255'])
cbar_snr.set_over('#251255')
cbar_snr.set_under('#070B4F')
cbar_snr.set_bad(color='white')

scale = [-300,-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250,300]
scale_cont= [-300,-250,-200,-150,-100,-50,-25,25,50,100,150,200,250,300]

for s in seasons:
    neutro = xr.open_dataset(cases_dir + 'hgt_neutros_' + s + '.nc').rename({'hgt':'var'})
    c_count = 0
    neutro = Weights(neutro.__mul__(9.80665))
    for c in cases:
        case = xr.open_dataset(cases_dir + 'hgt_' + c + '_' + s + '.nc').rename({'hgt':'var'})
        case = Weights(case.__mul__(9.80665))

        try:
            num_case = len(case.time)
            comp = case.mean('time') - neutro.mean('time')
            comp_var = comp['var']
            Plot(comp, comp_var, levels=scale,
                 cmap=cbar, dpi=dpi, step=1, name_fig='hgt_' + c + '_' + s,
                 title='Mean Composite - CFSv2 - ' + s + '\n' + title_case[c_count] + '\n' + ' ' + 'HGT 200hPa'
                       + ' - ' + 'Cases: ' + str(num_case),
                 save=save)

            spread = case - comp
            spread = spread.std('time')
            snr = comp / spread

            Plot(snr, snr['var'], levels = [-1,-.8,-.6,-.5,-.1,0,0.1,0.5,0.6,0.8,1],
                 cmap=cbar_snr, dpi=dpi, step=1, name_fig='SNR_hgt_' + c + '_' + s,
                 title='Signal-to-noise ratio - CFSv2 - ' + s + '\n' + title_case[c_count] + '\n' + ' ' + 'HGT 200hPa'
                       + ' - ' + 'Cases: ' + str(num_case),
                 save=save)

        except:
            print('Error in ' + c + ' - ' + s)

        c_count += 1
########################################################################################################################