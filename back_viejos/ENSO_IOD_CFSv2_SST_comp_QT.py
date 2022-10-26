import xarray as xr
import numpy as np
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
#warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
warnings.filterwarnings('ignore')

########################################################################################################################
cases_dir = '/pikachu/datos/luciano.andrian/cases/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/SNR_CFSv2/SST/QT/'
save=True
dpi=200
# Funciones ############################################################################################################
def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt

    comp_var = comp['var']
    fig = plt.figure(figsize=(7, 2), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([50,270, -20,20], crs_latlon)

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(50, 270, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-20, 20, 5), crs=crs_latlon)
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
seasons = ('JJA', 'JAS', 'ASO', 'SON')
set_name= ('3')

cases = ['sim_pos', 'sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg', 'sim_DMIneg_N34pos', 'sim_DMIpos_N34neg']

title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI negative phase ',
              'DMI positive phase ',
              'DMI isolated positive phase ',
              'DMI isolated negative phase ',
              'ENSO positive phase ',
              'ENSO negative phase ',
              'ENSO isolated positive phase ',
              'ENSO isolated negative phase ',
              'DMI negative and ENSO positive',
              'DMI positive and ENSO negative']

scale = [np.linspace(-30,30,13), np.linspace(-1,1,13)]
scale_sd = [np.linspace(0,10,9), np.linspace(0,1,11)]

cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')

scales = np.linspace(-2.4, 2.4, 13)

# v = 'prec'
# fix_factor=30
# s='SON'
# c='N34_un_neg'
# s_n= '0'
v_count = 0
v = 'sst'
for s in seasons:
    print(s)
    for s_n in set_name:
        print(s_n)
        c_count = 0
        for c in cases:
            print(c)
            data_neutral = xr.open_dataset(cases_dir + v + '_QT' + s + '_Set' + s_n + '_NEUTRO.nc').drop(['L', 'r'])
            data_neutral = data_neutral.rename({v: 'var'})

            try:
                data_case = xr.open_dataset(cases_dir + v + '_QT' + s + '_Set' + s_n + '_' + c + '.nc').drop(['L', 'r'])
                data_case = data_case.rename({v: 'var'})

                num_case = len(data_case.time)

                # signal
                aux = data_case.mean('time') - data_neutral.mean('time')
                Plot(aux, levels=np.linspace(-1.2,1.2,13), cmap=cbar_t,
                     dpi=dpi, step=1,
                     name_fig=v + '_QT_set_' + s_n + '_' + c + '_' + s + '.nc',
                     title='Mean Composite - CFSv2 - ' + s + '\n' + title_case[c_count] + '\n' + ' ' + v + ' - ' +
                           'Leads: ' + str(s_n) + ' - ' + 'Casos: ' + str(num_case),
                     save=save)
            except:
                print('No' + c + ' in Set' + s_n + ' at ' + s)

            c_count += 1
v_count += 1

