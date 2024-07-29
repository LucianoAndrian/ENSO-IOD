"""
Composiciones de T a partir de los outputs de
ENSO_IOD_CFSv_fixSELECT_variables.py
"""
################################################################################
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

from ENSO_IOD_Funciones import MakeMask
################################################################################
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/T/'

save = True
dpi = 300
detrend = True
save_nc = True
plot = False
# Funciones ####################################################################
def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt

    comp_var = comp['var']
    fig = plt.figure(figsize=(5, 6), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270,330, -60,20], crs_latlon)

    im = ax.contourf(comp.lon[::step], comp.lat[::step],
                     comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap,
                     extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
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

def SpatialProbability(data, mask):
    prob = xr.where(np.isnan(mask), mask, 1)
    for ln in range(0, 56):
        for lt in range(0, 76):
            prob['var'][lt, ln] = \
                len(data['var'][:, lt, ln]
                    [~np.isnan(data['var'][:, lt, ln].values)].values) \
                / len(data['var'][:, lt, ln])
    return prob*mask
################################################################################
seasons = ['SON']
cases = ['dmi_puros_pos', 'dmi_puros_neg', 'n34_puros_pos', 'n34_puros_neg',
         'sim_pos', 'sim_neg',
         'n34_pos', 'n34_neg', 'dmi_pos', 'dmi_neg']

title_case = ['DMI pure - positive',
              'DMI pure - negative',
              'El Niño pure', 'La Niña pure',
              'DMI positive - El Niño',
              'DMI negative - La Niña',
              'El Niño', 'La Niña',
              'DMI positive', 'DMI negative']

# colorbars --------------------------------------------------------------------
cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89',
                                '#FFCECC', 'white', '#B3DBFF', '#83B9EB',
                                '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')

cbar_snr = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#52C39D',
                                  '#6FFE9B','#FFFFFF', '#FEB77E', '#FB8761',
                                  '#CA3E72','#782281','#251255'])
cbar_snr.set_over('#251255')
cbar_snr.set_under('#070B4F')
cbar_snr.set_bad(color='white')

#-------------------------------------------------------------------------------

scale_signal =  [-1, -.75, -.5, -.25, -.1, 0, .1, .25, .5, .75, 1]
scale_snr = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]
scale_prob = [.2,.3,.4,.45,.5,.55,.6,.7,.8]


for s in seasons:
    if detrend:
        file_name_end = '.nc'
        name_fig_end = ''
        title_trend = 'Detrend'
    else:
        file_name_end = '.nc'
        name_fig_end = '_NoDetrend'
        title_trend = 'No Detrend'


    neutro = xr.open_dataset(
        cases_dir + 'tref_neutros_' + s + name_fig_end + '.nc')\
        .rename({'tref':'var'})

    for c, c_count in zip(cases, range(0, len(cases))):
        case = xr.open_dataset(
            cases_dir + 'tref_' +  c + '_' + s + name_fig_end + '.nc')\
            .rename({'tref':'var'})
        try:
            num_case = len(case.time)
            # signal (comp)
            comp = case.mean('time') - neutro.mean('time')
            mask = MakeMask(comp, 'var')
            comp *= mask

            Plot(comp, levels=scale_signal, cmap=cbar_t, dpi=dpi, step=1,
                 name_fig='tref_comp_' + c + '_' + s,
                 title='Mean Composite - CFSv2 - ' + s + '\n'
                       + title_case[c_count] + '\n' + ' ' + 'TREF'
                       + ' - ' + 'Cases: ' + str(num_case),
                 save=save)

            # spread
            spread = case - comp
            spread = spread.std('time')
            spread *= mask

            # # snr
            snr = comp/spread
            Plot(snr, levels=scale_snr, cmap=cbar_snr, dpi=dpi, step=1,
                 name_fig='tref_SNR_' + c + '_' + s,
                 title='Singal-to-noise ratio - CFSv2 - ' + s + '\n'
                       + title_case[c_count] + '\n' + ' ' + 'TREF'
                       + ' - ' + 'Cases: ' + str(num_case),
                 save=save)

            # #prob
            # comp_prob = case - neutro.mean('time')
            # aux = xr.where(comp_prob > 0, comp_prob, np.nan)
            # prob = SpatialProbability(aux, mask)
            # Plot(prob, levels=scale_prob,cmap=cbar_t, dpi=dpi, step=1,
            #      name_fig='tref_prob_' + c + '_' + s,
            #      title='Probability of T>0' + '- CFSv2 - ' + s + '\n'
            #            + title_case[c_count] + '\n' + ' ' + 'TREF'
            #            + ' - ' + 'Cases: ' + str(num_case),
            #      save=save)

        except:
            print('Error in ' + c + ' - ' + s)
################################################################################