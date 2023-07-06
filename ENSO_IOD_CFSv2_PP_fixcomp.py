"""
Composiciones de PP a partir de los outputs de
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
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/PP/'
save = True
dpi = 300
detrend = True
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
                     comp_var[::step, ::step], levels=levels,
                     transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.BORDERS, facecolor='k')
    ax.add_feature(cartopy.feature.OCEAN, zorder=10, facecolor='white',
                   edgecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', zorder=12)
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    # al usar ocean como mascara
    # y las girdlanes encima
    # los bordes del plot quedan tapadaos por las gridlines
    for k, spine in ax.spines.items():
        spine.set_zorder(13)

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
                len(data['var'][:, lt, ln][~np.isnan(
                    data['var'][:, lt, ln].values)].values) \
                / len(data['var'][:, lt, ln])
    return prob*mask
################################################################################
seasons = ['JJA','SON']
cases = ['dmi_puros_pos', 'dmi_puros_neg', 'n34_puros_pos', 'n34_puros_neg',
         'sim_pos', 'sim_neg']

title_case = ['DMI pure - positive',
              'DMI pure - negative',
              'El Niño pure', 'La Niña pure',
              'DMI positive - El Niño',
              'DMI negative - La Niña']
# colorbars -------------------------------------------------------------------#
cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC',
                                 '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07',
                                 '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cbar_snr = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#52C39D',
                                  '#6FFE9B',
                                  '#FFFFFF',
                                  '#FEB77E', '#FB8761','#CA3E72','#782281',
                                  '#251255'][::-1])
cbar_snr.set_over('#251255')
cbar_snr.set_under('#070B4F')
cbar_snr.set_bad(color='white')

cbar_Prueba = colors.ListedColormap(['#002A3D','#074D4F', '#1E6D5A' ,'#52C39D',
                                     '#6FFE9B',
                                  '#FFFFFF',
                                  '#DCBC75', '#995D13','#6A3D07','#543005',
                                     '#3F2404'][::-1])
cbar_Prueba.set_under('#3F2404')
cbar_Prueba.set_over('#002A3D')
cbar_Prueba.set_bad(color='white')
#------------------------------------------------------------------------------#

scale_signal = np.linspace(-30, 30, 13)
scale_snr = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]
scale_prob = [.2,.3,.4,.45,.5,.55,.6,.7,.8]

for s in seasons:
    if detrend:
        file_name_end = '.nc'
        name_fig_end = ''
        title_trend = 'Detrend'
    else:
        file_name_end = '_nodetrend.nc'
        name_fig_end = 'NoDetrend'
        title_trend = 'No Detrend'

    neutro = xr.open_dataset(
        cases_dir + 'prec_neutros_' + s + file_name_end).rename({'prec':'var'})
    neutro *= 30 #mm/day -> mm/month
    c_count = 0
    for c, c_count in zip(cases, range(0, len(cases))):
        case = xr.open_dataset(
            cases_dir + 'prec_' +  c + '_' + s +
            file_name_end).rename({'prec':'var'})
        case *= 30
        try:
            num_case = len(case.time)

            # signal (comp)
            comp = case.mean('time') - neutro.mean('time')
            mask = MakeMask(comp, 'var')
            # comp *= mask

            Plot(comp, levels=scale_signal, cmap=cbar_pp, dpi=dpi, step=1,
                 name_fig='prec_comp_' + c + '_' + s + name_fig_end,
                 title='Mean Composite - CFSv2 - ' + s + '\n'
                       + title_case[c_count] + '\n' + ' ' + 'PREC'
                       + ' - ' + 'Cases: ' + str(num_case) + ' ' + name_fig_end,
                 save=save)

            # # spread
            spread = case - comp
            spread = spread.std('time')
            # spread *= mask

            # snr
            snr = comp/spread
            Plot(snr, levels=scale_snr, cmap=cbar_Prueba, dpi=dpi, step=1,
                 name_fig='prec_SNR_' + c + '_' + s + '_' + name_fig_end,
                 title='Singal-to-noise ratio - CFSv2 - ' + s + '\n'
                       + title_case[c_count] + '\n' + ' ' + 'PREC'
                       + ' - ' + 'Cases: ' + str(num_case) + ' ' + name_fig_end,
                 save=save)

            #prob
            comp_prob = case - neutro.mean('time')
            aux = xr.where(comp_prob > 0, comp_prob, np.nan)
            prob = SpatialProbability(aux, mask)
            Plot(prob, levels=scale_prob,cmap=cbar_Prueba, dpi=dpi, step=1,
                 name_fig='prec_prob_' + c + '_' + s + '_' + name_fig_end,
                 title='Probability of PP>0' + '- CFSv2 - ' + s + '\n'
                       + title_case[c_count] + '\n' + ' ' + 'PREC'
                       + ' - ' + 'Cases: ' + str(num_case) + ' ' + name_fig_end,
                 save=save)

        except:
            print('Error in ' + c + ' - ' + s)
################################################################################