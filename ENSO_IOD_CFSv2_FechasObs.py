"""
Composites de CFSv2 a partir de las fechas DMI N34 observado
"""
########################################################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings("ignore")
########################################################################################################################
SelectDates = False
save = True
dpi = 200
########################################################################################################################
obs_dates_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/' #fechas
cfsv2_cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'

seasons = ['JJA', 'SON']
cases = ['DMI_un_pos', 'DMI_un_neg', 'N34_un_pos', 'N34_un_neg', 'DMI_sim_pos', 'DMI_sim_neg',
         'Neutral']

def SelectEvents_ObsDates(variable, s):
    obs_dates = xr.open_dataset(obs_dates_dir + '1920_2020_' + s.upper() + '.nc')
    cfsv2_case = xr.open_dataset(cfsv2_cases_dir + variable + '_' + s.lower() + '.nc')

    for c in cases:
        cfsv2_select_case = cfsv2_case.sel(time=cfsv2_case.time.dt.year.isin(obs_dates[c].values))
        cfsv2_select_case.to_netcdf(out_dir + variable + '_' + c + '_' + s + '_CFSv2_obsDates.nc')

if SelectDates:
    # SST ##############################################################################################################
    SelectEvents_ObsDates('sst', 'JJA')
    SelectEvents_ObsDates('sst', 'SON')

    # HGT ##############################################################################################################
    SelectEvents_ObsDates('hgt', 'JJA')
    SelectEvents_ObsDates('hgt', 'SON')

    # PP ###############################################################################################################
    SelectEvents_ObsDates('prec', 'JJA')
    SelectEvents_ObsDates('prec', 'SON')

    # T ################################################################################################################

    def SelectEvents_ObsDates(variable, s):
        obs_dates = xr.open_dataset(obs_dates_dir + '1920_2020_' + s.upper() + '.nc')
        cfsv2_case = xr.open_dataset(cfsv2_cases_dir + variable + '_' + s.lower() + '_nodetrend.nc')

        for c in cases:
            cfsv2_select_case = cfsv2_case.sel(time=cfsv2_case.time.dt.year.isin(obs_dates[c].values))
            cfsv2_select_case.to_netcdf(out_dir + variable + '_' + c + '_' + s + '_CFSv2_obsDates.nc')

    SelectEvents_ObsDates('tref', 'JJA')
    SelectEvents_ObsDates('tref', 'SON')

########################################################################################################################
#-------------------------------------------------- Composites --------------------------------------------------------#
########################################################################################################################
seasons = ['JJA', 'SON']
cases = ['DMI_un_pos', 'DMI_un_neg', 'N34_un_pos', 'N34_un_neg', 'DMI_sim_pos', 'DMI_sim_neg']

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


cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')


cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')


cbar_snr = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#FEB77E', '#FB8761','#CA3E72','#782281','#251255'][::-1])
cbar_snr.set_over('#251255')
cbar_snr.set_under('#070B4F')
cbar_snr.set_bad(color='white')


cbar_Prueba = colors.ListedColormap(['#002A3D','#074D4F', '#1E6D5A' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#DCBC75', '#995D13','#6A3D07','#543005','#3F2404'][::-1])
cbar_Prueba.set_under('#3F2404')
cbar_Prueba.set_over('#002A3D')
cbar_Prueba.set_bad(color='white')


cbar_snr_t = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#FEB77E', '#FB8761','#CA3E72','#782281','#251255'])
cbar_snr_t.set_over('#251255')
cbar_snr_t.set_under('#070B4F')
cbar_snr_t.set_bad(color='white')
########################################################################################################################
########################################################################################################################
########################################################################################################################
#--- SST ---#
#----------------------------------------------------------------------------------------------------------------------#
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/CFSv2_ObsDate/Composites/SST/'
# Funciones -----------------------------------------------------------------------------------------------------------#
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
#----------------------------------------------------------------------------------------------------------------------#

scale = [-1.5, -1, -0.5, -0.2, 0, 0.2, 0.5, 1, 1.5]

for s in seasons:
    neutro = xr.open_dataset(cases_dir + 'sst_Neutral_' + s + '_CFSv2_obsDates.nc').rename({'sst':'var'})
    c_count = 0
    for c in cases:
        case = xr.open_dataset(cases_dir + 'sst_' + c + '_' + s + '_CFSv2_obsDates.nc').rename({'sst':'var'})

        try:
            num_case = len(case.time)*len(case.r)
            comp = case.mean(['time', 'r']) - neutro.mean(['time', 'r'])

            Plot(comp, levels=scale, cmap=cbar_sst, dpi=dpi, step=1,
                 name_fig='sst_' + c + '_' + s,
                 title='Mean Composite - CFSv2 - ' + s + '\n' + title_case[c_count] + '\n' + ' ' + 'SST'
                       + ' - ' + 'Cases: ' + str(num_case),
                 save=save)
        except:
            print('Error in ' + c + ' - ' + s)

        c_count += 1
########################################################################################################################
########################################################################################################################
########################################################################################################################
#--- HGT ---#
#----------------------------------------------------------------------------------------------------------------------#
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/CFSv2_ObsDate/Composites/HGT/'

scale = [-300,-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250,300]
scale_cont= [-300,-250,-200,-150,-100,-50,-25,25,50,100,150,200,250,300]
# Funciones -----------------------------------------------------------------------------------------------------------#
def Plot(comp, comp_var, levels = np.linspace(-1,1,11),
         cmap='RdBu', dpi=100, save=True, step=1, name_fig='fig', title='title', color_map='#4B4B4B'):

    import matplotlib.pyplot as plt
    levels_contour = levels.copy()
    levels_contour.remove(0)
    #comp_var = comp['var']
    fig = plt.figure(figsize=(8, 3), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([30,340, -80,20], crs_latlon)

    ax.contour(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step], linewidths=.8,
                     levels=levels_contour, transform=crs_latlon, colors='black')

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')

    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey', edgecolor=color_map)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(30, 340, 25), crs=crs_latlon)
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

#----------------------------------------------------------------------------------------------------------------------#

for s in seasons:
    neutro = xr.open_dataset(cases_dir + 'hgt_Neutral_' + s + '_CFSv2_obsDates.nc').rename({'hgt':'var'})
    c_count = 0
    neutro = neutro.__mul__(9.80665)
    for c in cases:
        case = xr.open_dataset(cases_dir + 'hgt_' + c + '_' + s + '_CFSv2_obsDates.nc').rename({'hgt':'var'})
        case = case.__mul__(9.80665)

        try:
            num_case = len(case.time)*len(case.r)
            comp = case.mean(['time', 'r']) - neutro.mean(['time', 'r'])
            comp_var = comp['var']
            Plot(comp, comp_var, levels=scale,
                 cmap=cbar, dpi=dpi, step=1, name_fig='hgt_' + c + '_' + s,
                 title='Mean Composite - CFSv2 - ' + s + '\n' + title_case[c_count] + '\n' + ' ' + 'HGT 200hPa'
                       + ' - ' + 'Cases: ' + str(num_case),
                 save=save)
        except:
            print('Error in ' + c + ' - ' + s)

        c_count += 1
########################################################################################################################
########################################################################################################################
########################################################################################################################
#--- PP ---#
#----------------------------------------------------------------------------------------------------------------------#
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/CFSv2_ObsDate/Composites/PP/'

# colorbars -----------------------------------------------------------------------------------------------------------#
scale_signal = np.linspace(-30, 30, 13)
scale_snr = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]
scale_prob = [.2,.3,.4,.45,.5,.55,.6,.7,.8]

# Funciones -----------------------------------------------------------------------------------------------------------#
def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt

    comp_var = comp['var']
    fig = plt.figure(figsize=(5, 6), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270,330, -60,20], crs_latlon)

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
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

detrend=True
from ENSO_IOD_Funciones import MakeMask

for s in seasons:
    if detrend:
        file_name_end = '_CFSv2_obsDates.nc'
        name_fig_end = ''
        title_trend = 'Detrend'
    else:
        file_name_end = '_nodetrend_CFSv2_obsDates.nc'
        name_fig_end = 'NoDetrend'
        title_trend = 'No Detrend'

    neutro = xr.open_dataset(cases_dir + 'prec_Neutral_' + s + file_name_end).rename({'prec':'var'})
    neutro *= 30 #mm/day -> mm/month
    c_count = 0
    for c in cases:
        case = xr.open_dataset(cases_dir + 'prec_' +  c + '_' + s + file_name_end).rename({'prec':'var'})
        case *= 30
        try:
            num_case = len(case.time)

            # signal (comp)
            num_case = len(case.time)*len(case.r)
            comp = case.mean(['time', 'r']) - neutro.mean(['time', 'r'])
            mask = MakeMask(comp, 'var')
            comp *= mask

            # Plot(comp, levels=scale_signal, cmap=cbar_pp, dpi=dpi, step=1,
            #      name_fig='prec_comp_' + c + '_' + s + name_fig_end,
            #      title='Mean Composite - CFSv2 - ' + s + '\n'
            #            + title_case[c_count] + '\n' + ' ' + 'PREC'
            #            + ' - ' + 'Cases: ' + str(num_case) + ' ' + name_fig_end,
            #      save=save)

            # # spread
            spread = case - comp
            spread = spread.std(['time', 'r'])
            spread *= mask

            # snr
            snr = comp/spread
            Plot(snr, levels=scale_snr, cmap=cbar_Prueba, dpi=dpi, step=1,
                 name_fig='prec_SNR_' + c + '_' + s + '_' + name_fig_end,
                 title='Singal-to-noise ratio - CFSv2 - ' + s + '\n'
                       + title_case[c_count] + '\n' + ' ' + 'PREC'
                       + ' - ' + 'Cases: ' + str(num_case) + ' ' + name_fig_end,
                 save=save)

            #prob
            # comp_prob = case - neutro.mean('time')
            # aux = xr.where(comp_prob > 0, comp_prob, np.nan)
            # prob = SpatialProbability(aux, mask)
            # Plot(prob, levels=scale_prob,cmap=cbar_Prueba, dpi=dpi, step=1,
            #      name_fig='prec_prob_' + c + '_' + s + '_' + name_fig_end,
            #      title='Probability of PP>0' + '- CFSv2 - ' + s + '\n'
            #            + title_case[c_count] + '\n' + ' ' + 'PREC'
            #            + ' - ' + 'Cases: ' + str(num_case) + ' ' + name_fig_end,
            #      save=save)

        except:
            print('Error in ' + c + ' - ' + s)

        c_count += 1

########################################################################################################################
########################################################################################################################
########################################################################################################################
#--- T ---#
#----------------------------------------------------------------------------------------------------------------------#
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/CFSv2_ObsDate/Composites/T/'

# colorbars -----------------------------------------------------------------------------------------------------------#
scale_signal =  np.linspace(-1.2,1.2,13)
scale_snr = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]
scale_prob = [.2,.3,.4,.45,.5,.55,.6,.7,.8]


for s in seasons:
    neutro = xr.open_dataset(cases_dir + 'tref_Neutral_' + s + '_CFSv2_obsDates.nc').rename({'tref':'var'})
    c_count = 0
    for c in cases:
        case = xr.open_dataset(cases_dir + 'tref_' +  c + '_' + s + '_CFSv2_obsDates.nc').rename({'tref':'var'})
        try:
            num_case = len(case.time)
            # signal (comp)
            # signal (comp)
            num_case = len(case.time)*len(case.r)
            comp = case.mean(['time', 'r']) - neutro.mean(['time', 'r'])
            mask = MakeMask(comp, 'var')
            comp *= mask

            # Plot(comp, levels=scale_signal, cmap=cbar_t, dpi=dpi, step=1,
            #      name_fig='tref_comp_' + c + '_' + s,
            #      title='Mean Composite - CFSv2 - ' + s + '\n'
            #            + title_case[c_count] + '\n' + ' ' + 'TREF'
            #            + ' - ' + 'Cases: ' + str(num_case),
            #      save=save)

            # spread
            spread = case - comp
            spread = spread.std(['time','r'])
            spread *= mask

            # # snr
            snr = comp/spread
            Plot(snr, levels=scale_snr, cmap=cbar_snr_t, dpi=dpi, step=1,
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

        c_count += 1
########################################################################################################################
########################################################################################################################