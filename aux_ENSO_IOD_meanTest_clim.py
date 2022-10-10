"""
t-test de medias para las climatologias de las variables HGT
"""
########################################################################################################################
import xarray as xr
import numpy as np
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
from ENSO_IOD_Funciones import SelectNMMEFiles
# Funciones ############################################################################################################
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

def TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, main_month_season):

    for l in [0,1,2,3]:
        season_1982_1998 = data_1982_1998.sel(time=data_1982_1998.time.dt.month.isin(main_month_season-l), L=l)
        season_1999_2011 = data_1999_2011.sel(time=data_1999_2011.time.dt.month.isin(main_month_season-l), L=l)

        if l==0:
            season_clim_1982_1998 = season_1982_1998.mean(['r', 'time'])
            season_clim_1999_2011 = season_1999_2011.mean(['r', 'time'])

            season_data_1982_1998 = season_1982_1998
            season_data_1999_2011 = season_1999_2011
            #
            # season_anom_1982_1998 = season_1982_1998 - season_clim_1982_1998
            # season_anom_1999_2011 = season_1999_2011 - season_clim_1999_2011
        else:
            season_clim_1982_1998 = xr.concat([season_clim_1982_1998, season_1982_1998.mean(['r', 'time'])], dim='L')
            season_clim_1999_2011 = xr.concat([season_clim_1999_2011, season_1999_2011.mean(['r', 'time'])], dim='L')

            season_data_1982_1998 = xr.concat([season_data_1982_1998, season_1982_1998], dim='L')
            season_data_1999_2011 = xr.concat([season_data_1999_2011, season_1999_2011], dim='L')
            # aux_1982_1998 = season_1982_1998 - season_1982_1998.mean(['r', 'time'])
            # aux_1999_2011 = season_1999_2011 - season_1999_2011.mean(['r', 'time'])
            #
            # season_anom_1982_1998 = xr.concat([season_anom_1982_1998, aux_1982_1998], dim='time')
            # season_anom_1999_2011 = xr.concat([season_anom_1999_2011, aux_1999_2011], dim='time')
    #
    # season_1982_1998 = data_1982_1998.sel(time=data_1982_1998.time.dt.month.isin(main_month_season))
    # season_1999_2011 = data_1999_2011.sel(time=data_1999_2011.time.dt.month.isin(main_month_season))
    #
    # season_clim_1982_1998 = season_1982_1998.mean(['r', 'time'])
    # season_clim_1999_2011 = season_1999_2011.mean(['r', 'time'])
    #
    # season_anom_1982_1998 = season_1982_1998 - season_clim_1982_1998
    # season_anom_1999_2011 = season_1999_2011 - season_clim_1999_2011

    return season_clim_1982_1998, season_clim_1999_2011,  season_data_1982_1998, season_data_1999_2011#, season_anom_1982_1998, season_anom_1999_2011

def Plot(data_plot, data_var, levels = np.linspace(-3,3,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt

    #comp_var = comp['var']
    fig = plt.figure(figsize=(8, 3), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([30,340, -80,20], crs_latlon)

    im = ax.contourf(data_plot.lon[::step], data_plot.lat[::step], data_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
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

########################################################################################################################
dir_hc = '/pikachu/datos/luciano.andrian/hindcast/'
dir_rt = '/pikachu/datos/luciano.andrian/real_time/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/HGT/diff_clim/'
v = 'hgt'

save = True
########################################################################################################################
# Climatologias 1982-1998 1999-2011 del Hindcast ----------------------------------------------------------------------#
# debido al salto en la climatologia en SST

# usando SelectNMMEFiles con All=True,
# abre TODOS los archivos .nc de la ruta en dir

### HINDCAST ###
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_hc, All=True)
files = sorted(files, key=lambda x: x.split()[0])

#abriendo todos los archivos
data = xr.open_mfdataset(files, decode_times=False) #xr no entiende la codificacion de Leads, r y las fechas
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data = data.sel(L=[0.5, 1.5, 2.5, 3.5]) # Solo leads 0 1 2 3
data['L'] = [0,1,2,3]
data = xr.decode_cf(fix_calendar(data)) # corrigiendo fechas
data = data.sel(lon=slice(30, 340), lat=slice(-80, 20), P=200)
data = data.drop('P')

#media movil de 3 meses para separar en estaciones
data = data.rolling(time=3, center=True).mean()

# 1982-1998, 1999-2011
data_1982_1998 = data.sel(time=data.time.dt.year.isin(np.linspace(1982,1998,17)))
data_1999_2011 = data.sel(time=data.time.dt.year.isin(np.linspace(1999,2011,13)))

n1 = 1998-1982+1
n2 = 2011-1999+1

#---------------------- Climatologias, anomalias y detrend por ESTACIONES JJA, JAS, ASO y SON -------------------------#
#--- JJA --------------------------------------------------------------------------------------------------------------#
print('JJA')
jja_clim_82_98, jja_clim_99_11, jja_data_82_98, jja_data_99_11 = TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 7)
jja_clim_82_98.load()
jja_clim_99_11.load()
jja_data_99_11.load()
jja_data_82_98.load()

aux_82_98_var = jja_data_82_98.var(['time', 'r'])
aux_99_11_var = jja_data_99_11.var(['time', 'r'])

t_jja = (jja_clim_82_98 - jja_clim_99_11)/np.sqrt((aux_82_98_var/n1)+(aux_99_11_var/n2))
aux_df = (((aux_82_98_var/n1)+(aux_99_11_var/n2))**2)/((((aux_82_98_var/n1)**2)/(n1-1))+(((aux_99_11_var/n2)**2)/(n2-1)))

from scipy.stats import t
t_cr=t.ppf(0.95, aux_df.hgt)
aux = t_jja.where(abs(t_jja) > t_cr)

for l in [0,1,2,3]:
    Plot(aux, aux.hgt[l, :, :], cmap='RdBu_r', title='t JJA - Lead ' + str(l) + ' - HGT200',
         name_fig='HGT200_t_test_jja_' + str(l), save=save)

del jja_data_82_98, jja_data_99_11, jja_clim_82_98, jja_clim_99_11
#--- JAS --------------------------------------------------------------------------------------------------------------#
print('JAS')
jas_clim_82_98, jas_clim_99_11, jas_data_82_98, jas_data_99_11 = TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 8)
jas_clim_82_98.load()
jas_clim_99_11.load()
jas_data_99_11.load()
jas_data_82_98.load()

aux_82_98_var = jas_data_82_98.var(['time', 'r'])
aux_99_11_var = jas_data_99_11.var(['time', 'r'])

t_jas = (jas_clim_82_98 - jas_clim_99_11)/np.sqrt((aux_82_98_var/n1)+(aux_99_11_var/n2))
aux_df = (((aux_82_98_var/n1)+(aux_99_11_var/n2))**2)/((((aux_82_98_var/n1)**2)/(n1-1))+(((aux_99_11_var/n2)**2)/(n2-1)))

from scipy.stats import t
t_cr=t.ppf(0.95, aux_df.hgt)
aux = t_jas.where(abs(t_jas) > t_cr)

for l in [0,1,2,3]:
    Plot(aux, aux.hgt[l, :, :], cmap='RdBu_r', title='t JAS- Lead ' + str(l) + ' - HGT200',
         name_fig='HGT200_t_test_jas_' + str(l), save=save)


del jas_data_82_98, jas_data_99_11, jas_clim_82_98, jas_clim_99_11

#--- ASO --------------------------------------------------------------------------------------------------------------#
print('ASO')
aso_clim_82_98, aso_clim_99_11, aso_data_82_98, aso_data_99_11 = TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 9)
aso_clim_82_98.load()
aso_clim_99_11.load()
aso_data_99_11.load()
aso_data_82_98.load()

aux_82_98_var = aso_data_82_98.var(['time', 'r'])
aux_99_11_var = aso_data_99_11.var(['time', 'r'])

t_aso = (aso_clim_82_98 - aso_clim_99_11)/np.sqrt((aux_82_98_var/n1)+(aux_99_11_var/n2))
aux_df = (((aux_82_98_var/n1)+(aux_99_11_var/n2))**2)/((((aux_82_98_var/n1)**2)/(n1-1))+(((aux_99_11_var/n2)**2)/(n2-1)))

from scipy.stats import t
t_cr=t.ppf(0.95, aux_df.hgt)
aux = t_aso.where(abs(t_aso) > t_cr)

for l in [0,1,2,3]:
    Plot(aux, aux.hgt[l, :, :], cmap='RdBu_r', title='t ASO- Lead ' + str(l) + ' - HGT200',
         name_fig='HGT200_t_test_aso_' + str(l), save=save)


del aso_data_82_98, aso_data_99_11, aso_clim_82_98, aso_clim_99_11
#--- SON --------------------------------------------------------------------------------------------------------------#
print('SON')
son_clim_82_98, son_clim_99_11, son_data_82_98, son_data_99_11 = TwoClim_Anom_Seasons(data_1982_1998, data_1999_2011, 10)
son_clim_82_98.load()
son_clim_99_11.load()
son_data_99_11.load()
son_data_82_98.load()

aux_82_98_var = son_data_82_98.var(['time', 'r'])
aux_99_11_var = son_data_99_11.var(['time', 'r'])

t_son = (son_clim_82_98 - son_clim_99_11)/np.sqrt((aux_82_98_var/n1)+(aux_99_11_var/n2))
aux_df = (((aux_82_98_var/n1)+(aux_99_11_var/n2))**2)/((((aux_82_98_var/n1)**2)/(n1-1))+(((aux_99_11_var/n2)**2)/(n2-1)))

from scipy.stats import t
t_cr=t.ppf(0.95, aux_df.hgt)
aux = t_son.where(abs(t_son) > t_cr)

for l in [0,1,2,3]:
    Plot(aux, aux.hgt[l, :, :], cmap='RdBu_r', title='t SON- Lead ' + str(l) + ' - HGT200',
         name_fig='HGT200_t_test_son_' + str(l), save=save)


del son_data_82_98, son_data_99_11, son_clim_82_98, son_clim_99_11