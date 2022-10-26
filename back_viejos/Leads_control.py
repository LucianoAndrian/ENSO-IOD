import xarray as xr
import numpy as np
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

########################################################################################################################
models_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/'
data_ref_dir = '/datos/luciano.andrian/ncfiles/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/'
save=True
# Funciones ############################################################################################################
def OpenDataSet(name, interp=False, lat_interp=None, lon_interp=None):
    from ENSO_IOD_Funciones import ChangeLons

    if name == 'pp_gpcc':
        # GPCC2018
        aux = xr.open_dataset(data_ref_dir + 'pp_gpcc.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        pp_gpcc = pp_gpcc.rename({'precip': 'var'})
        pp_gpcc = pp_gpcc.sel(time=slice('1982-01-01','2020-12-01'))

        return pp_gpcc

    elif name == 't_cru':
        # CRU
        aux = xr.open_dataset(data_ref_dir + 't_cru.nc')

        if interp:
            t_cru = aux.interp(lon=lon_interp, lat=lat_interp)

        t_cru = ChangeLons(aux)
        t_cru = t_cru.sel(lon=slice(270, 330), lat=slice(-60, 15),
                          time=slice('1982-01-01', '2020-12-31'))
        t_cru = t_cru.rename({'tmp': 'var'})
        t_cru = t_cru.drop('stn')

        return t_cru

def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         SA=False, dpi=100, save=True, step=1,
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
########################################################################################################################

variables = ('prec', 'tref')
variables_ref = ('pp_gpcc', 't_cru')
seasons = ('JJA', 'JAS', 'ASO', 'SON')
seasons_months = ([6,7,8],[7,8,9],[8,9,10],[9,10,11])
conjuntos = ([0],[0,1],[1,2,3],[4,5,6,7],[0,1,2,3,4,5,6,7])
conj_name= ('0','1','3','7', 'todo')
leads = (0,1,2,3,4,5,6,7)
colormap_mean = ['PuBuGn', 'YlOrRd']
scale_mean = [np.linspace(0,400,9), np.linspace(270,300,16)]
scale_sd = [np.linspace(0,70,9), np.linspace(0,1,11)]

"""
por conjuntos
0-1
1,2,3
4,5,6,7
y todos juntos (~Kumar)
"""

for v in variables:
    if v == variables[0]:
        v2 = variables_ref[0]
        v_name = 'precipitation'
        v_change_name = 'prec'
        factor=30
        sc_count=0
        factor_ref = 0
    else:
        v2 = variables_ref[1]
        v_name = 'temperature'
        v_change_name = 'tref'
        factor = 1
        sc_count = 1
        factor_ref=273

    s_count = 0
    for s in seasons:
        print(s)
        ms = seasons_months[s_count]

        data_mod = xr.open_dataset(models_dir + v + '_CFSv2_' + s + '_Lead_0.nc') # solo para la interpolacion de lo obs

        # Observed
        data_ref = OpenDataSet(v2, interp=True,
                               lat_interp=data_mod.lat.values,
                               lon_interp=data_mod.lon.values)
        data_ref += factor_ref
        data_ref = data_ref.sel(time=data_ref.time.dt.month.isin(ms))
        data_ref_mean_years = data_ref.groupby('time.year').mean()
        # SD
        data_ref_SD = data_ref_mean_years.std('year')
        # mean
        data_ref_mean = data_ref_mean_years.mean('year')

        Plot(data_ref_mean,
             levels=scale_mean[sc_count], cmap=colormap_mean[sc_count],
             dpi=200, step=1,
             name_fig=v + 'OBS_ref_set_' + s + '.nc',
             title='Obs - ' + s + '\n' + 'mean ' + v_name,
             save=save)

        Plot(data_ref_SD,
             levels=scale_sd[sc_count], cmap='Spectral_r',
             dpi=200, step=1,
             name_fig=v + 'OBS_SD_set_' + s + '.nc',
             title='Obs - ' + s + '\n' + 'SD ' + v_name,
             save=save)

        c_count = 0
        for c in conjuntos:
            c_name = conj_name[c_count]
            print(c_name)

            l_count = 0
            for l in c:
                # CFSv2
                data_mod = xr.open_dataset(models_dir + v + '_CFSv2_' + s + '_Lead_' + str(l) + '.nc').__mul__(factor)#mm/day
                data_mod = data_mod.rename({v_change_name: 'var'})
                data_mod_mean_years = data_mod.mean(['r'])
                # SD
                data_mod_SD = data_mod_mean_years.std('Year')
                # mean
                data_mod_mean = data_mod_mean_years.mean('Year')


                # mascara solo continental
                if v == variables[0]:
                    mask = xr.where(np.isnan(data_ref_mean), data_ref_mean, 1)
                data_mod_SD *= mask
                data_mod_mean *= mask

                if len(c) != 1:

                    data_mod_SD = data_mod_SD.expand_dims({'l': [l]})
                    # data_ref_SD = data_ref_SD.expand_dims({'l': [l]})
                    data_mod_mean = data_mod_mean.expand_dims({'l': [l]})
                    # data_ref_mean = data_ref_mean.expand_dims({'l': [l]})

                    if l_count == 0:
                        # data_ref_mean_f = data_ref_mean
                        data_mod_mean_f = data_mod_mean
                        # data_ref_SD_f = data_ref_SD
                        data_mod_SD_f = data_mod_SD
                    else:
                        # data_ref_mean_f = xr.merge([data_ref_mean_f, data_ref_mean])
                        data_mod_mean_f = xr.merge([data_mod_mean_f, data_mod_mean])
                        # data_ref_SD_f = xr.merge([data_ref_SD_f, data_ref_SD])
                        data_mod_SD_f = xr.merge([data_mod_SD_f, data_mod_SD])
                else:

                    # data_ref_mean_f = data_ref_mean
                    data_mod_mean_f = data_mod_mean
                    # data_ref_SD_f = data_ref_SD
                    data_mod_SD_f = data_mod_SD


                l_count += 1

            if len(c) != 1:
                # data_ref_mean_f = data_ref_mean_f.mean('l')
                data_mod_mean_f = data_mod_mean_f.mean('l')
                # data_ref_SD_f = data_ref_SD_f.mean('l')
                data_mod_SD_f = data_mod_SD_f.mean('l')

            # plot
            # MEAN
            Plot(data_mod_mean_f,
                 levels=scale_mean[sc_count], cmap=colormap_mean[sc_count],
                 dpi=200, step=1,
                 name_fig=v + 'MOD_mean_set_' + c_name + '_' + s + '.nc',
                 title='CFSv2 - ' + s + '\n' + 'mean ' + v_name + '\n' + 'Leads: ' + str(c),
                 save=save)

            # SD
            Plot(data_mod_SD_f,
                 levels=scale_sd[sc_count], cmap='Spectral_r',
                 dpi=200, step=1,
                 name_fig=v + 'MOD_SD_set_' + c_name + '_' + s + '.nc',
                 title='CFSv2 - ' + s + '\n' + 'SD ' + v_name + '\n' + 'Leads: ' + str(c),
                 save=save)

            c_count += 1
        s_count += 1
