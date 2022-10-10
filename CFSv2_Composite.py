import time
import xarray as xr
import numpy as np
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
"""
Primero, meterse en la climatología del modelo!!!
es necesario descargar las SST del modelo
verificar DMI en cada miembro y luego hacer la seleccion 

Quizas haya q alejarse del DMI de saji y ver uno mas de modelos.
"""
########################################################################################################################
models_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/'
data_ref_dir = '/datos/luciano.andrian/ncfiles/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/'
nc_date_dir = '/datos/luciano.andrian/ncfiles/nc_composites_dates/'
save=False
# Funciones ############################################################################################################
def CompositeSimpleCFSv2_r(original_data, index, l, name_case):
    """
    A diferencia del los demás, el promedio trimestral ya está hecho.

    """

    if len(index) != 0:
        comp_field = original_data.sel(Year=index)
        dims= comp_field['var'].shape
        r = 24
        y = dims[0]
        num_case = r * y
        if len(dims)>4:
            l = dims[1]
            num_case *=  l


        if len(comp_field.Year) != 0:
            if l > 1:
                comp_field = comp_field.mean(['r', name_case, 'l'], skipna=True)
            else:
                comp_field = comp_field.mean(['r', name_case], skipna=True)
        else:  # si sólo hay un año
            comp_field = comp_field.mean(['r'])

            comp_field = comp_field.drop_dims(['Year'])
            if l > 1:
                comp_field = comp_field.mean(['l'])


        return comp_field, num_case
    else:
        print(' len index = 0')

def CompositeSimple_SD_CFSv2_r(original_data, index, l, name_case):

    if len(index) != 0:
        comp_field = original_data.sel(Year=index)
        dims= comp_field['var'].shape
        r = 24
        y = dims[0]
        num_case = r * y
        if len(dims)>4:
            l = dims[1]
            num_case *=  l


        if len(comp_field.Year) != 0:
            if l > 1:
                comp_field = comp_field.std(['r', name_case, 'l'], skipna=True)
            else:
                comp_field = comp_field.std(['r', name_case], skipna=True)
        else:  # si sólo hay un año
            comp_field = comp_field.mean(['r'])

            comp_field = comp_field.drop_dims(['Year'])
            if l > 1:
                comp_field = comp_field.mean(['l'])


        return comp_field, num_case
    else:
        print(' len index = 0')

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

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI negative phase ',
              'DMI positive phase ',
              'DMI isolated positive phase ',
              'DMI isolated negative phase ',
              'ENSO positive phase ',
              'ENSO negative phase ',
              'ENSO isolated positive phase ',
              'ENSO isolated negative phase ']

scale = [np.linspace(-30,30,13), np.linspace(-1,1,13)]
scale_sd = [np.linspace(0,70,9), np.linspace(0,1,11)]

cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')

colormap= [cbar_pp, cbar_t]
colormap_sd = ['PuBuGn', 'YlOrRd']

cbar_snr = colors.ListedColormap(['#FFFFFF', '#FDE7A9', '#FEB77E', '#FB8761', '#CA3E72', '#A1307E',
                              '#782281',
                              '#4F127B', '#251255'])
cbar_snr.set_over('#251255')
cbar_snr.set_under('#FFFFFF')
cbar_snr.set_bad(color='white')

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

        data_mod = xr.open_dataset(models_dir + v + '_CFSv2_SON'+ '_Lead_0.nc').__mul__(factor) # solo para la interpolacion de lo obs

        # Observed
        mask = OpenDataSet('pp_gpcc', interp=True,
                               lat_interp=data_mod.lat.values,
                               lon_interp=data_mod.lon.values)
        mask = mask.mean('time')
        mask = xr.where(np.isnan(mask), mask, 1)

        c_count = 0
        for c in conjuntos:
            c_name = conj_name[c_count]
            print(c_name)

            l_count = 0
            for l in c:
                # CFSv2
                data_mod = xr.open_dataset(models_dir + v + '_CFSv2_' + s + '_Lead_' + str(l) + '.nc').__mul__(factor)#mm/day
                data_mod = data_mod.rename({v_change_name: 'var'})

                if len(c) != 1:
                    data_mod_aux = data_mod.expand_dims({'l': [l]})
                    if l_count == 0:
                        data_mod_i = data_mod_aux
                    else:
                        data_mod_l = xr.merge([data_mod_i, data_mod_aux])
                else:
                    data_mod_l = data_mod
                l_count += 1

            cs_count=0
            for cs in cases:
                aux = xr.open_dataset(nc_date_dir + 'Composite_1920_2020' + '_' + s + '.nc')
                neutro = aux.Neutral
                neutro = neutro.where(neutro>=1982,drop=True)
                try:
                    case = aux[cs]
                    case = case.where(case>=1982,drop=True)
                    aux.close()

                    neutro_comp, aux = CompositeSimpleCFSv2_r(data_mod_l, neutro, len(c), 'Neutral')
                    data_comp, num_case = CompositeSimpleCFSv2_r(data_mod_l, case, len(c), cs)

                    comp = data_comp - neutro_comp
                    # mascara solo continental
                    comp *= mask

                    Plot(comp,
                         levels=scale[sc_count], cmap=colormap[sc_count],
                         dpi=200, step=1,
                         name_fig=v + '_set_' + c_name + '_' + cs + '_' + s + '.nc',
                         title='CFSv2 - ' + s + '\n' + title_case[cs_count] + '\n' + ' ' + v_name + '\n' +
                               'Leads: ' + str(c) + ' - ' + 'Casos: ' + str(num_case),
                         save=save)

                    # SD
                    comp_SD, num_case = CompositeSimple_SD_CFSv2_r(data_mod_l, case, len(c), cs)
                    comp_SD *= mask

                    Plot(comp_SD,
                         levels=scale_sd[sc_count], cmap=colormap_sd[sc_count],
                         dpi=200, step=1,
                         name_fig=v + 'SD_set_' + c_name + '_' + cs + '_' + s + '.nc',
                         title='SD - CFSv2 - ' + s + '\n' + title_case[cs_count] + '\n' + ' ' + v_name + '\n' +
                               'Leads: ' + str(c) + ' - ' + 'Casos: ' + str(num_case),
                         save=save)

                    SNR = comp.__abs__()/comp_SD
                    Plot(SNR,
                         levels=np.linspace(0.1,1,10), cmap=cbar_snr,
                         dpi=200, step=1,
                         name_fig=v + 'SNR_set_' + c_name + '_' + cs + '_' + s + '.nc',
                         title='SNR - CFSv2 - ' + s + '\n' + title_case[cs_count] + '\n' + ' ' + v_name + '\n' +
                               'Leads: ' + str(c) + ' - ' + 'Casos: ' + str(num_case),
                         save=save)

                except:
                    print('No ANDAAAAA')

                cs_count +=1
            c_count += 1
        s_count += 1
