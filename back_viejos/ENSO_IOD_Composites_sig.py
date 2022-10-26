"""
Composites observado con filtro de tendencia corregida en las variables T, PP y las de ERA5
Por ahora sin SIGNIFICANCIA, HAY QUE VOLVER A HACER EL TEST MC
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
#from ENSO_IOD_Funciones import OpenDatasets
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import warnings
warnings.filterwarnings( "ignore" )
########################################################################################################################
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates/'
data_dir_t_pp = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/'
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/w_sig/'
save=True
dpi = 300

sig_dir = 'pikachu/datos/luciano.andrian/ncfiles/observado/nc_quantiles/'
sig_t_f = False
sig_c = False
########################################################################################################################

seasons = ("Full_Season", 'JJA', 'JAS', 'ASO', 'SON')
min_max_months = [[7,11], [6,8],[7,9],[8,10],[9,11]]

variables_t_p = ['t_cru_d_w_c', 'pp_gpcc_d_w_c']
variables_ERA5 = ['hgt200_mer_d_w', 'div200_mer_d_w', 'vp200_mer_d_w']

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']


scales = [np.linspace(-1, 1 ,17),  #t
          [-300, -250, -200, -150, -100, -50, -25, 0, 25, 50, 100, 150, 200, 250, 300],  # hgt
          np.linspace(-45, 45, 13),  #pp
          [-300,-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250,300],  #hgt
          np.linspace(-0.5e-5, 0.5e-5, 13),  #div200
          np.linspace(-4.5e6, 4.5e6, 13)]  #vp

          # np.linspace(-2e-10, 2e-10, 13), #RWS
          # np.linspace(-9, 9, 13), #u200
          # np.linspace(-3e-11, 3e-11, 13), #etay
          # np.linspace(-1, 1 ,17),  #t
          # np.linspace(-30, 30, 13)] #pp
          # np.linspace(-3, 3, 13),  #psl
          # np.linspace(-4.5e6, 4.5e6, 13),  #sf
          # np.linspace(-4.5e6, 4.5e6, 13),  #vp
          # np.linspace(-1, 1 ,17),  #t
          # np.linspace(-1, 1 ,17),  #t
          # np.linspace(-30, 30, 13)] #pp

SA = [True, True,False,False,False]
step = [1,1,1,1,4]
text = True

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
#
# v_name = ['HGT 200hPa', 'RWS', 'Div', 'U200'
#           'Stream Function', 'Potential Velocity',
#           'Temperature - Cru', 'Temperature - BEIC', 'Precipitation - GPCC']

v_name = [ 'Temperature - Cru',None, 'Precipitation - GPCC', 'HGT200hPa','Div200hPa', 'Potential Velocity']


cbar_r = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                                'white',
                                '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'])
cbar_r.set_under('#9B1C00')
cbar_r.set_over('#014A9B')
cbar_r.set_bad(color='white')

cbar = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar.set_over('#9B1C00')
cbar.set_under('#014A9B')
cbar.set_bad(color='white')

scale = [np.linspace(-30,30,13), np.linspace(-1,1,13)]
scale_sd = [np.linspace(0,70,9), np.linspace(0,1,11)]

cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')


cbar_t = colors.ListedColormap(['#9B1C00' ,'#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar_t.set_over('#691800')
cbar_t.set_under('#013774')
cbar_t.set_bad(color='white')

cbar_sst = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_sst.set_over('#9B1C00')
cbar_sst.set_under('#014A9B')
cbar_sst.set_bad(color='white')


cmap = [cbar_t, None, cbar_pp, cbar, cbar, cbar, cbar, 'BrBG']

## Functions ###########################################################################################################
def CompositeSimple(original_data, index, mmin, mmax):
    def is_months(month, mmin, mmax):
        return (month >= mmin) & (month <= mmax)

    if len(index) != 0:
        comp_field = original_data.sel(time=original_data.time.dt.year.isin([index]))
        comp_field = comp_field.sel(
            time=is_months(month=comp_field['time.month'], mmin=mmin, mmax=mmax))
        if len(comp_field.time) != 0:
            comp_field = comp_field.mean(['time'], skipna=True)
        else:  # si sólo hay un año
            comp_field = comp_field.drop_dims(['time'])

        return comp_field
    else:
        print(' len index = 0')

def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1, SA=False,
         name_fig='fig', title='title',
         two_variables=False, comp2=None, levels2=np.linspace(-1,1,13)):

    import matplotlib.pyplot as plt
    levels_contour = levels.copy()
    if isinstance(levels, np.ndarray):
        levels_contour = levels_contour[levels_contour != 0]
    else:
        levels_contour.remove(0)

    comp_var = comp['var']
    if SA:
        fig = plt.figure(figsize=(5, 6), dpi=dpi)
    else:
        fig = plt.figure(figsize=(8, 3), dpi=dpi)

    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    if SA:
        ax.set_extent([270,330, -60,20],crs=crs_latlon)
    else:
        ax.set_extent([30, 330, -80, 20], crs=crs_latlon)

    if two_variables:
        levels_contour2 = levels2.copy()
        comp_var2=comp2['var']
        if isinstance(levels2, np.ndarray):
            levels_contour2 = levels2[levels2 != 0]
        else:
            levels_contour2.remove(0)

        ax.contour(comp2.lon[::step], comp2.lat[::step], comp_var2[::step, ::step], linewidths=.8,
                   levels=levels_contour2, transform=crs_latlon, colors='black')

    else:
        ax.contour(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step], linewidths=.8,
                   levels=levels_contour, transform=crs_latlon, colors='black')


    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    if SA:
        ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
        ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    else:
        ax.set_xticks(np.arange(30, 330, 60), crs=crs_latlon)
        ax.set_yticks(np.arange(-80, 20, 20), crs=crs_latlon)

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

def CasesPlot(data, s, s_count, c_var, dpi=100,title='title', name_fig='name_fig',
              SA=True, save=False, two_variables=False, data2=None):

    # Las fechas se toman del periodo 1920-2020 basados en el DMI y N34 con ERSSTv5
    # Cuando se toman los periodos 1920-1949 y 1950_2020 las fechas que no pertencen
    # se excluyen de los composites en CompositeSimple()
    aux = xr.open_dataset(nc_date_dir + 'Composite_1920_2020' + '_' + s + '.nc')
    neutro = aux.Neutral

    try:
        case = aux[c]
        aux.close()
        mmonth = min_max_months[s_count]
        mmin = mmonth[0]
        mmax = mmonth[-1]

        neutro_comp = CompositeSimple(original_data=data, index=neutro, mmin=mmin, mmax=mmax)
        data_comp = CompositeSimple(original_data=data, index=case, mmin=mmin, mmax=mmax)

        comp = data_comp - neutro_comp

        if two_variables:
            neutro_comp2 = CompositeSimple(original_data=data2, index=neutro, mmin=mmin, mmax=mmax)
            data_comp2 = CompositeSimple(original_data=data2, index=case, mmin=mmin, mmax=mmax)

            comp2 = data_comp2 - neutro_comp2
        else:
            comp2 = None

        if len(comp) != 0:
            # if sig_c:
            #     sig = comp.where((comp < data_sig['var'][0]) | (comp > data_sig['var'][1]))
            #     try:
            #         sig = sig.where(np.isnan(sig['var']), 0)
            #     except:
            #         print('Sin valores significativos')
            # else:
            sig = None
            sig_t_f = False
            case = case[c].where(case[c], drop=True)

            Plot(comp=comp, levels=scales[c_var], cmap=cmap[c_var],
                 dpi=dpi, step=step[c_var], SA=SA,
                 name_fig=name_fig, title=title + ' - Events: '+ str(len(case)), save=save,
                 two_variables=two_variables, comp2=comp2, levels2=scales[c_var+1])
    except:
        print('No hay eventos ' + c)


def PlotSST(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt
    levels_contour = levels.copy()
    if isinstance(levels, np.ndarray):
        levels_contour = levels_contour[levels_contour != 0]
    else:
        levels_contour.remove(0)

    comp_var = comp['var']

    fig = plt.figure(figsize=(7, 2), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()

    ax.set_extent([50, 270, -20, 20], crs=crs_latlon)

    ax.contour(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step], linewidths=.8,
                     levels=levels_contour, transform=crs_latlon, colors='black')

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)

    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')

    ax.set_xticks(np.arange(50, 270, 30), crs=crs_latlon)
    ax.set_yticks(np.arange(-20, 20, 10), crs=crs_latlon)

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
v_count=0
for v in variables_t_p:
    print('Variable: ' + v)
    p = '1950-2020'

    print(p)
    data = xr.open_dataset(data_dir_t_pp + v + '_' + p + '_0.25.nc')
    data2 = xr.open_dataset(data_dir_era5 + variables_ERA5[0] + '.nc')
    data2 = data2.sel(lon=slice(270, 330), lat=slice(15, -60))
    data2 = data2.interp(lon=data.lon.values, lat=data.lat.values)
    tw = True

    c_cases = 0
    for c in cases:
        print(c)
        s_count = 0
        for s in seasons:
            print(s)
            CasesPlot(data, s, s_count, c_var=v_count, dpi=dpi, two_variables=tw, data2=data2,
                      title=v_name[v_count] + '\n' + title_case[c_cases] + '\n' + s + ' ' + p,
                      name_fig=v + '_' + cases[c_cases] + '_' + s + '_' + p, save=save, SA=True)

            s_count += 1
        c_cases += 1
    v_count += 2


tw=[None, None, None, False,False]
v_count=3
for v in variables_ERA5:
    print('Variable: ' + v)
    data = xr.open_dataset(data_dir_era5 + v + '.nc')
    if tw[v_count]:
        data2 = xr.open_dataset(data_dir_era5 + variables_ERA5[v_count + -2] + '.nc')
    else:
        data2 = None

    c_cases = 0
    for c in cases:
        print(c)
        s_count = 0
        for s in seasons:
            CasesPlot(data, s, s_count, c_var=v_count, dpi=dpi, two_variables=tw[v_count], data2=data2,
                      title=v_name[v_count] + '\n' + title_case[c_cases] + '\n' + s,
                      name_fig=v + '_' + cases[c_cases] + '_' + s, save=save, SA=False)

            s_count += 1
        c_cases += 1

    v_count += 1
    if v_count==5:
        break


# SST ##################################################################################################################
"""
Se filtra la tendencia lineal similar a lo que se hace en el DMI
En el N34 se muve la climatologia cada 5 años. 
Puede haber algunos problemas.
"""
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

data = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
data = data.sel(lon=slice(50, 270),lat=slice(20, -20))
data = data.rename({'sst':'var'})
data = data.drop('time_bnds')
data_20_20 = data.sel(time=slice('1920-01-01', '2020-12-01'))
data_50_20 = data.sel(time=slice('1950-01-01', '2020-12-01'))

data_20_20 = Detrend(data_20_20, 'time')
data_50_20 = Detrend(data_50_20, 'time')


c_cases = 0
for c in cases:
    print(c)
    s_count = 0
    for s in seasons:

        aux = xr.open_dataset(nc_date_dir + 'Composite_1920_2020' + '_' + s + '.nc')
        neutro = aux.Neutral

        case = aux[c]
        aux.close()
        mmonth = min_max_months[s_count]
        mmin = mmonth[0]
        mmax = mmonth[-1]

        neutro_comp = CompositeSimple(original_data=data, index=neutro, mmin=mmin, mmax=mmax)
        data_comp = CompositeSimple(original_data=data, index=case, mmin=mmin, mmax=mmax)

        comp = data_comp - neutro_comp

        if len(comp) != 0:
            # if sig_c:
            #     sig = comp.where((comp < data_sig['var'][0]) | (comp > data_sig['var'][1]))
            #     try:
            #         sig = sig.where(np.isnan(sig['var']), 0)
            #     except:
            #         print('Sin valores significativos')
            # else:
            sig = None
            sig_t_f = False
            case = case[c].where(case[c], drop=True)

            PlotSST(comp=comp, levels=np.linspace(-1.2,1.2,13), cmap=cbar_sst,
                 dpi=dpi, step=1,
                 name_fig='sst_' + cases[c_cases] + '_' + s,
                 title='SST - ' + title_case[c_cases] + '\n' + s + ' - Events: '+ str(len(case)),
                 save=save)

        s_count += 1
    c_cases += 1
########################################################################################################################

