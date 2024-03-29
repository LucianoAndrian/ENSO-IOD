import time
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
from ENSO_IOD_Funciones import OpenDatasets
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
########################################################################################################################
nc_date_dir = '/datos/luciano.andrian/ncfiles/nc_composites_dates/'
#data_dir = '/datos/luciano.andrian/ncfiles/'
data_dir1 = '/datos/luciano.andrian/ncfiles/ERA20C/'
data_dir2 = '/datos/luciano.andrian/ncfiles/ERA5/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/w_sig/'
sig_dir = '/datos/luciano.andrian/ncfiles/nc_quantiles/'

periodos = [['50','20']]
seasons = ("Full_Season", 'JJA', 'ASO', 'SON')
min_max_months = [[7,11], [6,8],[8,10],[9,11]]

variables = ['u']

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

SA = [False,False,False,False,False, True, True]
step = [1,4,4,1,4,1,1]
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

v_name = ['U200']

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

cmap = [cbar, cbar, cbar, cbar, cbar, cbar, 'BrBG']
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

def Plot(comp, comp_sig, significance=False,
         levels = np.linspace(-1,1,11), cmap='RdBu',
         two_variables=False, comp2=None, comp2_sig=None,
         significance2=False, levels2 = np.linspace(-1,1,11),
         number_events=None, text=True,SA=False, dpi=100,
         contour0=False,save=False, step=1,name_fig='fig', title='title'):

    from numpy import ma
    import matplotlib.pyplot as plt

    comp_var  = comp['var']
    number_events = len(number_events)

    if SA:
        fig = plt.figure(figsize=(5, 6), dpi=dpi)
    else:
        fig = plt.figure(figsize=(7, 3.5), dpi=dpi)

    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    if SA:
        ax.set_extent([270,330, -60,20],crs=crs_latlon)
    else:
        ax.set_extent([0, 359, -80, 40], crs=crs_latlon)


    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step,::step],
                     levels=levels,transform=crs_latlon, cmap=cmap, extend='both')

    if contour0:
        values = ax.contour(comp.lon, comp.lat, comp_var, levels=0,
                            transform=crs_latlon, colors='magenta', linewidths=1)
        ax.clabel(values, inline=1.5, fontsize=5, fmt='%1.1f')

    if significance:
        comp_sig_var = comp_sig['var']
        values = ax.contourf(comp_sig.lon, comp_sig.lat, comp_sig_var, levels=levels,
                            transform=crs_latlon, colors='k', hatches=["", "..."], alpha=0)

    if two_variables:
        print('comp2 as contour')
        comp2_var = comp2['var']
        values = ax.contour(comp2.lon, comp2.lat, comp2_var, levels=levels2,
                            transform=crs_latlon, colors='forestgreen', linewidths=1)
        ax.clabel(values, inline=1, fontsize=5, fmt='%1.1f')

        if significance2:
            comp2_sig_var = comp2_sig['var']
            values = ax.scatter(comp2_sig.lon, comp2_sig.lat, comp2_sig_var, levels=[-1, 1],
                                transform=crs_latlon, colors='magenta', linewidths=1)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
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
    if text:
        plt.figtext(0.5, 0.01, 'Number of events: ' + str(number_events), ha="center", fontsize=10,
                bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5})
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()
########################################################################################################################
c_var=0
v = 'HGT200'
sig_t_f=True
for v in variables:
    print('Variable: ' + v)
    if v != 'HGT200':
        sig_c = False
    else:
        sig_c = True
    for i in periodos:
        print('Período: ' + i[0] + '-' + i[1])
        if i[0] == '20':
            ruta = data_dir1
            name = 'ERA20'
        elif i[0] == '50':
            ruta = data_dir2
            name = 'ERA5'
        else:
            print('TA MAL!')

        # print('Período: ' + i + '- 2020')
        # print('Open ' + v + '.nc')
        # if (v == 't_cru') | (v == 't_BEIC') | (v == 'pp_gpcc'):
        #     print('Using OpenDatasets')
        #     data = OpenDatasets(name=v, interp=True)  # Mc se hizo en 1x1
        # else:
        #     data = xr.open_dataset(data_dir + v + '.nc')

        print('Open ' + ruta + name + '_' + v + '_' + i[0] + '-' + i[1] + '.nc')
        if v != 'HGT200':
            data = xr.open_dataset(ruta + v + '.nc')
            if v == 'RWS':
                data = data.rename({'__xarray_dataarray_variable__' : 'var'})
            elif v == 'div200_from_W':
                data = data.rename({'divergence' : 'var'})
            elif v == 'etay':
                data = data.rename({'meridional_gradient_of_absolute_vorticity':'var'})
        else:
            data = xr.open_dataset(ruta + name + '_' + v + '_' + i[0] + '-' + i[1] + '.nc')

        if v == 'psl':
            print('to hPa')
            data = data.__mul__(1 / 100)

        c_cases = 0
        for c in cases:
            print(c)

            count = 0
            for s in seasons:

                if sig_c:
                    # Significance # En caso de correr en conjunto con comp_sig_MC.py
                    for m in [1, 2, 3, 4, 5, 6, 8, 9, 10]:
                        try:
                            if m != 1:
                                print('Waiting..')
                                time.sleep(300)  # 5min
                            else:
                                print('first try')
                            if name == 'ERA20':
                                data_sig = xr.open_dataset(
                                    sig_dir + v + '_' + c + '_19' + i[0] + '_19' + i[1] + '_' + s + '.nc')
                            elif name == 'ERA5':
                                data_sig = xr.open_dataset(
                                    sig_dir + v + '_' + c + '_19' + i[0] + '_20' + i[1] + '_' + s + '.nc')

                            break
                        except:
                            try:
                                s2 = seasons[count + 1]
                                if name == 'ERA20':
                                    data_sig = xr.open_dataset(
                                        sig_dir + v + '_' + c + '_19' + i[0] + '_19' + i[1] + '_' + s2 + '.nc')
                                elif name == 'ERA5':
                                    data_sig = xr.open_dataset(
                                        sig_dir + v + '_' + c + '_19' + i[0] + '_20' + i[1] + '_' + s2 + '.nc')
                                print(v + '_' + c + '_' + s + 'No existe')
                                s = s2
                                print('pasando a siguiente archivo')
                                break
                            except:
                                print(v + '_' + c + '_' + s + 'aún no disponible')

                # Las fechas se toman del periodo 1920-2020 basados en el DMI y N34 con ERSSTv5
                # Cuando se toman los periodos 1920-1949 y 1950_2020 las fechas que no pertencen
                # se excluyen de los composites en CompositeSimple()
                aux = xr.open_dataset(nc_date_dir + 'Composite_1920_2020' + '_' + s + '.nc')

                neutro = aux.Neutral
                try:
                    case = aux[c]
                    aux.close()
                    mmonth = min_max_months[count]

                    mmin = mmonth[0]
                    mmax = mmonth[-1]

                    neutro_comp = CompositeSimple(original_data=data, index=neutro,
                                                  mmin=mmin, mmax=mmax)
                    data_comp = CompositeSimple(original_data=data, index=case,
                                                mmin=mmin, mmax=mmax)
                    comp = data_comp - neutro_comp

                    if len(comp) != 0:
                        if sig_c:
                            sig = comp.where((comp < data_sig['var'][0]) | (comp > data_sig['var'][1]))

                            try:
                                sig = sig.where(np.isnan(sig['var']), 0)
                            except:
                                print('Sin valores significativos')

                        else:
                            sig=None
                            sig_t_f=False

                        if name == 'ERA20':
                            title = v_name[c_var] + '\n' + title_case[c_cases] + '\n' + s + '  ' + '1920-1949_ERA-20C'
                            case = case[c].where(case[c] < 1950, drop=True)
                        elif name == 'ERA5':
                            title = v_name[c_var] + '\n' + title_case[c_cases] + '\n' + s + '  ' + '1950-2020_ERA5'
                            case = case[c].where(case[c] > 1949, drop=True)

                        print(title_case[c_cases])
                        Plot(comp=comp, comp_sig=sig, significance=sig_t_f, cmap=cbar,
                             levels=np.linspace(-9, 9, 13), number_events=case,
                             SA=False, step=1, text=text, dpi=200,
                             title=title, two_variables=True, comp2=data_comp, levels2=[30,40],
                             name_fig=v + '_Comp_' + c + '_' + s + '_' + i[0] + '-' + i[1],
                             save=True)

                except:
                    print('No hay eventos ' + c)

                count += 1
            c_cases += 1
    c_var += 1
