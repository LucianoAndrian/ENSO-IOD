import xarray as xr
import numpy as np
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
from ENSO_IOD_Funciones import OpenDatasets
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import warnings
warnings.filterwarnings("ignore")
save=True
########################################################################################################################
nc_date_dir = '/datos/luciano.andrian/ncfiles/nc_composites_dates/'
sig_dir = '/datos/luciano.andrian/ncfiles/nc_quantiles/'
data_dir = '/datos/luciano.andrian/ncfiles/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/w_sig/T_PP/'


seasons = ("Full_Season", 'JJA', 'JAS', 'ASO', 'SON')
min_max_months = [[7,11], [6,8],[7,9],[8,10],[9,11]]
variables = ['t_cru', 'pp_gpcc']

scales = [np.linspace(-1, 1 ,17),  #t
          np.linspace(-30, 30, 13)]
SA=True
step=1
text=True


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

v_name = [ 'T-Obs', 'PP-Obs']


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


cmap = [cbar_t, cbar_pp]

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
                            transform=crs_latlon, colors='darkgray', linewidths=1)
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
for v in variables:
    print(v)
    print('Using OpenDatasets')
    data = OpenDatasets(name=v, interp=True)  # Mc se hizo en 1x1

    c_cases = 0
    for c in cases:
        s_count = 0
        for s in seasons:
            print(s)
            if s != 'JAS':
                sig_t_f = True
            else:
                sig_t_f = False

            # fechas
            aux = xr.open_dataset(nc_date_dir + 'Composite_1920_2020' + '_' + s + '.nc')
            neutro = aux.Neutral

            try:
                case = aux[c]
                aux.close()
                mmonth = min_max_months[s_count]

                mmin = mmonth[0]
                mmax = mmonth[-1]

                neutro_comp = CompositeSimple(original_data=data, index=neutro,
                                              mmin=mmin, mmax=mmax)
                data_comp = CompositeSimple(original_data=data, index=case,
                                            mmin=mmin, mmax=mmax)
                comp = data_comp - neutro_comp

                if len(comp) != 0:
                    if sig_t_f:
                        data_sig = xr.open_dataset(sig_dir + v + '_' + c + '_1920_2020_' + s + '.nc')
                        sig = comp.where((comp < data_sig['var'][0]) | (comp > data_sig['var'][1]))

                        try:
                            sig = sig.where(np.isnan(sig['var']), 0)
                        except:
                            print('Sin valores significativos')
                            sig = None
                            sig_t_f = False

                    title = v_name[c_var] + '\n' + title_case[c_cases] + '\n' + s + '  ' + '1920-1920'
                    case = case[c]

                    print(title_case[c_cases])
                    Plot(comp=comp, comp_sig=sig, significance=sig_t_f, cmap=cmap[c_var],
                         levels=scales[c_var], number_events=case,
                         SA=SA, step=step, text=text, dpi=200,
                         title=title,
                         name_fig=v + '_Comp_' + c + '_' + s,
                         save=save)
            except:
                print(str(c) + '_' + s + ' no data')

            s_count += 1
        c_cases += 1
    c_var += 1