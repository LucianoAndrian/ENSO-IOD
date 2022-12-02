import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings("ignore")
########################################################################################################################
dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/'
out_dir = ['/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/no_sig/',
           '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/no_sig/DMIbase/']

nc_date_dir1 = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates/' #fechas
nc_date_dir2 = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/' #fechas

dir_dates = [nc_date_dir1, nc_date_dir2]

set = True
save = True
dpi = 200
#lon_cut=60
########################################################################################################################
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

def DetrendClim(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients[0])
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients[0])
    dt = xrda - trend
    return dt


def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180), (len(data.lon), 1)))
    data_w = data * weights
    return data_w


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


def CaseComp(data, s, mmonth, c, two_variables=False, data2=None, nc_date_dir=None, clim=False):
    """
    Las fechas se toman del periodo 1920-2020 basados en el DMI y N34 con ERSSTv5
    Cuando se toman los periodos 1920-1949 y 1950_2020 las fechas que no pertencen
    se excluyen de los composites en CompositeSimple()
    """
    mmin = mmonth[0]
    mmax = mmonth[-1]

    aux = xr.open_dataset(nc_date_dir + '1920_2020' + '_' + s + '.nc')
    neutro = aux.Neutral

    try:
        case = aux[c]
        case = case.where(case >= 1950)
        aux.close()

        case_num = len(case.values[np.where(~np.isnan(case.values))])

        neutro_comp = CompositeSimple(original_data=data, index=neutro, mmin=mmin, mmax=mmax)
        data_comp = CompositeSimple(original_data=data, index=case, mmin=mmin, mmax=mmax)

        comp = data_comp - neutro_comp

        if two_variables:
            neutro_comp2 = CompositeSimple(original_data=data2, index=neutro, mmin=mmin, mmax=mmax)
            data_comp2 = CompositeSimple(original_data=data2, index=case, mmin=mmin, mmax=mmax)

            comp2 = data_comp2 - neutro_comp2
        else:
            comp2 = None
    except:
        print('Error en ' + s + c)

    if two_variables:
        return comp, case_num, comp2
    else:
        if clim:
            return data_comp, case_num
        else:
            return comp, case_num, data_comp, neutro_comp

def Plot(comp, comp_var, levels = np.linspace(-1,1,11),
         cmap='RdBu', dpi=100, save=True, step=1,
         name_fig='fig', title='title', color_map='#4B4B4B',
         comp2=None, comp_var2=None, levels2=[30,40], colors2='k', lwd=1,
         out_dir='~/'):

    levels_contour = levels.copy()
    levels_contour.remove(0)
    #comp_var = comp['var']
    fig = plt.figure(figsize=(8, 3), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([30,340, -80,20], crs_latlon)

    # ax.contour(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step], linewidths=.8,
    #                  levels=levels_contour, transform=crs_latlon, colors='black')
    ax.contour(comp2.lon[::step], comp2.lat[::step], comp_var2[::step, ::step], linewidths=lwd,
               levels=levels2, transform=crs_latlon, colors=colors2)
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

########################################################################################################################
if set:
    try:
        u200 = xr.open_dataset(dir + 'ERA5_U_lvs_xrmer.nc').sel(level=200)
    except:
        u_50_78 = xr.open_dataset(dir + 'ERA5_U200_lvs_50_78.nc')
        u_79_20 = xr.open_dataset(dir + 'ERA5_U200_lvs_79_20.nc')
        print('ALGO TA MAAAAL!!!!')
        u = xr.merge([u_50_78, u_79_20])
        u.to_netcdf(dir + 'ERA5_U_lvs_xrmer.nc')

    u200 = u200.rename({'longitude':'lon', 'latitude':'lat', 'u':'var'})


    u200_clim = DetrendClim(u200, 'time')
    u200_clim = u200_clim.rolling(time=3, center=True).mean()
    u200 = Detrend(u200, 'time')
    #u200.to_netcdf('/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/ERA5_U_mer_d_w.nc')

seasons = ('JJA', 'SON')
min_max_months = [[6,8], [9,11]]
mm = [7, 10]
cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos',
         'DMI_un_neg','N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI negative phase ',
              'DMI positive phase ',
              'DMI pure positive phase ',
              'DMI pure negative phase ',
              'ENSO positive phase ',
              'ENSO negative phase ',
              'ENSO pure positive phase ',
              'ENSO pure negative phase ']

# colorbars
cbar_t = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar_t.set_over('#691800')
cbar_t.set_under('#013774')
cbar_t.set_bad(color='white')


data = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/ERA5_U_mer_d_w.nc')
data = data.sortby('level', ascending=False)
data = data.sel(level=850)
u200 = data
#----------------------------------------------------------------------------------------------------------------------#

date_dir_count = 0
for nc_date_dir in dir_dates:
    c_count = 0
    for c in cases:
        s_count = 0
        for s in seasons:
            #anom
            comp1, num_case, data_comp, neutro_comp = CaseComp(u200.sel(lat=slice(20, -80)), s, mmonth=min_max_months[s_count],
                                       c=c, two_variables=False, data2=None, nc_date_dir=nc_date_dir)

            comp2 = u200_clim.sel(lat=slice(20, -80), time=u200_clim.time.dt.month.isin(mm[s_count]))
            comp2 = comp2.mean(['time'])

            if date_dir_count == 1:
                add_to_name_fig = 'no_SSTanom'
            else:
                add_to_name_fig = ''

            Plot(comp=comp1, comp_var=comp1['var'], cmap=cbar_t,
                 levels=[-12, -10, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 10, 12],
                  title="U'[somb.] y Umed 20 y 30ms [Cont.] - 200hPa " + add_to_name_fig + ' ' +
                        '\n' + title_case[c_count] + '\n' + s + ' - Events: ' + str(num_case),
                  name_fig="U200_" + s + '_' + cases[c_count] + '_' + '_mer_d_w' + add_to_name_fig,
                  dpi=dpi, save=save, out_dir=out_dir[date_dir_count],
                 comp2=comp2, comp_var2=comp2['var'],
                 colors2=['#5A5A5A', '#313131', 'k'], lwd=1.5, levels2=[20, 30, 40])

            s_count += 1
        c_count += 1
    date_dir_count += 1

#----------------------------------------------------------------------------------------------------------------------#
########################################################################################################################

plt.plot(data_comp.sel(lon=slice(130,290)).mean('lon')['var'], label='El Niño - Pacifico');plt.plot(neutro_comp.mean('lon')['var'], label='neutro', color='k');plt.xticks(np.arange(0, len(data_comp.lat.values),20),np.linspace(data_comp.lat[0], data_comp.lat[-1],11));plt.legend();plt.show()
