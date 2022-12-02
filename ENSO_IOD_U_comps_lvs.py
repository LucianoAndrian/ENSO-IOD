import xarray as xr
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from cartopy.mpl.ticker import LatitudeFormatter
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings("ignore")
########################################################################################################################
dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/'
out_dir = ['/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/no_sig/U_zonal/',
           '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/no_sig/DMIbase/U_zonal/']

nc_date_dir1 = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates/' #fechas
nc_date_dir2 = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/' #fechas

dir_dates = [nc_date_dir1, nc_date_dir2]

set = False
save = True
step=1

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
    aux2 = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux2.var_polyfit_coefficients[0])
    except:
        trend = xr.polyval(xrda[dim], aux2.polyfit_coefficients[0])
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

def CaseComp(data, s, mmonth, c, two_variables=False, data2=None, nc_date_dir=None, clim=False,
             zonalA=False, ZA_type='vertical'):
    """
    Las fechas se toman del periodo 1920-2020 basados en el DMI y N34 con ERSSTv5
    Cuando se toman los periodos 1920-1949 y 1950_2020 las fechas que no pertencen
    se excluyen de los composites en CompositeSimple()
    """
    mmin = mmonth[0]
    mmax = mmonth[-1]

    if zonalA:
        two_variables=True
        if data2 is None:
            return print('Error en data2')
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

        if zonalA:
            if ZA_type == 'vertical':
                comp = data_comp.mean('lon') - neutro_comp2.mean('lon')
            else:
                comp = data_comp - neutro_comp2.mean('lon')
    except:
        print('Error en ' + s + c)

    if two_variables:
        if zonalA:
            return comp, case_num
        else:
            return comp, case_num, comp2
    else:
        if clim:
            return data_comp, case_num
        else:
            return comp, case_num

#Plot

########################################################################################################################
if set:
    try:
        u_lon = xr.open_dataset(dir + 'ERA5_U_lvs_xrmer.nc')
    except:
        u_50_78 = xr.open_dataset(dir + 'ERA5_U200_lvs_50_78.nc')
        u_79_20 = xr.open_dataset(dir + 'ERA5_U200_lvs_79_20.nc')
        print('ALGO TA MAAAAL!!!!')
        u = xr.merge([u_50_78, u_79_20])
        u.to_netcdf(dir + 'ERA5_U_lvs_xrmer.nc')
    u_lon = u_lon.rename({'longitude':'lon', 'latitude':'lat', 'u':'var'})
    u_lon = Weights(u_lon)
    u_lon = Detrend(u_lon, 'time')
    u_lon.to_netcdf('/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/ERA5_U_mer_d_w.nc')
########################################################################################################################
seasons = ('JJA', 'SON')
min_max_months = [[6,8], [9,11]]
mm = [7, 10]
cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_un_pos',
         'DMI_un_neg', 'N34_un_pos', 'N34_un_neg']

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

regiones = [[40,130], [130,290]]
regiones_name = ['Indian Ocean', 'Pacific Ocean']
mm = [7,10]
########################################################################################################################
data = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/ERA5_U_mer_d_w.nc')
data = data.sortby('level', ascending=False) # estan TODOS los niveles

data_dir_count=1
d_d = dir_dates[1]
c_count=4
c=cases[4]
s_count=1
s = 'SON'

lv = 850
date_dir_count = 0
for d_d in dir_dates:

    c_count = 0
    for c in cases:
        s_count = 0
        for s in seasons:

            fig = plt.figure(figsize=(8,4), dpi=400)
            gs = fig.add_gridspec(2, 2)

            # campos de anomalias zonales -----------------------------------------------------------------------------#
            lv_count = 0
            for lv in [850, 200]:
                ax = fig.add_subplot(gs[lv_count, 0],
                                     projection=ccrs.PlateCarree(central_longitude=180))
                crs_latlon = ccrs.PlateCarree()
                ax.set_extent([30, 340, -80, 10], crs=crs_latlon)

                #-- composicion de la anomalia zonal de U en lv --#
                aux_comp = data.sel(level=lv)
                comp, num_case = CaseComp(data=aux_comp, s=s, mmonth=min_max_months[s_count],
                                          c=c, two_variables=True, nc_date_dir=d_d,
                                          zonalA=True, ZA_type='campo', data2=aux_comp)


                #-- media climatologica de U en lv --#
                aux_clim = aux_comp.rolling(time=3, center=True).mean()
                aux_clim = aux_clim.sel(time=aux_clim.time.dt.month.isin(mm[s_count])).mean('time')

                # Plot ----------------------------------------------------------------------------------#
                ax.contour(aux_clim.lon[::step], aux_clim.lat[::step], aux_clim['var'][::step, ::step],
                           transform=crs_latlon, colors=['#5A5A5A', '#313131', 'k'],
                           levels=[20, 30, 40], linewidths=1.5)
                im = ax.contourf(comp.lon[::step], comp.lat[::step], comp['var'][::step, ::step],
                                 levels=[-12, -10, -8, -6, -4, -2, -1, 0, 1, 2, 4, 6, 8, 10, 12],
                                 transform=crs_latlon, cmap=cbar_t, extend='both')
                cb = plt.colorbar(im, fraction=0.018, pad=0.035, shrink=0.8)
                cb.ax.tick_params(labelsize=5)
                ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey', edgecolor='#4B4B4B')
                ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
                ax.coastlines(color='#4B4B4B', linestyle='-', alpha=1)
                ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
                ax.set_xticks(np.arange(30, 340, 25), crs=crs_latlon)
                ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
                lon_formatter = LongitudeFormatter(zero_direction_label=True)
                lat_formatter = LatitudeFormatter()
                ax.xaxis.set_major_formatter(lon_formatter)
                ax.yaxis.set_major_formatter(lat_formatter)
                ax.tick_params(labelsize=5)
                ax.set_title('U* y U clim. [cont.] - ' + str(lv) + 'hPa',  fontsize=10)
                plt.tight_layout()

                lv_count += 1
            #----------------------------------------------------------------------------------------------------------#
            # corte vertical de U en la region r  ---------------------------------------------------------------------#
            r_count = 0
            for r in regiones:
                ax = fig.add_subplot(gs[r_count, 1])

                aux_comp0 = data.sel(lat=slice(10, -80)) # sel en lat.
                aux_comp = aux_comp0.sel(lon=slice(r[0], r[1])) # selecciona la region

                #-- climatologia de U en r --#
                aux_clim = aux_comp.rolling(time=3, center=True).mean()
                aux_clim = aux_clim.sel(time=aux_clim.time.dt.month.isin(mm[s_count])).mean('time')
                aux_clim = aux_clim.sel(lat=slice(10, -80),lon=slice(r[0], r[1]))
                aux_clim = aux_clim.mean(['lon'])

                # -- composicion de la anomalia zonal de U en r --#
                comp, num_case = CaseComp(data=aux_comp, s=s, mmonth=min_max_months[s_count],
                                           c=c, two_variables=True, nc_date_dir=d_d,
                                           zonalA=True, ZA_type='vertical', data2=aux_comp0)

                # Plot ------------------------------------------------------------------------------------------------#
                ax.contour(aux_clim.lat[::step], aux_clim.level[::step], aux_clim['var'][::step, ::step], linewidths=1,
                           levels=np.arange(-10, 35, 5), colors='k')
                im = ax.contourf(comp.lat[::step], comp.level[::step], comp['var'][::step, ::step],
                                 levels=[-8, -6, -4, -2, -1, -.5, 0, .5, 1, 2, 4, 6, 8],
                                 cmap=cbar_t, extend='both')
                cb = plt.colorbar(im, fraction=0.042, pad=0.035, shrink=0.8)
                cb.ax.tick_params(labelsize=6)
                ax.set_xticks(np.arange(-80, 0, 10))
                lat_formatter = LatitudeFormatter()
                ax.xaxis.set_major_formatter(lat_formatter)
                ax.tick_params(labelsize=7)
                plt.yscale('log')
                ax.set_ylabel("Pressure [hPa]", fontsize=8)
                ax.set_yscale('log')
                ax.set_ylim(10. * np.ceil(comp.level.values.max() / 10.), 100)
                subs = [1, 2, 5]
                if comp.level.values.max() / 100 < 100.:
                    subs = [1, 2, 3, 4, 5, 6, 7, 8, 9]
                y1loc = matplotlib.ticker.LogLocator(base=10., subs=subs)
                ax.yaxis.set_major_locator(y1loc)
                fmt = matplotlib.ticker.FormatStrFormatter("%g")
                ax.yaxis.set_major_formatter(fmt)
                ax.set_title('[U*] - ' + regiones_name[r_count], fontsize=10)
                plt.tight_layout()

                r_count += 1

            #-- titulo y guardado del plot --#
            fig.suptitle(c + 'number cases:' + str(num_case))
            if save:
                plt.savefig(out_dir[date_dir_count] + 'U_zonal_' + c + '_' + s + '.jpg', d=400)
            else:
                plt.show()

            s_count += 1
        c_count += 1
    date_dir_count += 1
########################################################################################################################