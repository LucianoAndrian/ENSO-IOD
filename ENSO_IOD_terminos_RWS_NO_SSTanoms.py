import xarray as xr
import numpy as np
from matplotlib import colors
import warnings
import matplotlib
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
########################################################################################################################
dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/RWS_terms/no_sstanoms/'
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/'

save=True
dpi=200
# Functions ############################################################################################################
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
            print('1 year')
            comp_field = comp_field.drop_dims(['time'])

        return comp_field
    else:
        print(' len index = 0')


def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1, color_map='#4B4B4B',
         name_fig='fig', title='title', contour0=True,):

    import cartopy.feature
    from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
    import cartopy.crs as ccrs
    import numpy as np
    import matplotlib.pyplot as plt

    comp_var = comp['var']

    fig = plt.figure(figsize=(8, 3), dpi=dpi)

    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([30, 330, -80, 20], crs=crs_latlon)


    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step,::step],
                     levels=levels,transform=crs_latlon, cmap=cmap, extend='both')
    if contour0:
        ax.contour(comp.lon, comp.lat, comp_var, levels=0,
                   transform=crs_latlon, colors='green', linewidths=1)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE)
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
seasons = ['JJA', 'SON']
min_max_months = [[6, 8], [9, 11]]

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

title_case = ['DMI-ENSO simultaneous positive phase ','DMI-ENSO simultaneous negative phase ',
              'DMI negative phase ', 'DMI positive phase ',
              'DMI pure positive phase ', 'DMI pure negative phase ',
              'ENSO positive phase ', 'ENSO negative phase ',
              'ENSO pure positive phase ', 'ENSO pure negative phase ']

variables = ['RWS', 'u', 'etay', 'div200_from_W']

scales = [np.linspace(-2e-10, 2e-10, 13),
          np.linspace(-12, 12, 13),
          np.linspace(-3e-11, 3e-11, 13),
          np.linspace(-0.5e-5, 0.5e-5, 13)]

change_names = ['__xarray_dataarray_variable__', None,
                'meridional_gradient_of_absolute_vorticity',
                'divergence']

#----------------------------------------------------------------------------------------------------------------------#
cbar_r = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'])
cbar_r.set_under('#641B00')
cbar_r.set_over('#012A52')
cbar_r.set_bad(color='white')

cbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar.set_over('#641B00')
cbar.set_under('#012A52')
cbar.set_bad(color='white')
########################################################################################################################
#term1
aux = xr.open_dataset(dir + 'eta.nc')
aux = aux.rename({'absolute_vorticity':'var'})
aux2 = xr.open_dataset(dir + 'div_from_w.nc')
aux2 = aux2.rename({'divergence':'var'})

term1= -1*aux*aux2


#term2
aux = xr.open_dataset(dir + 'uchi.nc')
aux = aux.rename({'u_chi':'var'})
aux2 = xr.open_dataset(dir + 'vchi.nc')
aux2 = aux2.rename({'v_chi':'var'})

aux3 = xr.open_dataset(dir + 'etax.nc')
aux3 = aux3.rename({'zonal_gradient_of_absolute_vorticity':'var'})
aux4 = xr.open_dataset(dir + 'etay.nc')
aux4 = aux4.rename({'meridional_gradient_of_absolute_vorticity':'var'})
term2 = -1*(aux*aux3 + aux2*aux4)

del aux, aux2, aux3, aux4

variables = [term1, term2]
########################################################################################################################
step=4
i='1950'
c_count=0
for c in cases:
    print(c)
    s_count = 0
    for s in seasons:
        print(s)

        # Las fechas se toman del periodo 1920-2020 basados en el DMI y N34 con ERSSTv5
        # Cuando se toman los periodos 1920-1949 y 1950_2020 las fechas que no pertencen
        # se excluyen de los composites en CompositeSimple()
        aux = xr.open_dataset(nc_date_dir + '1920_2020_' + s + '.nc')
        neutro = aux.Neutral
        case = aux[c]
        aux.close()

        mmonth = min_max_months[s_count]
        mmin = mmonth[0]
        mmax = mmonth[-1]

        neutro_comp_term1 = CompositeSimple(original_data=term1, index=neutro,
                                            mmin=mmin, mmax=mmax)
        neutro_comp_term2 = CompositeSimple(original_data=term2, index=neutro,
                                            mmin=mmin, mmax=mmax)

        comp_case_term1 = CompositeSimple(original_data=term1, index=case,
                                          mmin=mmin, mmax=mmax)
        comp_case_term2 = CompositeSimple(original_data=term2, index=case,
                                          mmin=mmin, mmax=mmax)

        comp_term1 = comp_case_term1 - neutro_comp_term1
        comp_term2 = comp_case_term2 - neutro_comp_term2

        if len(comp_term1) != 0:
            comp_term1 = xr.Dataset(
                data_vars=dict(
                    var=(['y', 'x'], comp_term1['var'][0, :, :].values)

                ),
                coords=dict(
                    lon=(['x'], comp_term1.lon),
                    lat=(['y'], comp_term1.lat),
                )
            )

            Plot(comp=comp_term1, levels=np.linspace(-1.2e-10,1.2e-10, 15), cmap=cbar,
                 dpi=dpi, save=save, step=1,
                 name_fig='term1_' + c + '_' + s + '_NSA',
                 title=r'$-\eta$' + 'D' + ' -  ' + s + '\n' + title_case[c_count],
                 contour0=False)

        if len(comp_term2) != 0:
            comp_term2 = xr.Dataset(
                data_vars=dict(
                    var=(['y', 'x'], comp_term2['var'][0, :, :].values)

                ),
                coords=dict(
                    lon=(['x'], comp_term2.lon),
                    lat=(['y'], comp_term2.lat),
                )
            )
            Plot(comp=comp_term2, levels=np.linspace(-1.2e-10,1.2e-10, 15), cmap=cbar,
                 dpi=dpi, save=save, step=1,
                 name_fig='term2_' + c + '_' + s + '_NSA',
                 title=r'$-\vec V_{\chi}.\nabla\eta$' + ' -  ' + s + '\n' + title_case[c_count],
                 contour0=False)
        s_count += 1
    c_count += 1
########################################################################################################################