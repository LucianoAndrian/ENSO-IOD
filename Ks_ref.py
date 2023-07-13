import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
import cartopy.feature
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
from ENSO_IOD_Funciones import CompositeSimple
#from ENSO_IOD_Funciones import CompositeSimple

import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
################################################################################
data_dir2 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'


save = True
dpi = 300
step = 1

true_dipole = [False, True]
v = 'Ks'

# Functions ####################################################################
def WhatDMI(true_dipole):
    if true_dipole:
        out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/' \
                  'composite/dmi_true_dipole/'
        nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                      'nc_composites_dates/'
        dmi_index='DMI_true_dipole'
    else:
        out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/' \
                  'composite/dmi_standard/'
        nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                      'nc_composites_dates_no_ind_sst_anom/'
        dmi_index = 'DMI_standard'
    return out_dir, nc_date_dir, dmi_index

# Comp y Ks
def CompComp(eta_grad_y, u, index, mmin, mmax):
    data_y_eta = CompositeSimple(original_data=eta_grad_y, index=index,
                                 mmin=mmin, mmax=mmax)

    data_y_u = CompositeSimple(original_data=u, index=index,
                               mmin=mmin, mmax=mmax)

    Rt2 = np.transpose(np.tile(6.378e6 * np.cos(u.lat.values * np.pi / 180),
                               [359, 1]), [1, 0])

    comp = np.sqrt(data_y_eta['var'][0,:,:] / data_y_u['var']) * Rt2

    comp = xr.where(comp < 0, np.nan, comp)

    return comp

def extractlevel(x):
    aux = [int(s) for s in x if s.isdigit()]
    return''.join(str(num) for num in aux)


def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         SA=False, dpi=100, save=True, step=1,contourf=True,
         name_fig='fig', title='title', contour0=True, color_map='k',
         out_dir='out_dir'):

    from numpy import ma
    import matplotlib.pyplot as plt

    comp_var = comp
    fig = plt.figure(figsize=(8, 3), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([0, 360, -80, 20], crs=crs_latlon)

    if contourf:
        im = ax.contourf(comp.lon[::step], comp.lat[::step],
                         comp_var[::step, ::step],
                         levels=levels, transform=crs_latlon,
                         cmap=cmap, extend='both')
    else:
        im = ax.contour(comp.lon[::step], comp.lat[::step],
                        comp_var[::step, ::step],
                         levels=levels, transform=crs_latlon, cmap=cmap)

    if contour0:
        ax.contour(comp.lon, comp.lat, comp_var, levels=0,
                   transform=crs_latlon, colors='green', linewidths=1)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(0, 360, 30), crs=crs_latlon)
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

# cbar #########################################################################
cbar = colors.ListedColormap(['#FFE541', '#46B9F4'])
cbar.set_under('#E07700')
cbar.set_over('#016ABA')
cbar.set_bad(color='white')

################################################################################
seasons = ['SON']
min_max_months = [[9, 11]]

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_un_pos', 'DMI_un_neg',
         'N34_un_pos', 'N34_un_neg']

title_case = ['IOD-ENSO simultaneous positive phase ',
              'IOD-ENSO simultaneous negative phase ',
              'IOD pure positive phase ', 'IOD pure negative phase ',
              'ENSO pure positive phase ', 'ENSO pure negative phase ']


for td in true_dipole:
    out_dir, nc_date_dir, dmi_index = WhatDMI(td)

    for v in ['UV200']:
        # open data
        u = xr.open_dataset(data_dir2 + 'u_' + v + '_w_detrend.nc')
        eta_grad_y = xr.open_dataset(data_dir2 + 'etay_' + v + '.nc')
        eta_grad_y = eta_grad_y.rename(
            {'meridional_gradient_of_absolute_vorticity': 'var'})

        u = u.interp(lon=eta_grad_y.lon.values, lat=eta_grad_y.lat.values)
        ntime, nlats, nlons = u['var'].shape
        ########################################################################
        # Compute, plot and save ###############################################
        for c, c_cases in zip(cases, range(0, len(cases))):
            print(c)
            for s, count in zip(seasons, range(0, len(seasons))):
                print(s)

                aux = xr.open_dataset(nc_date_dir + '1920_2020_' + s + '.nc')

                case = aux[c]
                aux.close()
                mmonth = min_max_months[count]
                mmin = mmonth[0]
                mmax = mmonth[-1]

                # Ks del composite
                # usando estado basico de cada "case"
                comp = CompComp(eta_grad_y=eta_grad_y, u=u,
                                index=case, mmin=mmin, mmax=mmax)

                if len(comp) != 0:
                    Plot(comp=comp, levels=[2, 3, 4], cmap=cbar,
                         dpi=dpi, save=save, step=step,
                         name_fig='KS' + v + '_' + c + '_' + s + '_' +
                                  '1940_2020_' + dmi_index,
                         title= 'Ks of Composite - ' + extractlevel(v) + 'hPa'
                                + '\n' + dmi_index + '\n'
                                + title_case[c_cases] + ' - ' + s + '1940-2020',
                         contour0=False, contourf=True, color_map='grey',
                         out_dir=out_dir)
########################################################################################################################
