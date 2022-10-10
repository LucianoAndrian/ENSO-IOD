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
#from ENSO_IOD_Funciones import CompositeSimple

import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
########################################################################################################################
data_dir2 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/KS/'

save = True
step = 1
v = 'Ks'; print('\U0001F60D')
# Functions ############################################################################################################
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


def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         SA=False, dpi=100, save=True, step=1,contourf=True,
         name_fig='fig', title='title', contour0=True, color_map='k'):

    from numpy import ma
    import matplotlib.pyplot as plt

    comp_var = comp
    fig = plt.figure(figsize=(8, 3), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([30, 340, -80, 20], crs=crs_latlon)

    if contourf:
        im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                         levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    else:
        im = ax.contour(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                         levels=levels, transform=crs_latlon, cmap=cmap)  # , extend='max')

    if contour0:
        ax.contour(comp.lon, comp.lat, comp_var, levels=0,
                   transform=crs_latlon, colors='green', linewidths=1)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
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

def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients[0])
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients[0])
    dt = xrda - trend
    return dt
# cbar #################################################################################################################

cbar = colors.ListedColormap(['#FFE541', '#46B9F4'])
cbar.set_under('#E07700')
cbar.set_over('#016ABA')
cbar.set_bad(color='white')

########################################################################################################################
# seasons = ("Full_Season",'SON', 'ASO', 'JAS', 'JJA')
# min_max_months = [[7,11],[9,11], [8,10], [7,9], [6,8]]
seasons = ['JJA', 'SON']
min_max_months = [[6, 8], [9, 11]]

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

title_case = ['DMI-ENSO simultaneous positive phase ','DMI-ENSO simultaneous negative phase ',
              'DMI negative phase ', 'DMI positive phase ',
              'DMI pure positive phase ', 'DMI pure negative phase ',
              'ENSO positive phase ', 'ENSO negative phase ',
              'ENSO pure positive phase ', 'ENSO pure negative phase ']

# Open Data ############################################################################################################

u = xr.open_dataset(data_dir2 + 'u_mer_d_w_world.nc')

eta_grad_y = xr.open_dataset(data_dir2 + 'etay.nc')
eta_grad_y = eta_grad_y.rename({'meridional_gradient_of_absolute_vorticity':'var'})

u = u.interp(lon=eta_grad_y.lon.values, lat=eta_grad_y.lat.values)
ntime, nlats, nlons = u['var'].shape

########################################################################################################################
# Compute, plot and save ###############################################################################################
c_cases = 0
for c in cases:
    print(c)

    count = 0
    for s in seasons:
        print(s)

        aux = xr.open_dataset(nc_date_dir + '1920_2020_' + s + '.nc')

        case = aux[c]
        aux.close()
        mmonth = min_max_months[count]
        mmin = mmonth[0]
        mmax = mmonth[-1]

        #Ks del composite
        # usando estado basico de cada "case"
        comp = CompComp(eta_grad_y=eta_grad_y, u=u,
                        index=case, mmin=mmin, mmax=mmax)

        if len(comp) != 0:
            Plot(comp=comp, levels=[2,3,4], cmap=cbar,
                 dpi=200, save=save, step=step,
                 name_fig=v + '_' + c + '_' + s + '_' + '1950_2020',
                 title=v + ' - ' +  title_case[c_cases] + '- ' + '\n' + s + ' ' +
                       '1950-2020' + '\n' + 'Ks of Composite',
                 contour0=False, contourf=True, color_map='k')

        count += 1
    c_cases += 1
########################################################################################################################