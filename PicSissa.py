"""
Composites OBSERVADOS
DMI standard
"""
################################################################################
import xarray as xr
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings("ignore")
from ENSO_IOD_Funciones import CaseComp
################################################################################
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
              'nc_composites_dates_no_ind_sst_anom/' #fechas
data_dir_t_pp = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'data_obs_d_w_c/' #T y PP ya procesados
out_dir = '/home/luciano.andrian/doc/plotsissa/'

#Plot
save = True
dpi = 300
sig_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_quantiles/'
################################################################################
scale_pp = np.linspace(-45, 45, 13)  #pp

# colorbar
cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC',
                                 '#B4E2DB',
                                 'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07',
                                 '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cases = ['DMI_un_neg', 'N34_un_neg', 'DMI_sim_neg']
title_case = ['IOD negativo',
              'La Niña',
              'IOD negativo + La Niña']

seasons = ['SON']
min_max_months = [[9,11]]

variables_t_p = ['pp_gpcc_d_w_c_1920-2020_0.25.nc']
v_count = 0
plt.rcParams['hatch.linewidth'] = 2

levels = scale_pp
cmap = cbar_pp
name_fig = 'name',
significance = True
color_sig = 'k'

v = variables_t_p[0]
c_count = 0
s_count = 0
s = 'SON'

data = xr.open_dataset(data_dir_t_pp + v)
fig_size = (9, 3.5)
extent = [270, 330, -58, 15]
xticks = np.arange(270, 340, 20)
yticks = np.arange(-55, 15, 15)

plt.rcParams['hatch.linewidth'] = 2
fig, axs = plt.subplots(nrows=1, ncols=3,
                        subplot_kw={'projection': ccrs.PlateCarree(
                            central_longitude=180)}, figsize=fig_size,
                        dpi=dpi)
crs_latlon = ccrs.PlateCarree()
for c in cases:
    comp1, num_case = \
        CaseComp(data, s, mmonth=min_max_months[s_count], c=c,
                 two_variables=False, nc_date_dir=nc_date_dir)

    data_sig = xr.open_dataset(
        sig_dir + v.split('_')[0] + '_' + v.split('_')[1] + '_' + c +
        '1950_2020_' + s + '_DMIbase.nc')

    comp1_i = comp1.interp(
        lon=data_sig.lon.values, lat=data_sig.lat.values)
    sig = comp1_i.where(
        (comp1_i < data_sig['var'][0]) | (comp1_i > data_sig['var'][1]))
    sig = sig.where(np.isnan(sig['var']), 0)

    axs[c_count].set_extent(extent, crs=crs_latlon)
    im=axs[c_count].contourf(comp1.lon, comp1.lat,
                    comp1['var'],
                    levels=levels, transform=crs_latlon,
                    cmap=cmap, extend='both')
    comp_sig = sig
    if significance:
        colors_l = [color_sig, color_sig]
        comp_sig_var = comp_sig['var']
        cs = axs[c_count].contourf(comp_sig.lon, comp_sig.lat, comp_sig_var,
                                   transform=crs_latlon, colors='none',
                                   hatches=["..", ".."], extend='lower')
        for i, collection in enumerate(cs.collections):
            collection.set_edgecolor(colors_l[i % len(colors_l)])

        for collection in cs.collections:
            collection.set_linewidth(0.)

        axs[c_count].contour(comp_sig.lon, comp_sig.lat, comp_sig_var,
                             linewidths=.8, levels=[0, 1],
                             transform=crs_latlon, colors='black')
    color_map = '#4B4B4B'
    axs[c_count].add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                             edgecolor=color_map)
    axs[c_count].add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    axs[c_count].coastlines(color=color_map, linestyle='-', alpha=1)
    axs[c_count].gridlines(linewidth=0.3, linestyle='-')
    axs[c_count].set_xticks(xticks, crs=crs_latlon)
    axs[c_count].set_yticks(yticks, crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    axs[c_count].xaxis.set_major_formatter(lon_formatter)
    axs[c_count].yaxis.set_major_formatter(lat_formatter)
    axs[c_count].tick_params(labelsize=8)
    title = title_case[c_count]
    axs[c_count].set_title(title, fontsize=15)
    plt.tight_layout()

    c_count += 1

fig.subplots_adjust(right=0.925)
pos = fig.add_axes([0.935, 0.2, 0.012, 0.6])
cbar = fig.colorbar(im, cax=pos, pad=0.1)
if save:
    plt.savefig('/home/luciano.andrian/doc/plotsissa/figura6.tiff', dpi=300)
    plt.close()
else:
    plt.show()

