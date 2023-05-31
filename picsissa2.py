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
################################################################################
#Plot
save = True
dpi = 300
################################################################################
# colorbar
cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC',
                                 '#B4E2DB',
                                 'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07',
                                 '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')
title_case = ['MAM', 'JJA', 'SON', 'DEF']
levels = np.linspace(0,1,11)
cmap = cbar_pp
name_fig = 'name',
significance = True
color_sig = 'k'

#
pred = xr.open_dataset('/home/luciano.andrian/doc/plotsissa/pred.nc')

fig_size = (9, 2.7)
extent = [270, 330, -58, 15]
xticks = np.arange(270, 340, 20)
yticks = np.arange(-55, 15, 15)

plt.rcParams['hatch.linewidth'] = 2
fig, axs = plt.subplots(nrows=1, ncols=4,
                        subplot_kw={'projection': ccrs.PlateCarree(
                            central_longitude=180)}, figsize=fig_size,
                        dpi=dpi)
crs_latlon = ccrs.PlateCarree()

for s in [1,2,3,4]:
    aux = pred.sel(estaciones=s)

    axs[s-1].set_extent(extent, crs=crs_latlon)
    im= axs[s-1].contourf(aux.lon, aux.lat,
                    aux.pred,
                    levels=levels, transform=crs_latlon,
                    cmap='YlGnBu')


    color_map = '#4B4B4B'
    axs[s-1].add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                             edgecolor=color_map)
    axs[s-1].add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    axs[s-1].coastlines(color=color_map, linestyle='-', alpha=1)

    axs[s-1].set_xticks(xticks, crs=crs_latlon)
    axs[s-1].set_yticks(yticks, crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    axs[s-1].xaxis.set_major_formatter(lon_formatter)
    axs[s-1].yaxis.set_major_formatter(lat_formatter)
    axs[s - 1].gridlines(linewidth=0.3, linestyle='-')
    axs[s-1].tick_params(labelsize=8)
    title = title_case[s-1]
    axs[s-1].set_title(title, fontsize=14)
    plt.tight_layout()

fig.subplots_adjust(right=0.925)
pos = fig.add_axes([0.935, 0.170, 0.012, 0.65])
cbar = fig.colorbar(im, cax=pos, pad=0.1)


if save:
    plt.savefig('/home/luciano.andrian/doc/plotsissa/figura8.tiff', dpi=300)
    plt.close()
else:
    plt.show()
