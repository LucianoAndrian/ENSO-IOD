########################################################################################################################
import xarray as xr
import numpy as np
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
import cartopy.feature as cfeature
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import cartopy.io.img_tiles as cimgt
import warnings
warnings.filterwarnings('ignore')
import matplotlib.patches as mpatches

########################################################################################################################

out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/CONGREMET/'
save = False
dpi = 100
def Plot(dpi=100, save=True,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(8.5, 5), dpi=dpi)
    #ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    stamen_terrain = cimgt.Stamen('terrain-background')
    ax = plt.axes(projection=stamen_terrain.crs)
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([30, 130, -30, 30], crs_latlon)

    ax.text(53, -1, 'Oeste', size=25, transform=crs_latlon)
    ax.text(94.5, -6, 'Este', size=25, transform=crs_latlon)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    bodr = cartopy.feature.NaturalEarthFeature(category='cultural',
                                               name='admin_0_boundary_lines_land', scale='50m', facecolor='none',
                                               alpha=0.7)
    ocean = cartopy.feature.NaturalEarthFeature('physical', 'ocean',
                                                scale='50m', edgecolor='none', facecolor=cfeature.COLORS['water'])
    land = cartopy.feature.NaturalEarthFeature('physical', 'land', \
                                               scale='50m', edgecolor='k', facecolor=cfeature.COLORS['land'])
    lakes = cartopy.feature.NaturalEarthFeature('physical', 'lakes', \
                                                scale='50m', edgecolor='b', facecolor=cfeature.COLORS['water'])
    rivers = cartopy.feature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', \
                                                 scale='50m', edgecolor='b', facecolor='none')
    # ax.add_feature(land, facecolor='beige')
    # ax.add_feature(ocean, linewidth=0.2)
    #ax.add_feature(lakes, alpha=0.5)
    #ax.add_feature(rivers, linewidth=0.5)
    #ax.add_feature(bodr, linestyle='--', edgecolor='k', alpha=0.5)
    ax.add_patch(mpatches.Rectangle(xy=[50, -10], width=20, height=20,
                                    facecolor='gray',
                                    alpha=0.5, edgecolor='black', linewidth=2,
                                    transform=ccrs.PlateCarree())
                 )

    ax.add_patch(mpatches.Rectangle(xy=[90, -10], width=20, height=10,
                                    facecolor='gray',
                                    alpha=0.5, edgecolor='black', linewidth=2,
                                    transform=ccrs.PlateCarree())
                 )
    ax.gridlines(crs=crs_latlon, linewidth=0.4, linestyle='-', color='gray', alpha=1)
    ax.set_xticks(np.arange(30, 130, 20), crs=crs_latlon)
    ax.set_yticks(np.arange(-30, 30, 10), crs=crs_latlon)
    ax.add_image(stamen_terrain, 8)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=10)

    #plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg',dpi=dpi)
        plt.close()
    else:
        plt.show()

Plot(save=True, dpi=300, name_fig='mapa_DMI')