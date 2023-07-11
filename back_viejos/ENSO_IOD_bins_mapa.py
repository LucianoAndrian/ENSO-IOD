import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
#warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
warnings.filterwarnings('ignore')


def OpenDataSet(name, interp=False, lat_interp=None, lon_interp=None):

    if name == 'pp_gpcc':
        # GPCC2018
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_viejo/' + 'pp_gpcc.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        pp_gpcc = pp_gpcc.rename({'precip': 'var'})
        pp_gpcc = pp_gpcc.sel(time='1982-01-01')
        pp_gpcc = xr.where(np.isnan(pp_gpcc), np.nan, 1)

        return pp_gpcc

mask = OpenDataSet('pp_gpcc', interp=True,
                   lat_interp=np.linspace(-60,15,76),
                   lon_interp=np.linspace(275,330,56))

# box_name = ['S_SESA', 'S_SESA_exp', 'N_SESA', 'C_Brazil', 'Chile', 'Chile_sur']
# box_color =['Spectral_r', 'afmhot_r', 'BrBG', 'BrBG_r', 'RdBu', 'RdBu_r', 'Accent']
# box_lats = [[-35,-29],[-35,-20],[-29,-20],[-20,0],[-40,-30], [-60,-43]]
# box_lons = [[300, 310],[300,320],[300,320],[300,325],[285,290], [285,288]]
box_name = ['SESA', "N-SESA", 'Patagonia']
box_color =['Spectral_r', 'RdBu', 'Spectral']
box_lats = [[-39,-29], [-29,-17], [-55,-40]]
box_lons = [[296, 315], [296, 315], [288,300]]
#----------------------------------------------------------------------

import matplotlib.patches as patches
fig = plt.figure(figsize=(5, 6), dpi=300)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([270, 330, -60, 20], crs_latlon)

mask_var = mask['var']

for i in range(0,len(box_name)):
    aux = mask.sel(lon=slice(box_lons[i][0],box_lons[i][1]),
                   lat=slice(box_lats[i][0],box_lats[i][1]))
    ax.contourf(aux.lon, aux.lat, aux['var'],
                levels=[-0.5, 0.5], transform=crs_latlon, cmap=box_color[i],
                extend='both', alpha=0.5, label=box_name[i])


ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.LAND, facecolor='white')
ax.add_feature(cartopy.feature.BORDERS, facecolor='white')
ax.add_feature(cartopy.feature.COASTLINE)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labelsize=7)
#plt.title('cajas', fontsize=10)
plt.tight_layout()
plt.show()



