import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
#-----------------------------------------------------------------------------#
data_dir = '/home/luciano.andrian/doc/plotsissa/Marzo2023/'
out_dir = '/home/luciano.andrian/doc/plotsissa/'
save=True
save_tiff=True
dpi=300
#-----------------------------------------------------------------------------#


cbar = colors.ListedColormap(['#FFAE0D', '#FFCC00', '#FFD249',
                              '#FFE72B', '#FFED74',
                              'white',
                               '#FF74E7', '#F02BC0', '#D76272',
                              '#B24576','#6F3460'][::-1])
cbar.set_over('#CB8A0D')
cbar.set_under('#552A46')
cbar.set_bad(color='white')

mycbar = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838',
                                '#E25E55', '#F28C89', '#FFCECC',
                                'white',
                                '#B3DBFF', '#83B9EB', '#5E9AD7',
                                '#3C7DC3', '#2064AF', '#014A9B'][::-1])
mycbar.set_over('#641B00')
mycbar.set_under('#012A52')
mycbar.set_bad(color='white')

cbar_snr = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4',
                                  '#6FFE9B',
                                  '#FFFFFF',
                                  '#FEB77E','#CA3E72','#782281',
                                  '#251255'])
cbar_snr.set_over('#251255')
cbar_snr.set_under('#070B4F')
cbar_snr.set_bad(color='white')

#-----------------------------------------------------------------------------#
# for v, v_n in zip(['hgt_1000', 'hgt_500', 'hgt_200', 'sf_200'],
#                   ['hgt', 'hgt', 'hgt', 'psi']):
for v, v_n in zip(['hgt_200'],['hgt']):
    if v_n == 'hgt':
        levels =  [-150, -75, -50, -25, -15, 0, 15, 25, 50, 75, 150]
    else:
        levels = [-2e7,-1.5e7, -1e7, -0.5e7, -0.25e7,
                  0, 0.25e7, 0.5e7, 1e7, 1.5e7, 2e7]
    # for c_bar, c_bar_name in zip([cbar, mycbar, cbar_snr],
    #                              ['propuesta', 'rb_propia', 'snr_propia']):
    for c_bar, c_bar_name in zip([cbar],
                                 ['propuesta']):

        data = xr.open_dataset(data_dir + v + '.nc')



        crs_latlon = ccrs.PlateCarree()

        fig = plt.figure(figsize=(9, 3.5), dpi=dpi)
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        ax.set_extent([70, 355, -75, 10], crs=crs_latlon)

        ax.contourf(data.lon, data.lat, data[v_n][0, :, :], levels=levels,
                    transform=crs_latlon, cmap=c_bar, extend='both',
                    zorder=2, alpha=0.7,
                    antialiased=True)

        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.LAND, color='#989898', alpha=0.7)
        ax.add_feature(cfeature.LAKES, color='skyblue')
        ax.add_feature(cfeature.BORDERS, linestyle='-')
        ax.add_feature(cfeature.OCEAN, color='skyblue', alpha=0.7)
        ax.add_feature(cfeature.RIVERS, edgecolor='skyblue', alpha=1)
        ax.add_feature(cfeature.STATES)

        ax.set_xticks(np.arange(80, 355, 15), crs=crs_latlon)
        ax.set_yticks(np.arange(-70, 20, 10), crs=crs_latlon)

        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)

        ax.tick_params(labelsize=7)
        # plt.title(title, fontsize=10)
        plt.tight_layout()
        if save:
            if save_tiff:
                plt.savefig(out_dir + v + '_' + c_bar_name + '.tiff', dpi=300)
                plt.close()
            else:
                plt.savefig(out_dir + v + '_' + c_bar_name + '.jpg', dpi=300)
                plt.close()

        else:
            plt.show()




