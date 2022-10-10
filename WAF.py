import xarray as xr
import numpy as np


def manipular_nc(archivo,variable):
    #con xarray
    dataset = xr.open_dataset(archivo)
    var_out = dataset[variable].values
    lon = dataset.lon.values
    lat = dataset.lat.values
    dataset.close()
    return var_out, lon, lat

def c_diff(arr, h, dim, cyclic=False):  # compute derivate of array variable respect to h associated to dim
    # adapted from kuchaale script
    ndim = arr.ndim
    lst = [i for i in range(ndim)]

    lst[dim], lst[0] = lst[0], lst[dim]
    rank = lst
    arr = np.transpose(arr, tuple(rank))

    if ndim == 3:
        shp = (arr.shape[0] - 2, 1, 1)
    elif ndim == 4:
        shp = (arr.shape[0] - 2, 1, 1, 1)

    d_arr = np.copy(arr)
    if not cyclic:
        d_arr[0, ...] = (arr[1, ...] - arr[0, ...]) / (h[1] - h[0])
        d_arr[-1, ...] = (arr[-1, ...] - arr[-2, ...]) / (h[-1] - h[-2])
        d_arr[1:-1, ...] = (arr[2:, ...] - arr[0:-2, ...]) / np.reshape(h[2:] - h[0:-2], shp)

    elif cyclic:
        d_arr[0, ...] = (arr[1, ...] - arr[-1, ...]) / (h[1] - h[-1])
        d_arr[-1, ...] = (arr[0, ...] - arr[-2, ...]) / (h[0] - h[-2])
        d_arr[1:-1, ...] = (arr[2:, ...] - arr[0:-2, ...]) / np.reshape(h[2:] - h[0:-2], shp)

    d_arr = np.transpose(d_arr, tuple(rank))

    return d_arr


# [psiclm, lon, lat] = manipular_nc(ruta + out_var[0] + '.nc', nc_var)
#
# [psiaa, lonpsi, latpsi] = manipular_nc(ruta + out_var[1] + '.nc', nc_var)  # psia [1 nlat nlon]

WAF(psiclm=All_neutral['var'].values.reshape(1,len(All_neutral.lat),len(All_neutral.lon)),
    psiaa=aux['var'].values.reshape(1,len(All_neutral.lat),len(All_neutral.lon)),
    lon=aux.lon.values,
    lat=aux.lat.values)


def WAF(psiclm, psiaa, lon, lat,reshape=True, variable='var'):
    #agregar xr=True

    if reshape:
        psiclm=psiclm[variable].values.reshape(1,len(psiclm.lat),len(psiclm.lon))
        psiaa = psiaa[variable].values.reshape(1, len(psiaa.lat), len(psiaa.lon))

    lon=lon.values
    lat=lat.values

    [xxx, nlats, nlons] = psiaa.shape  # get dimensions
    a = 6400000
    coslat = np.cos(lat * np.pi / 180)

    # climatological wind at psi level
    dpsiclmdlon = c_diff(psiclm, lon, 2)
    dpsiclmdlat = c_diff(psiclm, lat, 1)

    uclm = -1 * dpsiclmdlat
    vclm = dpsiclmdlon
    magU = np.sqrt(np.add(np.power(uclm, 2), np.power(vclm, 2)))

    dpsidlon = c_diff(psiaa, lon, 2)
    ddpsidlonlon = c_diff(dpsidlon, lon, 2)
    dpsidlat = c_diff(psiaa, lat, 1)
    ddpsidlatlat = c_diff(dpsidlat, lat, 1)
    ddpsidlatlon = c_diff(dpsidlat, lon, 2)

    termxu = dpsidlon * dpsidlon - psiaa * ddpsidlonlon
    termxv = dpsidlon * dpsidlat - ddpsidlatlon * psiaa
    termyv = dpsidlat * dpsidlat - psiaa * ddpsidlatlat

    # 0.2101 is the scale of p VER!!!
    coeff1 = np.transpose(np.tile(coslat, (nlons, 1))) * (0.2101) / (2 * magU)
    # x-component
    px = coeff1 / (a * a * np.transpose(np.tile(coslat, (nlons, 1)))) * (
            uclm * termxu / np.transpose(np.tile(coslat, (nlons, 1))) + vclm * termxv)
    # y-component
    py = coeff1 / (a * a) * (uclm / np.transpose(np.tile(coslat, (nlons, 1))) * termxv + vclm * termyv)
    return px, py

from numpy import ma
sdadas
fig = plt.figure(figsize=(7, 3.5), dpi=200)
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
crs_latlon = ccrs.PlateCarree()
ax.set_extent([0, 359, -90, 30], crs=crs_latlon)

im = ax.contourf(aux.lon, aux.lat, aux['var'].values, levels=np.arange(-1e7, 1e7, 0.25e7),
                 transform=crs_latlon, cmap='RdBu_r', extend='both')


# plot plumb fluxes and save again
# mask wind data to only show the 40% stronger fluxes.
Q60 = np.percentile(np.sqrt(np.add(np.power(px, 2), np.power(py, 2))), 60)
M = np.sqrt(np.add(np.power(px, 2), np.power(py, 2))) < Q60
# mask array
px_mask = ma.array(px, mask=M)
py_mask = ma.array(py, mask=M)
# plot vectors
lons, lats = np.meshgrid(aux.lon.values, aux.lat.values)
ax.quiver(lons[::20, ::20], lats[::20, ::20], px_mask[0,::20, ::20],
          py_mask[0, ::20, ::20], width=15e-4, headwidth=3,  # headwidht (default3)
          headlength=2.2, transform=crs_latlon)  # (default5))


cb = plt.colorbar(im, fraction=0.042, pad=0.035, shrink=0.8)
cb.ax.tick_params(labelsize=8)
ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
ax.add_feature(cartopy.feature.COASTLINE)
# ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
ax.set_xticks(np.arange(30, 330, 60), crs=crs_latlon)
ax.set_yticks(np.arange(-90, 30, 20), crs=crs_latlon)
lon_formatter = LongitudeFormatter(zero_direction_label=True)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.tick_params(labelsize=7)
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
plt.show()


plt.title(str(title) + ' - ' + str(main_month_name) + '  ' + str(fase.split(' ', 1)[1]), fontsize=10)
plt.figtext(0.5, 0.01, number_events, ha="center", fontsize=10,
            bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5})
plt.tight_layout()

if save:
    plt.savefig(name_fig + str(main_month_name) + '_' + str(fase.split(' ', 1)[1]) + '.jpg')
    plt.close()
else:
    plt.show()

