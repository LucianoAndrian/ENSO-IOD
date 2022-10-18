import xarray as xr
from ENSO_IOD_Funciones import SelectNMMEFiles
########################################################################################################################
dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
out_dir = '/pikachu/datos/luciano.andrian/cases_fields/'

v = 'prec'
# funciones ############################################################################################################
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds
########################################################################################################################
# HINDCAST ------------------------------------------------------------------------------------------------------------#
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_hc, All=True)
files = sorted(files, key=lambda x: x.split()[0])

#abriendo todos los archivos
data = xr.open_mfdataset(files, decode_times=False) #xr no entiende la codificacion de Leads, r y las fechas
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data = data.sel(L=[0.5, 1.5, 2.5, 3.5]) # Solo leads 0 1 2 3
data['L'] = [0,1,2,3]
data = xr.decode_cf(fix_calendar(data)) # corrigiendo fechas
data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
data = data.sel(r=slice(1,24))

seasons_data = data.rolling(time=3, center=True).mean()
#
# son = seasons_data.sel(time=seasons_data.time.dt.month.isin(10))*30
# son = son.mean(['time', 'r', 'L'])
# son = son.compute()

# jja
jja = seasons_data.sel(time=seasons_data.time.dt.month.isin(7))*30
jja = jja.mean(['time', 'r', 'L'])
jja=jja.compute()
jja.to_netcdf(out_dir + v + '_mask_to_dry_season.nc')

# K-MEANS ##############################################################################################################

import numpy as np
from threadpoolctl import threadpool_limits
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from ENSO_IOD_Funciones import SelectNMMEFiles
import xarray as xr

def KM(data, n_clusters, centroids=False, reshape_dim=[56,76]):
    with threadpool_limits(limits=1):#con 1 mas rapido q con 40...
        kmeans = KMeans(n_clusters=n_clusters, n_init=100, random_state=666).fit(data)

    y_pred = kmeans.predict(data)
    if centroids:
        return y_pred.reshape(reshape_dim[0], reshape_dim[1]).T, kmeans.cluster_centers_
    else:
        return y_pred.reshape(reshape_dim[0], reshape_dim[1]).T


def MakeMask(DataArray, dataname='mask'):
    import regionmask
    mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(DataArray)
    mask = xr.where(np.isnan(mask), mask, 1)
    mask = mask.to_dataset(name=dataname)
    return mask


def CFSv2Data(v):
    out_data_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
    dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
    try:
        data_cfsv2 = xr.open_dataset(out_data_dir + v + '_data_cfsv2_noRoll.nc')
    except:
        files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                                dir=dir_hc, All=True)
        files = sorted(files, key=lambda x: x.split()[0])

        # abriendo todos los archivos
        data_cfsv2 = xr.open_mfdataset(files, decode_times=False).sel(
            L=[0.5, 1.5, 2.5, 3.5])  # xr no entiende la codificacion de Leads, r y las fechas
        data_cfsv2 = data_cfsv2.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
        data_cfsv2['L'] = [0, 1, 2, 3]
        data_cfsv2 = xr.decode_cf(fix_calendar(data_cfsv2))  # corrigiendo fechas
        data_cfsv2 = data_cfsv2.sel(lon=slice(275, 330), lat=slice(-60, 15))
        data_cfsv2 = data_cfsv2.sel(r=slice(1, 24))
        data_cfsv2.to_netcdf(out_data_dir + v + '_data_cfsv2_noRoll.nc')
        data_cfsv2 = data_cfsv2.compute()
    return data_cfsv2

def KmMaks(data_cfsv2, v, l, pb):
    data2 = data_cfsv2.sel(L=l)
    # data2 = data2.sel(lon=slice(290, 330), lat=slice(-60, 0))
    data_clim = data2.groupby('time.month').mean()
    data_clim = data_clim.mean(['r'])
    data_clim = data_clim.rename({'month': 'time'})
    data_clim.load()
    data_clim *= 30
    # data_clim *= mask
    # if l == 0:
    #     pb = 0
    # else:
    #     pb = 2

    # pp2 = data_clim.sel(lon=slice(pre_box_lon[pb][0], pre_box_lon[pb][1]),
    #                     lat=slice(pre_box_lat[pb][0], pre_box_lat[pb][1]))
    pp2 = data_clim
    #### k-means #####
    X = pp2.stack(new=['lon', 'lat'])[v].T
    X[np.where(np.isnan(X))] = -99

    for n_c in [15]:
        reshape_lat_dim = len(range(pre_box_lat[pb][0], pre_box_lat[pb][1])) + 1
        reshape_lon_dim = len(range(pre_box_lon[pb][0], pre_box_lon[pb][1])) + 1

        clusters, centroids = KM(data=X, n_clusters=n_c, centroids=True)

        xr_cluster_label = MakeMask(pp2, dataname='cluster')
        xr_cluster_label['cluster'].values = clusters

        xr_cluster_label *= MakeMask(xr_cluster_label, 'cluster')

    return xr_cluster_label

#----------------------------------------------------------------------------------------------------------------------#
# Region seleccionada con k-means a paritr del CFSv2 ------------------------------------------------------------------#
pre_box = [0, 1, 2]
pre_box_lat = [[-60, 0], [-60, -20], [-60, -20]]
pre_box_lon = [[290, 330], [290, 330], [280, 290]]
#----------------------------------------------------------------------------------------------------------------------#

v = 'prec'
data_cfsv2 = CFSv2Data(v)
xr_cluster_label = KmMaks(data_cfsv2, v=v, l=0, pb=0)
plt.imshow(xr_cluster_label.cluster);plt.colorbar();plt.show()

mask_cluster = xr.where(xr_cluster_label == 7, 1, np.nan)#.sel(lat=slice(-40,-20), lon=slice(290, 320))
plt.imshow(mask_cluster.cluster);plt.show()
mask_cluster = mask_cluster.where(mask_cluster.lat>-20)
mask_cluster = xr.where(np.isnan(mask_cluster), 1, np.nan)

#primero, usando 15 cluster y mask el 7
#otro, mas grande pero menos preciso, usando 5 cluster y mask el 4