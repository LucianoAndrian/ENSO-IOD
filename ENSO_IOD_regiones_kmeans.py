"""
Probando K-means para delimitar regiones (EOF no funcionó)
"""
"""
Notas: 
Funciona mejor que EOF, mucho mejor. 
Complicaciones: 
- Las métricas no son del todo apropiadas para este caso, sirve más mirar los campos de clusters
- Las zonas con topografía alta complican mucho en el modelo. En cierto punto agregar clusters no ayuda a definir 
mas regiones sino sólo dividir más y más las regiones montañosas.
- No funciona con NAs, x lo q usar metricas para gpcc no es eficiente ya que previo se debe dar un valor numerico
a los Nas que serán puestos en un cluster... y cubren gran parte del dominio

Agregado:
Nueva selección a partir de "pre Box", regiones más chicas que ayudan a que el algoritmo funcione mejor
y permite delimitar mejor regiones como S-SESA. Aún es algo problematica la región de Chile

Nueva mascara de oceanos, en graficos sólo hacia falta usar cartopy, para calculos ahora se pueden enmascarar 
usando regionmask. Funcion MakeMask en ENSO_IOD_Funciones.py
"""
########################################################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits
#from ENSO_IOD_Funciones import SelectNMMEFiles, MakeMask
########################################################################################################################
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Regiones/'
save=True
dpi=150
########################################################################################################################
# Funciones ############################################################################################################
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

def xrFieldTimeDetrend_sst(xrda, dim, deg=1):
    # detrend along a single dimension
    aux = xrda.polyfit(dim=dim, deg=deg)
    trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    dt = xrda - trend
    return dt, aux

def OpenDataSet(name, interp=False, lat_interp=None, lon_interp=None):

    if name == 'pp_gpcc':
        # GPCC2018
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' +
                              'pp_gpcc_d_w_c_1950-2020_1.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        #pp_gpcc = pp_gpcc.rename({'precip': 'var'})

        return pp_gpcc
    elif name == 'mask':
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_viejo/' + 'pp_gpcc.nc')
        mask = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            mask = aux.interp(lon=lon_interp, lat=lat_interp)

        mask = mask.mean('time')
        mask = xr.where(np.isnan(mask), mask, 1)
        mask = mask.rename({'precip': 'mask'})
        return mask
    elif name == 'pp_gpcc_or':
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_no_detrend/' +
                              'pp_gpcc_v2020_0.25.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        return pp_gpcc
    elif name == 'cmap':
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_no_detrend/' +
                              'precip.mon.mean.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        return pp_gpcc

def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, OceanMask=False,
         name_fig='fig', title='title'):

    comp_var = comp
    fig = plt.figure(figsize=(5, 6), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([275,330, -60,15], crs_latlon)

    im = ax.contourf(comp.lon, comp.lat, comp_var['cluster'],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    if OceanMask:
        ax.add_feature(cartopy.feature.OCEAN, zorder=100, edgecolor='k', facecolor='white')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(275, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 15, 10), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)

    plt.title(title, fontsize=12)
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

def KM(data, n_clusters, centroids=False, reshape_dim=[56,76]):
    with threadpool_limits(limits=1):#con 1 mas rapido q con 40...
        kmeans = KMeans(n_clusters=n_clusters, n_init=100, random_state=666).fit(data)

    y_pred = kmeans.predict(data)
    if centroids:
        return y_pred.reshape(reshape_dim[0], reshape_dim[1]).T, kmeans.cluster_centers_
    else:
        return y_pred.reshape(reshape_dim[0], reshape_dim[1]).T

def PlotCentroids(serie, dpi, title, label, name_fig, save, ymax=400):
    fig = plt.figure(figsize=(6, 3), dpi=dpi)
    ax = fig.add_subplot(111)
    ax.plot(serie, linewidth=2, color='dodgerblue', label='Cluster ' + label)
    #ax.plot(serie2, linewidth=2, color='indianred', label='CFSv2')
    ax.tick_params(labelsize=7)
    ax.set_xlabel('Months')
    ax.set_xticks(range(0,12), ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], size=10)
    plt.ylim(0,ymax)
    plt.title(title, fontsize=10)
    plt.legend()
    plt.grid()
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

def SilhScore(X, title='title', name_fig='fig', save=False):
    from sklearn.metrics import silhouette_score

    with threadpool_limits(limits=1):
        kmeans_per_k = [KMeans(n_clusters=k, n_init=100, random_state=666).fit(X) # misma config. q KM
                        for k in range(1, 15)]

        silhouette_scores = [silhouette_score(X, model.labels_)
                             for model in kmeans_per_k[1:]]

    plt.figure(figsize=(8, 3))
    plt.plot(range(2, 15), silhouette_scores, "bo-")
    plt.xlabel("$k$", fontsize=14)
    plt.ylabel("Silhouette score", fontsize=14)
    plt.title(title)

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

def MakeMask(DataArray, dataname='mask'):
    import regionmask
    mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(DataArray)
    mask = xr.where(np.isnan(mask), mask, 1)
    mask = mask.to_dataset(name=dataname)
    return mask
########################################################################################################################
color=[plt.cm.tab20b(0),plt.cm.tab20b(1),plt.cm.tab20b(2),plt.cm.tab20b(3),
       plt.cm.tab20c(0),plt.cm.tab20c(1),plt.cm.tab20c(2),plt.cm.tab20c(3),
       '#8CEFDE','white','white', '#F9FD7C',
       plt.cm.tab20c(7), plt.cm.tab20c(6), plt.cm.tab20c(5), plt.cm.tab20c(4),
       plt.cm.tab20b(15), plt.cm.tab20b(14), plt.cm.tab20b(13), plt.cm.tab20b(12)]
cmap=colors.ListedColormap(color)

pre_box = [0, 1, 2]
pre_box_lat = [[-60, 0], [-60, -20], [-60, -20]]
pre_box_lon = [[290, 330], [290, 330], [280, 290]]

########################################################################################################################
# Observado -----------------------------------------------------------------------------------------------------------#
for b in ['pp_gpcc_or', 'cmap']:

    pp = OpenDataSet(b, interp=True,
                     lat_interp=np.linspace(-60, 15, 76),
                     lon_interp=np.linspace(275, 330, 56))

    # Mismo período que el CFSv2 (hindcast)
    pp = pp.sel(time=pp.time.dt.year.isin(range(1982, 2012)))

    # Climatologia
    pp = pp.groupby('time.month').mean()
    pp = pp.rename({'month': 'time'})
    pp = pp.rename({'precip': 'var'})
    #pp = pp.sel(lon=slice(290,330), lat=slice(-60,0))
    if b == 'cmap':
        pp *= 30

    #### k-means #####
    X = pp.stack(new=['lon', 'lat'])['var'].T
    X[np.where(np.isnan(X))] = -99

    SilhScore(X, title=b, name_fig= b + 'SilhScore', save=True)

    for n_c in [7, 8, 9, 10]:

        clusters, centroids = KM(data=X, n_clusters=n_c, centroids=True)
        xr_cluster_label = MakeMask(pp, dataname='cluster')
        xr_cluster_label['cluster'].values = clusters

        #xr_cluster_label *= MakeMask(pp, dataname='cluster') # No es necesario ahora
        xr_cluster_label = xr_cluster_label.interp(lat=np.linspace(-60, 15, 76 * 8),
                                                   lon=np.linspace(275, 330, 56 * 8),
                                                   method='nearest')

        Plot(comp=xr_cluster_label, levels=np.arange(0, n_c, 1), cmap='tab20',
             OceanMask=True,
             title=b + ' - Obs ' + '- n_clusters: ' + str(n_c),
             name_fig=b +'_Obs' + '_n_c' + str(n_c) + 'masked',
             dpi=dpi, save=save)

        for n_c2 in np.arange(0, n_c, 1):
            PlotCentroids(serie=centroids[n_c2, :], label=str(n_c2),
                          title=b + ' - Centroids cluster ' + str(n_c2) + '\n' +
                                '- total Clusters' + str(n_c),
                          name_fig=b + '_centroid_n_c2' + str(n_c2) + '_n_c' + str(n_c),
                          dpi=100, save=save, ymax=500)

#----------------------------------------------------------------------------------------------------------------------#
# CFSv2 ---------------------------------------------------------------------------------------------------------------#
#from ENSO_IOD_Funciones import SelectNMMEFiles

dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
dir_rt = '/pikachu/datos/osman/nmme/monthly/real_time/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Regiones/'
v = 'prec'

files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_hc, All=True)
files = sorted(files, key=lambda x: x.split()[0])

#abriendo todos los archivos
data = xr.open_mfdataset(files, decode_times=False).sel(L=[0.5, 1.5, 2.5, 3.5]) #xr no entiende la codificacion de Leads, r y las fechas
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data['L'] = [0,1,2,3]
data = xr.decode_cf(fix_calendar(data)) # corrigiendo fechas
data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
data = data.sel(r=slice(1,24))
aux_data = data.mean('r').polyfit(dim='time', deg=1)
# aux_trend = xr.polyval(data['time'], aux_data.prec_polyfit_coefficients)
# data_d = data - aux_trend

# Por leadtimes -------------------------------------------------------------------------------------------------------#
for l in [0,1,2,3]:
    data2 = data.sel(L=l)
    #data2 = data2.sel(lon=slice(290, 330), lat=slice(-60, 0))
    data_clim = data2.groupby('time.month').mean()
    data_clim = data_clim.mean(['r'])
    data_clim = data_clim.rename({'month': 'time'})
    data_clim.load()
    data_clim *= 30
    #data_clim *= mask

    #### k-means #####
    X = data_clim.stack(new=['lon', 'lat'])['prec'].T
    SilhScore(X, title=b, name_fig=str(l) + 'SilhScore', save=True)
    for n_c in [7, 8, 9, 10]:
        clusters, centroids = KM(data=X, n_clusters=n_c, centroids=True)
        xr_cluster_label = MakeMask(data_clim, dataname='cluster')
        xr_cluster_label['cluster'].values = clusters

        #xr_cluster_label *= MakeMask(data_clim, dataname='cluster')
        xr_cluster_label = xr_cluster_label.interp(lat=np.linspace(-60, 15, 76 * 8),
                                                   lon=np.linspace(275, 330, 56 * 8),
                                                   method='nearest')

        Plot(comp=xr_cluster_label, levels=np.arange(0, n_c, 1), cmap='tab20',
             title='CFSv2 Leadtime: ' + str(l) + ' - n_clusters: ' + str(n_c),
             name_fig='CFSv2_L' + str(l) + '_n_c' + str(n_c) + 'masked',
             dpi=dpi, save=save, OceanMask=True)

        for n_c2 in np.arange(0, n_c, 1):
            PlotCentroids(serie=centroids[n_c2, :], label=str(n_c2),
                          title='CFSv2 Leadtime '  + str(l) + 'Centroids cluster ' + str(n_c2) + '\n' +
                                '- total Clusters' + str(n_c),
                          name_fig='CFSv2_L'+ str(l) + '_centroid_n_c2' + str(n_c2) + '_n_c' + str(n_c),
                          dpi=100, save=save, ymax=500)
########################################################################################################################
# pre Box
########################################################################################################################
# Observado -----------------------------------------------------------------------------------------------------------#
for b in ['pp_gpcc_or', 'cmap']:

    pp = OpenDataSet(b, interp=True,
                     lat_interp=np.linspace(-60, 15, 76),
                     lon_interp=np.linspace(275, 330, 56))

    # Mismo período que el CFSv2 (hindcast)
    pp = pp.sel(time=pp.time.dt.year.isin(range(1982, 2012)))

    # Climatologia
    pp = pp.groupby('time.month').mean()
    pp = pp.rename({'month': 'time'})
    pp = pp.rename({'precip': 'var'})
    #pp = pp.sel(lon=slice(290,330), lat=slice(-60,0))
    if b == 'cmap':
        pp *= 30

    for pb in pre_box:
        pp2 = pp.sel(lon=slice(pre_box_lon[pb][0], pre_box_lon[pb][1]),
                     lat=slice(pre_box_lat[pb][0], pre_box_lat[pb][1]))
        #### k-means #####
        X = pp2.stack(new=['lon', 'lat'])['var'].T
        X[np.where(np.isnan(X))] = -99

        SilhScore(X, title=b + ' - pre Box: ' + str(pb), name_fig= b + 'SilhScore_Box' + str(pb), save=True)

        for n_c in [4, 5, 6, 7, 8]:
            reshape_lat_dim = len(range(pre_box_lat[pb][0], pre_box_lat[pb][1])) + 1
            reshape_lon_dim = len(range(pre_box_lon[pb][0], pre_box_lon[pb][1])) + 1

            clusters, centroids = KM(data=X, n_clusters=n_c, centroids=True,
                                     reshape_dim=[reshape_lon_dim, reshape_lat_dim])

            xr_cluster_label = MakeMask(pp, dataname='cluster').sel(lon=slice(pre_box_lon[pb][0], pre_box_lon[pb][1]),
                                                                    lat=slice(pre_box_lat[pb][0],
                                                                              pre_box_lat[pb][1])).copy()
            xr_cluster_label['cluster'].values = clusters

            xr_cluster_label = xr_cluster_label.interp(
                lat=np.linspace(pre_box_lat[pb][0], pre_box_lat[pb][1], reshape_lat_dim * 8),
                lon=np.linspace(pre_box_lon[pb][0], pre_box_lon[pb][1], reshape_lon_dim * 8),
                method='nearest')

            Plot(comp=xr_cluster_label, levels=np.arange(0, n_c, 1), cmap='tab20',
                 title=b + ' - Obs ' + '- n_clusters: ' + str(n_c) + '\n' + 'pre Box ' + str(pb),
                 name_fig=b + '_Obs' + '_n_c' + str(n_c) + '_preBox-' + str(pb),
                 dpi=dpi, save=save, OceanMask=True)

            for n_c2 in np.arange(0, n_c, 1):
                PlotCentroids(serie=centroids[n_c2, :], label=str(n_c2),
                              title=b + ' - Centroids cluster ' + str(n_c2) + '\n' +
                                    '- total Clusters' + str(n_c),
                              name_fig=b + '_centroid_n_c2' + str(n_c2) +
                                       '_n_c' + str(n_c) + '_preBox-' + str(pb),
                              dpi=100, save=save, ymax=500)

for l in [0,1,2,3]:
    data2 = data.sel(L=l)
    #data2 = data2.sel(lon=slice(290, 330), lat=slice(-60, 0))
    data_clim = data2.groupby('time.month').mean()
    data_clim = data_clim.mean(['r'])
    data_clim = data_clim.rename({'month': 'time'})
    data_clim.load()
    data_clim *= 30
    #data_clim *= mask

    for pb in pre_box:
        pp2 = data_clim.sel(lon=slice(pre_box_lon[pb][0], pre_box_lon[pb][1]),
                     lat=slice(pre_box_lat[pb][0], pre_box_lat[pb][1]))
        #### k-means #####
        X = pp2.stack(new=['lon', 'lat'])['prec'].T
        X[np.where(np.isnan(X))] = -99

        SilhScore(X, title='CFSv2' + ' - pre Box: ' + str(pb), name_fig='CFSv2_L' + str(l) + '_SilhScore_Box' + str(pb), save=True)

        for n_c in [4, 5, 6, 7, 8]:
            reshape_lat_dim = len(range(pre_box_lat[pb][0], pre_box_lat[pb][1])) + 1
            reshape_lon_dim = len(range(pre_box_lon[pb][0], pre_box_lon[pb][1])) + 1

            clusters, centroids = KM(data=X, n_clusters=n_c, centroids=True,
                                     reshape_dim=[reshape_lon_dim, reshape_lat_dim])

            xr_cluster_label = MakeMask(pp2, dataname='cluster')
            xr_cluster_label['cluster'].values = clusters

            # xr_cluster_label *= MakeMask(data_clim, dataname='cluster')
            xr_cluster_label = xr_cluster_label.interp(
                lat=np.linspace(pre_box_lat[pb][0], pre_box_lat[pb][1], reshape_lat_dim * 8),
                lon=np.linspace(pre_box_lon[pb][0], pre_box_lon[pb][1], reshape_lon_dim * 8),
                method='nearest')

            Plot(comp=xr_cluster_label, levels=np.arange(0, n_c, 1), cmap='tab20',
                 title='CFSv2 Leadtime: ' + str(l) + ' - n_clusters: ' + str(n_c),
                 name_fig='CFSv2_L' + str(l) + '_n_c' + str(n_c) + '_preBox-' + str(pb),
                 dpi=dpi, save=save, OceanMask=True)

            for n_c2 in np.arange(0, n_c, 1):
                PlotCentroids(serie=centroids[n_c2, :], label=str(n_c2),
                              title='CFSv2 Leadtime ' + str(l) + 'Centroids cluster ' + str(n_c2) + '\n' +
                                    '- total Clusters' + str(n_c) + ' - pre Box ' + str(pb),
                              name_fig='CFSv2_L' + str(l) + '_centroid_n_c2' + str(n_c2) + '_n_c'
                                       + str(n_c) + '_preBox-' + str(pb),
                              dpi=100, save=save, ymax=500)

########################################################################################################################