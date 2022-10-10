"""
Probando EOF para delimitar regiones
"""
########################################################################################################################
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from eofs.xarray import Eof
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
#from ENSO_IOD_Funciones import SelectNMMEFiles
########################################################################################################################
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Regiones/'
plot=False
save=False
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
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' + 'pp_gpcc_d_w_c_1950-2020_1.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        #pp_gpcc = pp_gpcc.rename({'precip': 'var'})

        return pp_gpcc


def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt

    comp_var = comp
    fig = plt.figure(figsize=(5, 6), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270,330, -60,20], crs_latlon)

    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step, ::step],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
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

def PlotCp(serie, serie2, dpi, title, name_fig, save):
    fig = plt.figure(figsize=(6, 3), dpi=dpi)
    ax = fig.add_subplot(111)
    ax.plot(serie, linewidth=2, color='dodgerblue', label='Obs')
    ax.plot(serie2, linewidth=2, color='indianred', label='CFSv2')
    ax.tick_params(labelsize=7)
    ax.set_xlabel('Months')
    ax.set_xticks(range(0,12), ['J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'], size=10)
    plt.title(title, fontsize=10)
    plt.legend()
    plt.grid()
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

########################################################################################################################
# Observado -----------------------------------------------------------------------------------------------------------#
pp = OpenDataSet('pp_gpcc', interp=True,
                 lat_interp=np.linspace(-60,15,76),
                 lon_interp=np.linspace(275,330,56))

# Mismo período que el CFSv2 (hindcast)
pp = pp.sel(time=pp.time.dt.year.isin(range(1982,2012)))

# Climatologia
pp = pp.groupby('time.month').mean()
pp = pp.rename({'month':'time'})

# EOF
# solver = Eof(pp['var'], center=False)
# eof_obs = solver.eofs(neofs=2)
# cps_obs = solver.pcs(npcs=6)
# eof_cor_obs = solver.eofsAsCorrelation(neofs=6)
# var_obs = np.around(solver.varianceFraction(neigs=6).values*100,1)
####
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

# mascara para el oceano ----------------------------------------------------------------------------------------------#
mask = OpenDataSet('pp_gpcc', interp=True,
                   lat_interp=np.linspace(-60, 15, 76),
                   lon_interp=np.linspace(275, 330, 56))
mask = mask.mean('time')
mask = xr.where(np.isnan(mask), mask, 1)
mask = mask.rename({'var': 'prec'})

color=[plt.cm.tab20b(0),plt.cm.tab20b(1),plt.cm.tab20b(2),plt.cm.tab20b(3),
       plt.cm.tab20c(0),plt.cm.tab20c(1),plt.cm.tab20c(2),plt.cm.tab20c(3),
       '#8CEFDE','white','white', '#F9FD7C',
       plt.cm.tab20c(7), plt.cm.tab20c(6), plt.cm.tab20c(5), plt.cm.tab20c(4),
       plt.cm.tab20b(15), plt.cm.tab20b(14), plt.cm.tab20b(13), plt.cm.tab20b(12)]
cmap=colors.ListedColormap(color)

# Por leadtimes -------------------------------------------------------------------------------------------------------#
for l in [0,1,2,3]:
    data2 = data.sel(L=l)
    data_clim = data2.groupby('time.month').mean()
    data_clim = data_clim.mean(['r'])
    data_clim = data_clim.rename({'month': 'time'})
    data_clim.load()
    data_clim *= 30
    data_clim *= mask

    # # Eof -------------------------------------------------------------------------------------------------------------#
    # solver = Eof(data_clim['prec'], center=False)
    # #eof_cfsv2 = solver.eofs(neofs=6)
    # cps_cfsv2 = solver.pcs(npcs=5)
    # eof_cor_cfsv2 = solver.eofsAsCorrelation(neofs=5)
    # #var_cfsv2 = np.around(solver.varianceFraction(neigs=6).values * 100, 1)
    #
    # if l==0:
    #     cps_cfsv2_l = cps_cfsv2
    #     cps_corr_l = eof_cor_cfsv2
    # else:
    #     # concatenando por leadtime
    #     # CPs aún desplazadas con los leadtimes!
    #     cps_cfsv2_l = xr.concat([cps_cfsv2_l, cps_cfsv2], dim='L')
    #     cps_corr_l = xr.concat([cps_corr_l, eof_cor_cfsv2], dim='L')
    #
    # # Plot ------------------------------------------------------------------------------------------------------------#
    # if plot:
    #     for i in range(0, 4):
    #         PlotCp(serie=cps_obs.sel(mode=i),
    #                serie2=np.hstack([cps_cfsv2.sel(mode=i)[-l::], cps_cfsv2.sel(mode=i)[0:-l]]),
    #                title='CP ' + str(i + 1) + ' Leadtime: ' + str(l),
    #                name_fig='series_CP_' + str(i + 1) + '_leadtime' + str(l),
    #                dpi=200, save=save)
    #
    #         Plot(comp=eof_cor_obs.sel(mode=i), cmap=cmap, levels=np.linspace(-1, 1, 21),
    #              title='Obs - Correlation between CP ' + str(i + 1) + ' and variable',
    #              name_fig='corr_obs_CP' + str(i + 1),  # se sobreescribe... que genialidad.
    #              dpi=300, save=save)
    #
    #         Plot(comp=eof_cor_cfsv2.sel(mode=i), cmap=cmap, levels=np.linspace(-1, 1, 21),
    #              title='CFSv2 - Correlation between CP ' + str(i + 1) + ' and variable - Leadtime: ' + str(l),
    #              name_fig='corr_CFSv2_CP' + str(i + 1) + '_leadtime' + str(l),
    #              dpi=300, save=save)
# ---------------------------------------------------------------------------------------------------------------------#
# # Invirtiendo cps y corr
# # CP 3 (mode=2) y leadtime = 3
# for l in [1,2,3]:
#     for cp in [0,1,2,3,4]:
#         cps_cfsv2_l.loc[dict(mode=cp, L=l)] = \
#             np.hstack([cps_cfsv2_l.sel(mode=cp, L=l)[-l::], cps_cfsv2_l.sel(mode=cp, L=l)[0:-l]])
#
#
# cps_cfsv2_l.loc[dict(mode=2, L=3)] *= -1
# cps_corr_l.loc[dict(mode=2, L=3)] *= -1
#
# # CP 4 (mode=3) y Leatime = 2
# cps_corr_l.loc[dict(mode=3, L=2)] *= -1
# cps_cfsv2_l.loc[dict(mode=3, L=2)] *= -1
#
# #CP 5 (mode=4) y Leadtime = 0, 2 y 3
# cps_cfsv2_l.loc[dict(mode=4, L=0)] *= -1
# cps_cfsv2_l.loc[dict(mode=4, L=2)] *= -1
# cps_cfsv2_l.loc[dict(mode=4, L=3)] *= -1
#
# cps_corr_l.loc[dict(mode=1)] *= -1
#
# mode=1
# plt.plot(cps_cfsv2_l.sel(mode=mode, L=0), label='Lead = 0')
# plt.plot(cps_cfsv2_l.sel(mode=mode, L=1), label='Lead = 1')
# plt.plot(cps_cfsv2_l.sel(mode=mode, L=2), label='Lead = 2')
# plt.plot(cps_cfsv2_l.sel(mode=mode, L=3), label='Lead = 3')
# plt.title('CP ' + str(mode + 1))
# plt.legend()
# plt.show()
#
#
# cps_corr_f = cps_corr_l.mean(['L'])

# for rc in [0.7,0.6,0.5,0.4]:
#     for i in range(0, 5):
#         aux = cps_corr_f.sel(mode=i).where(cps_corr_f.sel(mode=i).__abs__() > rc)
#         Plot(comp=aux, cmap=cmap, levels=np.linspace(-1, 1, 21),
#              title='CFSv2 - Correlation between CP ' + str(i + 1) + ' and variable' +
#                    '\n' + r"$\bf{" + '|r|>' + str(rc) + "}$" + ' ' +
#                    r"$\bf{" + ' (R2>' + str(np.around(rc**2,2)) + ')' + "}$",
#              name_fig='corr_CFSv2_CP' + str(i + 1) + '_rc_' + str(rc),
#              dpi=100, save=save)
#
# # ---------------------------------------------------------------------------------------------------------------------#
#eofAsCorrelation verification:
# es necesario no usar mask

def corr(x,y):
    from scipy.stats import pearsonr
    return pearsonr(x,y)[0]

def CorrCP(x, y, dim='time'):
    return xr.apply_ufunc(
        corr, x, y,
        input_core_dims=[[dim], [dim]],
        vectorize=True,  # !Important!
        output_dtypes=[float])

#
# aux = CorrCP(cps_cfsv2.sel(mode=0), data_clim) - eof_cor_cfsv2.sel(mode=0)
# plt.imshow(aux['prec']);plt.colorbar();plt.show()
# print(aux.mean(['lon', 'lat']).prec) # ~e-9
# #max y min ~e-7
# ---------------------------------------------------------------------------------------------------------------------#


k = 8

X = pp.stack(new=['lon','lat'])['var'].T
X = data_clim.stack(new=['lon', 'lat'])['prec'].T
X[np.where(np.isnan(X))]=-99

from threadpoolctl import threadpool_limits
with threadpool_limits(limits=1): #con 1 mas rapido q con 40...
    kmeans = KMeans(n_clusters=k, n_init=100).fit(X)

y_pred = kmeans.predict(X)

xr_cluster_label = mask.copy()
xr_cluster_label['prec'].values = y_pred.reshape(56,76).T
xr_cluster_label *= mask
xr_cluster_label = xr_cluster_label.rename({'prec':'cluster'})



Plot(comp=xr_cluster_label, levels=np.arange(0,9,1), save=False, cmap='Set1', dpi=300)
def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         dpi=100, save=True, step=1,
         name_fig='fig', title='title'):

    import matplotlib.pyplot as plt

    comp_var = comp
    fig = plt.figure(figsize=(5, 6), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270,330, -60,20], crs_latlon)

    im = ax.contourf(comp.lon, comp.lat, comp_var['cluster'],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')
    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey')
    ax.add_feature(cartopy.feature.COASTLINE)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
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
