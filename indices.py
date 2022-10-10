
import xarray as xr
import numpy as np
from statistics import stdev
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl
from itertools import groupby

import os
from Funciones import MovingBasePeriodAnomaly
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

w_dir = '/home/luciano.andrian/doc/salidas/'

sst = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")

plt.contourf(sst.sst[1])
nsamples, nx, ny = sst.sst.shape
sst_r = sst.sst.values.reshape((nsamples,nx*ny))



#N34
ninio34 = sst.sel(lat=slice(4.0,-4.0), lon=slice(190, 240), time=slice('1920-01-01', '2020-12-01'))
ninio34 = ninio34.sst.mean(['lon', 'lat'], skipna=True)

ninio3 = sst.sel(lat=slice(4.0,-4.0), lon=slice(210, 270), time=slice('1920-01-01', '2020-12-01'))
ninio3 = ninio3.sst.mean(['lon', 'lat'], skipna=True)

#trend
trend = np.polyfit(range(0,len(ninio34)), ninio34, deg=1)
ninio34 = ninio34 - trend[0]*range(0, len(ninio34))

#compute monthly anomalies
ninio34 = ninio34.groupby('time.month') - ninio34.groupby('time.month').mean('time', skipna=True)

#compute 3-month running mean
#ninio34_filtered = np.convolve(ninio34, np.ones((3,))/3, mode = 'same')
ninio34_filtered = np.convolve(ninio34, [.25,.5,.25], mode='same') # filtro binomial 1,2,1 /4 ??
ninio34_f = xr.DataArray(ninio34_filtered, coords=[ninio34.time.values], dims=['time'])
#ninio34_f.to_netcdf('ninio34.nc4')


seasons = [12,1]
titulos = [0,'DJF',None, 'FMA', 'MAM',None, 'MJJ', 'JJA',None,'ASO', 'SON',None,  'NDJ']

cmap = colors.ListedColormap(['blue', 'paleturquoise', 'white', 'navajowhite', 'red'])
cmap.set_bad(color="white")

sd = np.round(stdev(ninio34_f.values), 2)  # SD solo de la seasons, criterio deser?

for i in seasons:

    aux = ninio34_f.sel(time=ninio34.month == i)
    #sd = np.round(stdev(aux.values), 2)  # SD solo de la seasons, criterio deser?

    en = str(len(aux[np.where(aux > .5)]))
    en_st = str(len(aux[np.where(aux > sd)]))
    print(en)
    ln = str(len(aux[np.where(aux < -.5)]))
    ln_st = str(len(aux[np.where(aux < -sd)]))

    if i == 1:
        aux = np.concatenate((aux[1:101], np.full((10,), 0)))
        aux = aux.reshape((11, 10))
        aux = np.around(aux, 2)
    else:
        aux = np.concatenate((aux, np.full((9,), 0)))
        aux = aux.reshape((11, 10))
        aux = np.around(aux, 2)

    decadas = [f'+ {num}' for num in range(0, 10, 1)]
    anios = [f'{num}' for num in range(1920, 2030, 10)]

    aux2 = np.ma.masked_where(np.abs(aux) < 0.5, aux) #mi criterio
    bounds = [-2,-sd,0,sd,2]
    cmap = colors.ListedColormap(['dodgerblue', 'paleturquoise', 'navajowhite', 'firebrick'])
    cmap.set_bad(color="white")

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots()
    im = plt.imshow(aux2, cmap=cmap, norm=norm)#norm=norm, cmap=cmap,interpolation='nearest',vmax=2 + 1)
    #plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds)
    #plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)
    ax.set_xticks(np.arange(len(decadas)))
    ax.set_yticks(np.arange(len(anios)))
    ax.set_xticklabels(decadas)
    ax.set_yticklabels(anios)


    plt.title(titulos[i] + '  -  SD serie'+  str(sd), fontsize = 12)

    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")
    for K in range(len(anios)):
        for j in range(len(decadas)):
            text = ax.text(j, K, aux[K, j],
                           ha="center", va="center", color="black")
    ax.xaxis.set_ticks_position('top')
    #fig.tight_layout()

    plt.figtext(0.5, 0.01,
                en + 'EN & ' + ln + 'LN |T|>0.5;  '+ en_st + 'EN & ' + ln_st + ' LN |T|>SD'
                , ha = 'center', va='bottom',  bbox={"facecolor":"green", "alpha":0.3, "pad":5},
                fontsize = 12)
    plt.savefig(w_dir + 'n34_deser_sdtotal' + titulos[i] + '.jpg', figsize=(5, 5), dpi=200)
    plt.show()
    plt.close()

########################################################################################################################

#trend
trend = np.polyfit(range(0,len(ninio3)), ninio3, deg=1)
ninio3 = ninio3 - trend[0]*range(0, len(ninio3))

#compute monthly anomalies
ninio3 = ninio3.groupby('time.month') - ninio3.groupby('time.month').mean('time', skipna=True)

#compute 3-month running mean
#ninio34_filtered = np.convolve(ninio34, np.ones((3,))/3, mode = 'same')
ninio3_filtered = np.convolve(ninio3, [.25,.5,.25], mode='same') # filtro binomial 1,2,1 /4 ??
ninio3_f = xr.DataArray(ninio3_filtered, coords=[ninio3.time.values], dims=['time'])
#ninio34_f.to_netcdf('ninio34.nc4')


seasons = [12,1]
titulos = [0,'DJF',None, 'FMA', 'MAM',None, 'MJJ', 'JJA',None,'ASO', 'SON',None,  'NDJ']

cmap = colors.ListedColormap(['blue', 'paleturquoise', 'white', 'navajowhite', 'red'])
cmap.set_bad(color="white")


for i in seasons:

    aux = ninio3_f.sel(time=ninio3.month == i)
    sd = np.round(stdev(aux.values), 2)  # SD solo de la seasons, criterio deser?

    en = str(len(aux[np.where(aux > .5)]))
    en_st = str(len(aux[np.where(aux > sd)]))
    print(en)
    ln = str(len(aux[np.where(aux < -.5)]))
    ln_st = str(len(aux[np.where(aux < -sd)]))

    if i == 1:
        aux = np.concatenate((aux[1:101], np.full((10,), 0)))
        aux = aux.reshape((11, 10))
        aux = np.around(aux, 2)
    else:
        aux = np.concatenate((aux, np.full((9,), 0)))
        aux = aux.reshape((11, 10))
        aux = np.around(aux, 2)

    decadas = [f'+ {num}' for num in range(0, 10, 1)]
    anios = [f'{num}' for num in range(1920, 2030, 10)]

    aux2 = np.ma.masked_where(np.abs(aux) < 0.5, aux) #mi criterio
    bounds = [-2,-sd,0,sd,2]
    cmap = colors.ListedColormap(['dodgerblue', 'paleturquoise', 'navajowhite', 'firebrick'])
    cmap.set_bad(color="white")

    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    fig, ax = plt.subplots()
    im = plt.imshow(aux2, cmap=cmap, norm=norm)#norm=norm, cmap=cmap,interpolation='nearest',vmax=2 + 1)
    #plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds)
    #plt.colorbar(im, cmap=cmap, norm=norm, boundaries=bounds, ticks=bounds)
    ax.set_xticks(np.arange(len(decadas)))
    ax.set_yticks(np.arange(len(anios)))
    ax.set_xticklabels(decadas)
    ax.set_yticklabels(anios)


    plt.title(titulos[i] + '  -  SD'+  str(sd), fontsize = 12)

    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")
    for K in range(len(anios)):
        for j in range(len(decadas)):
            text = ax.text(j, K, aux[K, j],
                           ha="center", va="center", color="black")
    ax.xaxis.set_ticks_position('top')
    #fig.tight_layout()

    plt.figtext(0.5, 0.01,
                en + 'EN & ' + ln + 'LN |T|>0.5;  '+ en_st + 'EN & ' + ln_st + ' LN |T|>SD'
                , ha = 'center', va='bottom',  bbox={"facecolor":"green", "alpha":0.3, "pad":5},
                fontsize = 12)
    plt.savefig(w_dir + 'n3_deser' + titulos[i] + '.jpg', figsize=(5, 5), dpi=200)
    plt.show()
    plt.close()
########################################################################################################################

# CPC criteria. |sst'| > 0.5 for a minimum of 5 consecutive overlapping seasons


def Nino34CPC(data):

    from Funciones import MovingBasePeriodAnomaly

    sst = data
    # N34
    ninio34 = sst.sel(lat=slice(4.0, -4.0), lon=slice(190, 240), time=slice('1906-01-01', '2020-12-31'))
    ninio34 = ninio34.sst.mean(['lon', 'lat'], skipna=True)

    # compute monthly anomalies
    ninio34 = MovingBasePeriodAnomaly(ninio34)

    # compute 5-month running mean
    ninio34_filtered = np.convolve(ninio34, np.ones((3,)) / 3, mode='same')  #
    ninio34_f = xr.DataArray(ninio34_filtered, coords=[ninio34.time.values], dims=['time'])

    aux = abs(ninio34_f) > 0.5
    results = []
    for k, g in groupby(enumerate(aux.values), key=lambda x: x[1]):
        if k:
            g = list(g)
            results.append([g[0][0], len(g)])

    n34 = []
    for m in range(0, len(results)):
        # True values
        len_true = results[m][1]

        # True values for at least 5 consecutive seasons
        if len_true >= 5:
            a = results[m][0]
            n34.append([np.arange(a, a + results[m][1]), ninio34_f[np.arange(a, a + results[m][1])].values])

    return ninio34_f, n34

n34 = Nino34CPC(sst)[1]

aux = np.full_like(ninio34_f, 0)
for i in range(0, len(n34)):
    for m in range(0,len(n34[i][0])):
        aux[n34[i][0][m]] = n34[i][1][m]


##### graficar ####
bounds = [-2, 0, 2]
cmap = colors.ListedColormap(['dodgerblue', 'firebrick'])
cmap.set_bad(color="white")
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
mes = ['DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA', 'JAS', 'ASO', 'SON', 'OND', 'NDJ']

for v in range(5):
    aux1 = ninio34_f.values
    aux2 = ninio34_f.values

    if v == 4:
        anios = [f'{num}' for num in range(2000, 2021, 1)]
        start = 80
        end = 101
    else:
        anios = [f'{num}' for num in range(1920 + (20 * v), 1940 + (20 * v), 1)]
        start = 0 + (20 * v)
        end = 20 + (20 * v)

    # colors mask
    aux2 = np.ma.masked_where(aux2 != aux, aux)
    #aux2 = np.ma.masked_where(abs(aux2) < .4, aux2)

    aux1 = aux1.reshape((101, 12))
    aux2 = aux2.reshape((101, 12))

    aux1 = np.around(aux1, 1)

    aux1 = aux1[start:end, ]
    aux2 = aux2[start:end, ]

    fig, ax = plt.subplots()
    im = plt.imshow(aux2, cmap=cmap, norm=norm)


    # major ticks
    ax.set_xticks(np.arange(len(mes)))
    ax.set_yticks(np.arange(len(anios)))
    ax.set_xticklabels(mes)
    ax.set_yticklabels(anios)

    # minor ticks
    ax.set_yticks(np.arange(len(anios)) + 0.5, minor=True)
    ax.set_yticklabels([], minor=True)
    ax.tick_params(axis='y', which="minor", length=1)

    ax.set_xticks(np.arange(len(mes)) + 0.5, minor=True)
    ax.set_xticklabels([], minor=True)
    ax.tick_params(axis='x', which="minor", length=0)

    ax.grid(which='minor', alpha=1)
    ax.grid(which='major', alpha=0)
    plt.grid(True)

    fig.set_size_inches(8, 9)

    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")
    for K in range(len(anios)):
        for j in range(len(mes)):
            text = ax.text(j, K, aux1[K, j],
                           ha="center", va="center", color="black")
    ax.xaxis.set_ticks_position('top')
    # fig.tight_layout()

    plt.title('Niño3.4 |SST*| > 0.5 for a minimum of 5 consecutive overlapping seasons')
#    plt.savefig(w_dir + 'N34_CPC_' + anios[0] + "_" + anios[len(anios) - 1] + '.jpg', dpi=200)
    plt.show()
    plt.close()












#
# # NOAA20cr hasta 2014
# psl_noaa20 = xr.open_dataset("/pikachu/datos3/reanalysis/NCEP20CR/prmsl.mon.mean.nc")
# psl_noaa20 = psl_noaa20.sel(time=slice('1920-01-01', '2014-12-01'))
# psl_noaa20 = psl_noaa20.rename({'prmsl':'psl'})
#
# sam_noaa20 = SamDetrend(psl_noaa20)
#
# sam_n20_seasonal = xr.DataArray(np.convolve(sam_noaa20, np.ones((3,))/3, mode='same')
#                                 , coords=[sam_noaa20.time.values], dims=['time'])
#
# sam_n20_anual = xr.DataArray(np.convolve(sam_noaa20, np.ones((12,))/12, mode='same')
#                          , coords=[sam_noaa20.time.values], dims=['time'])
#
#
# #Series
# ####
# #ERA20c vs ERA5
# # 1979-2010
# era20c = psl_era20c.sel(time=slice('1979-01-01', '2010-12-01'))
# era5 = psl_era5.sel(time=slice('1979-01-01', '2010-12-01'))
# from scipy import stats
#
# stats.ttest_ind([1,1,1,1], [0.5,1,1,1])
#
# data=era20c
# p40 = data.sel(lat=-40, method='nearest').mean(dim='lon')/100
# trend = np.polyfit(range(0, len(p40.psl)), p40.psl, deg=1)
# p40_era20 = p40.psl - trend[0] * range(0, len(p40.psl))
# p65 = data.sel(lat=-65, method='nearest').mean(dim='lon')/100
# trend = np.polyfit(range(0, len(p65.psl)), p65.psl, deg=1)
# p65_era20 = p65.psl - trend[0] * range(0, len(p65.psl))
#
# data=era5
# p40 = data.sel(lat=-40, method='nearest').mean(dim='lon')/100
# trend = np.polyfit(range(0, len(p40.psl)), p40.psl, deg=1)
# p40_era5 = p40.psl - trend[0] * range(0, len(p40.psl))
# p65 = data.sel(lat=-65, method='nearest').mean(dim='lon')/100
# trend = np.polyfit(range(0, len(p65.psl)), p65.psl, deg=1)
# p65_era5 = p65.psl - trend[0] * range(0, len(p65.psl))
#
#
#
# p40_era20.plot(label='ERA20C')
# p40_era5.plot(label='ERA5')
#
# era20_mean = np.around(p40_era20.mean('time').values,2)
# era20_sd = np.around(p40_era20.std('time').values,2)
#
# era5_mean = np.around(p40_era5.mean('time').values,2)
# era5_sd = np.around(p40_era5.std('time').values,2)
#
# plt.title('40ºS Mean y SD [hPa] \n ERA20c: ' + str(era20_mean) +
#           ' ' + str(era20_sd) + '\n ERA5: ' + str(era5_mean) +
#           ' ' + str(era5_sd))
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()
#
# p65_era20.plot(label='ERA20C')
# p65_era5.plot(label='ERA5')
#
# era20_mean = np.around(p65_era20.mean('time').values,2)
# era20_sd = np.around(p65_era20.std('time').values,2)
#
# era5_mean = np.around(p65_era5.mean('time').values,2)
# era5_sd = np.around(p65_era5.std('time').values,2)
#
# plt.title('65ºS Mean y SD [hPa] \n ERA20c: ' + str(era20_mean) +
#           ' ' + str(era20_sd) + '\n ERA5: ' + str(era5_mean) +
#           ' ' + str(era5_sd))
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()
#
# import matplotlib.dates as mdates
#
# def PlotSAM (data1, data2, titulo, nombre, start=0, label1="", label2=""):
#     fig, ax = plt.subplots()
#     im = data1.plot(color='orangered', label=label1)
#     data2.plot(color='royalblue', label=label2)
#     plt.legend()
#     ax.grid(True)
#     fig.set_size_inches(9, 4)
#     plt.tight_layout()
#     ax.set_yticks(range(-5, 5))
#     ax.xaxis_date()
#     plt.axhline(y=0, color='black', linestyle='-')
#     datemin = np.datetime64(sam_n20_seasonal.time.values[start] - 1, 'Y')
#     datemax = np.datetime64(sam_n20_seasonal.time.values[-1] + 1, 'Y') + np.timedelta64(1, 'Y')
#     ax.set_xlim(datemin, datemax)
#     ax.format_xdata = mdates.DateFormatter('%Y-%m')
#     plt.title(titulo)
#     plt.tight_layout()
#     plt.savefig(w_dir + nombre + '.jpg', dpi=200)
#     plt.show()
#     plt.close()
#
# PlotSAM(data1=sam_n20_seasonal, data2=sam_era20c_seasonal,
#         titulo='SAM - 3-Month Running Mean', nombre='SAM_mon_n20-e20',
#         label1='NOAA20CR', label2="ERA-20C")
#
# PlotSAM(data1=sam_n20_anual, data2=sam_era20c_anual,
#         titulo='SAM - 12-Month Running Mean', nombre='SAM_anu_n20-e20',
#         label1='NOAA20CR', label2="ERA-20C")
#
# PlotSAM(data1=sam_era20c_seasonal, data2=sam_era5_seasonal,
#         titulo='SAM 1979-2020 - 3-Month Running Mean', nombre='SAM_mon_e20-e5',
#         start=708,
#         label1='ERA-20C', label2="ERA5")
#
# PlotSAM(data1=sam_era20c_anual, data2=sam_era5_anual,
#         titulo='SAM 1979-2020 - 12-Month Running Mean', nombre='SAM_anu_e20-e5',
#         start=708,
#         label1='ERA-20C', label2="ERA5")
#
#
# def PlotSAMSeasonal (data, titulo, nombre):
#     seasons = [0, 'DJF', None, 'FMA', 'MAM', None, 'MJJ', 'JJA', None, 'ASO', 'SON', None, 'NDJ']
#
#     for m in [3,4,6,7,9,10,12,1]:
#         crit = 0
#         z1 = np.array(data.groupby('time.month')[m].values)
#         z2 = np.array([crit] * len(z1))
#
#         plt.figure(figsize=(8, 4))
#         plt.plot(data.groupby('time.month')[m].time.values,
#                  data.groupby('time.month')[m].values)
#         data.groupby('time.month')[m].plot(color='black')
#
#         plt.fill_between(data.groupby('time.month')[m].time.values,
#                          data.groupby('time.month')[m].values,
#                          where=(z1 < z2),
#                          color='royalblue', interpolate=True)
#
#         plt.fill_between(data.groupby('time.month')[m].time.values,
#                          data.groupby('time.month')[m].values,
#                          where=(z1 > z2),
#                          color='red', interpolate=True)
#
#
#         plt.grid()
#         plt.ylim((-5, 5))
#
#         pos = len(np.where(data.groupby('time.month')[m] > 0)[0])
#         neg = len(np.where(data.groupby('time.month')[m] < 0)[0])
#         plt.title('SAM - ' + titulo + '  ' + seasons[m] + '  +' + str(pos) + ' -' + str(neg))
#         plt.savefig(w_dir + nombre + '_' + seasons[m] + '.jpg', dpi=200)
#         plt.tight_layout()
#         plt.show()
#         plt.close()
#
#
# PlotSAMSeasonal(sam_n20_seasonal, titulo="NOAA 20CR", nombre='n20')
# PlotSAMSeasonal(sam_era20c_seasonal, titulo="ERA-20C", nombre='era20c')
# PlotSAMSeasonal(sam_era5_seasonal, titulo="ERA5", nombre='era5')
#
# data=sam_era20c_seasonal
#
#         crit = 0
#         z1 = np.array(data.groupby('time.month')[4,7,10,1].values)
#         z2 = np.array([crit] * len(z1))
#
#         plt.figure(figsize=(8, 4))
#         plt.plot(data.groupby('time.month')[4,7,10,1].time.values,
#                  data.groupby('time.month')[4,7,10,1].values)
#         data.groupby('time.month')[m].plot(color='black')
#
#         plt.fill_between(data.groupby('time.month')[4,7,10,1].time.values,
#                          data.groupby('time.month')[4,7,10,1].values,
#                          where=(z1 < z2),
#                          color='royalblue', interpolate=True)
#
#         plt.fill_between(data.groupby('time.month')[4,7,10,1].time.values,
#                          data.groupby('time.month')[4,7,10,1].values,
#                          where=(z1 > z2),
#                          color='red', interpolate=True)
#
#         plt.grid()
#         plt.ylim((-5, 5))
#
#         pos = len(np.where(data.groupby('time.month')[[4,7,10,1] > 0)[0])
#         neg = len(np.where(data.groupby('time.month')[4,7,10,1] < 0)[0])
#         #plt.title('SAM - ' + titulo + '  ' + seasons[4,7,10,1] + '  +' + str(pos) + ' -' + str(neg))
#         #plt.savefig(w_dir + nombre + '_' + seasons[4,7,10,1] + '.jpg', dpi=200)
#         plt.tight_layout()
#         plt.show()
#         plt.close()
#
#
#
# data = sam_era20c_seasonal
#
# def prueba(month):
#     return((month == 1) | (month == 4) | (month == 7) | (month == 10))
#
# aux = data.sel(time=prueba(data['time.month']))
# #aux = aux[236::]
# crit = 0
# z1 = np.array(aux.values)
# z2 = np.array([crit] * len(z1))
#
#
# fig, ax = plt.subplots()
# im = plt.plot(aux.time.values,
#          aux.values, color = 'black', alpha=1)
#
# plt.figure(figsize=(8, 4))
#
#
# plt.fill_between(aux.time.values,
#                  aux.values,
#                  where=(z1 < z2),
#                  color='royalblue', interpolate=True)
#
# plt.fill_between(aux.time.values,
#                  aux.values,
#                  where=(z1 > z2),
#                  color='red', interpolate=True)
# ax.grid(True)
# fig.set_size_inches(9, 4)
#
# ax.set_yticks(range(-5, 5))
# ax.xaxis_date()
# plt.axhline(y=0, color='black', linestyle='-')
# datemin = np.datetime64(aux.time.values[0] - 1, 'Y')
# datemax = np.datetime64(aux.time.values[-1] + 1, 'Y') + np.timedelta64(1, 'Y')
# ax.set_xlim(datemin, datemax)
# ax.format_xdata = mdates.DateFormatter('%Y-%m')
#
# plt.grid()
# plt.ylim((-5, 5))
#
# #pos = len(np.where(data.groupby('time.month')[m] > 0)[0])
# #neg = len(np.where(data.groupby('time.month')[m] < 0)[0])
# plt.title('SAM Seasonal - ERA20C')
# #plt.savefig(w_dir + nombre + '_' + seasons[m] + '.jpg', dpi=200)
# plt.tight_layout()
# plt.show()
# plt.close()
#
#
#
#
# x =open('/home/luciano.andrian/doc/scrips/newsam.1957.2007.txt','r')
