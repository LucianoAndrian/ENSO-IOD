"""
Calculo del DMI y Nino3.4 en cada miembro de ensamble y lead del NMME-CFSv2
"""
import matplotlib.pyplot as plt
import pandas as pd
import xarray as xr
import numpy as np
from ENSO_IOD_Funciones import SelectNMMEFiles


# Funciones ############################################################################################################

def TwoClimEnsemble(x, output=1):
    #clim1 1982-1998
    aux1 = x.sel(time=x.time.dt.year.isin(np.linspace(1982, 1998, 17)))
    #aux1 = aux1.mean(['r'])

    # clim1 1999-2011
    aux2 = x.sel(time=x.time.dt.year.isin(np.linspace(1999, 2011, 13)))
    #aux2 = aux2.mean(['r'])
    if output==1:
        return aux1.mean(['r']), aux2.mean(['r'])
    else:
        return aux1, aux2


def CheckTrend(x):
    from scipy.stats import kendalltau
    x['L']=[0,1,2,3]
    for m in [7, 8, 9, 10]:
        l_count = 0
        for l in [0, 1, 2, 3]:
            aux4 = x.sel(time=x.time.dt.month.isin(m - l), L=l)
            p = kendalltau(aux4.sst, range(0, len(aux4.sst)))[1]
            if p < 0.05:
                print('Tendencia en ' + str(m) + ' con lead ' + str(l))


def DetrendTwoClimEnsemble(x,x1,x2, path):

    if path == dir_hc:
        aux1 = x.sel(time=x.time.dt.year.isin(np.linspace(1982, 1998, 17)))
        aux1_anom = aux1.groupby('time.month') - x1.groupby('time.month').mean()
        aux1_anom = aux1_anom.mean(['r'])

        aux2 = x.sel(time=x.time.dt.year.isin(np.linspace(1999, 2011, 13)))
        aux2_anom = aux2.groupby('time.month') - x2.groupby('time.month').mean()
        aux2_anom = aux2_anom.mean(['r'])

        aux3 = xr.concat([aux1_anom, aux2_anom], dim='time')
        #aux4_trend = np.polyfit(range(0, len(aux4.sst)), aux4.sst, deg=1)


    elif path == dir_rt:
        print('Real_time')
        aux1 = x
        aux1_anom = aux1.groupby('time.month') - x1.groupby('time.month').mean()
        aux1_anom = aux1_anom.mean(['r'])
        aux3 = aux1_anom
        #x_trend = np.polyfit(range(0, len(aux3.sst)), aux3.sst, deg=1)

    return aux3

def PreProc(x, clim1, clim2, path, detrend=False, trend=None):
    # x = np.apply_along_axis(lambda m: np.convolve(m, np.ones((3,)) / 3, mode='same'), axis=0,
    #                                  arr=x.sst)
    x = x.rolling(time=3, center=True).mean()

    if path==dir_hc:
        # clim1
        x1 = x.sel(time=x.time.dt.year.isin(np.linspace(1982, 1998, 17)))
        x_anom_clim1 = x1.groupby('time.month') - clim1.groupby('time.month').mean()

        # clim2
        x2 = x.sel(time=x.time.dt.year.isin(np.linspace(1999, 2011, 13)))
        x_anom_clim2 = x2.groupby('time.month') - clim2.groupby('time.month').mean()

        x3 = xr.concat([x_anom_clim1, x_anom_clim2], dim='time')

    elif path==dir_rt:
        x1 = x
        x_anom_clim2 = x1.groupby('time.month') - clim2.groupby('time.month').mean()

        x3=x_anom_clim2

    # if detrend:
    #     x_detrend = x3.sst[:, :, 0] - \
    #                 (trend[0] * np.repeat([np.arange(1, len(x3.sst) + 1)], 4, axis=0).T + trend[1])
    #
    #     x_filtered = np.apply_along_axis(lambda m: np.convolve(m, np.ones((3,)) / 3, mode='same'), axis=0,
    #                                      arr=x_detrend)
    # else:
    #     x_filtered = np.apply_along_axis(lambda m: np.convolve(m, np.ones((3,)) / 3, mode='same'), axis=0,
    #                                      arr=x3.sst)

    return x3

def fix_calendar(ds, timevar='time'):
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds
########################################################################################################################
dir_hc = '/pikachu/datos/luciano.andrian/hindcast/'
dir_rt = '/pikachu/datos/luciano.andrian/real_time/'
out_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/'

leads = [0, 1, 2, 3]
seasons = ['JJA', 'JAS', 'ASO', 'SON']
mmonth_seasons = [7, 8, 9, 10]

v = 'sst'

anios_fill = [np.arange(1982, 2012), np.arange(2011, 2021)]

# Climatologias separadas ##############################################################################################
#### HINDCAST ####
files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                        dir=dir_hc, All=True)
files = sorted(files, key=lambda x: x.split()[0])

data = xr.open_mfdataset(files, decode_times=False)
data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
data = data.sel(L=[0.5, 1.5, 2.5, 3.5])
data = xr.decode_cf(fix_calendar(data))

iodw = data.sel(lon=slice(50, 70), lat=slice(-10, 10)).mean(['lon', 'lat'])
iode = data.sel(lon=slice(90, 110), lat=slice(-10, 0)).mean(['lon', 'lat'])
n34 = data.sel(lat=slice(-4, 4), lon=slice(190, 240)).mean(['lon', 'lat'])

#rolling
iodw = iodw.rolling(time = 3, center = True).mean()
iode = iode.rolling(time = 3, center = True).mean()
n34 = n34.rolling(time = 3, center = True).mean()

#Tomando climatologia y tendencia en el caso de IODW y IODE
# n34 no significativa
iodw_clim1, iodw_clim2 = TwoClimEnsemble(iodw)
iodw_clim1 = iodw_clim1.load()
iodw_clim2 = iodw_clim2.load()
#iodw_trend_hc = DetrendTwoClimEnsemble(iodw, iodw_clim1, iodw_clim2, dir_hc)
#print("#################################################################")
#print("iodw")
#CheckTrend(iodw_trend_hc)

iode_clim1, iode_clim2 = TwoClimEnsemble(iode)
iode_clim1 = iode_clim1.load()
iode_clim2 = iode_clim2.load()
# iode_trend_hc = DetrendTwoClimEnsemble(iode, iode_clim1, iode_clim2, dir_hc)
# print("#################################################################")
# print("iode")
# CheckTrend(iode_trend_hc)

n34_clim1, n34_clim2 = TwoClimEnsemble(n34)
n34_clim1 = n34_clim1.load()
n34_clim2 = n34_clim2.load()

#iode_trend_hc = DetrendTwoClimEnsemble(n34, n34_clim1, n34_clim2, dir_hc)
# print("#################################################################")
# print("n34")
# CheckTrend(iode_trend_hc)

#### REAL_TIME ####
#usa *clim2 del hindcast (1999-2011)
# sólo es necesario tomar la tendencia de IODW y IODE para filtrarla
# n34 no significativa

# files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
#                         dir=dir_rt, All=True)
# files = sorted(files, key=lambda x: x.split()[0])
#
# data = xr.open_mfdataset(files, decode_times=False)
# data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
# data = data.sel(L=[0.5, 1.5, 2.5, 3.5])
# data = xr.decode_cf(fix_calendar(data))
#
# iodw = data.sel(lon=slice(50, 70), lat=slice(-10, 10)).mean(['lon', 'lat'])
# iode = data.sel(lon=slice(90, 110), lat=slice(-10, 0)).mean(['lon', 'lat'])
#
# iodw_trend_rt = DetrendTwoClimEnsemble(iodw, iodw_clim1, iodw_clim2, dir_rt)
# iode_trend_rt = DetrendTwoClimEnsemble(iode, iodw_clim1, iodw_clim2, dir_rt)

#detrend_iods=[False, False]
# iodw_trends = [iodw_trend_hc, None]
# iode_trends = [iode_trend_hc, None]
########################################################################################################################

for ms in mmonth_seasons: # loop en las estaciones (mes central)
    print(seasons[ms-7])
    r_count = 0
    for r in range(1, 25):  # loop en los miembros de ensamble
        print(r)
        count_path = 0
        for path in [dir_hc, dir_rt]: # loop hindcast y real_time
            # Fechas
            time = []
            for num in anios_fill[count_path]:
                for m in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
                    if m < 10:
                        time.append(str(num) + '-0' + str(m) + '-01')
                    else:
                        time.append(str(num) + '-' + str(m) + '-01')
            if count_path == 0:
                dates = pd.to_datetime(time[0:-9], format='%Y/%m/%d')
                anios = np.arange(1982, 2011)
            else:
                dates = pd.to_datetime(time[3:], format='%Y/%m/%d')
                anios = np.arange(2011, 2021)

            # Abriendo archivos nc por r ###############################################################################
            files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                                    dir=path,
                                    by_r=True, r=str(r))
            # ordenando por anios
            files = sorted(files, key=lambda x: x.split()[0])

            data = xr.open_mfdataset(files, decode_times=False)  # se pudre con los tiempos.... ***
            data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
            data = data.sel(L=[0.5, 1.5, 2.5, 3.5])
            data = xr.decode_cf(fix_calendar(data)) #***

            # DMI y N34 ################################################################################################
            iodw = data.sel(lon=slice(50, 70), lat=slice(-10, 10)).mean(['lon', 'lat'])
            iode = data.sel(lon=slice(90, 110), lat=slice(-10, 0)).mean(['lon', 'lat'])
            n34 = data.sel(lat=slice(-4, 4), lon=slice(190, 240)).mean(['lon', 'lat'])

            # ---------------------- control NANs (a ojo)  ------------------------ #
            # print('#############' + str(r) + '###############' )
            # for L in [0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5]:
            #     for t in range(0,351):
            #         if np.isnan(n34.load().sel(L=L).sst[t].values):
            #             print(L)
            #             print(n34.time[t].values)

            # PrePROC ##################################################################################################
            """
            PreProc toma la anomalia mensual en funcion de las climatologias.
            No hay tendencia en las series de cada serie para la media del ensamble -> no filtro
            """
            iodw_f = PreProc(iodw, iodw_clim1, iodw_clim2, path,
                             detrend=False)
            iode_f = PreProc(iode, iode_clim1, iode_clim2, path,
                             detrend=False)
            n34_f = PreProc(n34, n34_clim1, n34_clim2, path,
                            detrend=False)
            dmi = iodw_f - iode_f

            # if detrend_iods[count_path]==False:
            #     dmi = dmi[:,:,0]

            # dmi = xr.Dataset({'index': (('time', 'L'), dmi)},
            #                  coords={'time': iodw.time.values, 'L': [0, 1, 2, 3]})

            # n34 = xr.Dataset({'index': (('time', 'L'), n34_f[:,:,0])},
            #                  coords={'time': iodw.time.values, 'L': [0, 1, 2, 3]})

            count_path = 1
            # Selección estación #######################################################################################
            dmi['L'] = [0, 1, 2 ,3]
            dmi = dmi.sel(r=r)
            dmi = dmi.drop('r')

            n34 = n34_f
            n34['L'] = [0, 1, 2, 3]
            n34 = n34.sel(r=r)
            n34 = n34.drop('r')

            l_count = 0
            for l in leads:

                ic_month = ms - l  # tomar el mes del centro del trimestre debido al 3rm
                if ic_month == 0:  # JJA y Lead 7
                    ic_month = 12

                dmi_season = dmi.sel(time=dmi.time.dt.month.isin(ic_month), L=l)
                n34_season = n34.sel(time=n34.time.dt.month.isin(ic_month), L=l)

                if l_count == 0:
                    dmi_season_f = dmi_season
                    n34_season_f = n34_season
                    l_count = 1
                else:
                    dmi_season_f = xr.concat([dmi_season_f, dmi_season], dim='time')
                    n34_season_f = xr.concat([n34_season_f, n34_season], dim='time')

            dmi_season_f = dmi_season_f.expand_dims({'r': [r]})
            n34_season_f = n34_season_f.expand_dims({'r': [r]})

            if r_count == 0:
                dmi_season_r_f = dmi_season_f
                n34_season_r_f = n34_season_f
                r_count = 1
            else:
                dmi_season_r_f = xr.merge([dmi_season_r_f, dmi_season_f])
                n34_season_r_f = xr.merge([n34_season_r_f, n34_season_f])



    dmi_season_r_f.to_netcdf(out_dir + seasons[ms - 7] + '_DMI_Leads_r_CFSv2.nc')
    n34_season_r_f.to_netcdf(out_dir + seasons[ms - 7] + '_N34_Leads_r_CFSv2.nc')
    print('Save' + seasons[ms - 7])
###