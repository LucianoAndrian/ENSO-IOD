"""
Selecciona los conjuntos de leads de NMME CFSv2 en JJA, JAS, ASO y SON en hindcast y realtime para hgt
"""
import numpy as np
import xarray as xr
import pandas as pd
from ENSO_IOD_Funciones import SelectNMMEFiles

########################################################################################################################
dir_hc = '/pikachu/datos/luciano.andrian/hindcast/'
dir_rt = '/pikachu/datos/luciano.andrian/real_time/'
out_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/by_months/'

# Functions ############################################################################################################
def open_mfdataset_merge_only(paths):
    """
    Los .nc de real_time con xr.open_mfdataset: Resulting object does not have monotonic global indexes along dimension S
    esto es lo mismo pero con mas ram...
    """
    return xr.merge([xr.open_dataset(path,decode_times=False).sel(X=slice(275, 330), Y=slice(-60, 15)) for path in paths])

def fix_calendar(ds, timevar='time'):
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds
########################################################################################################################

variables = ['hgt']
anios = np.arange(1982, 2021)

leads = [0,1,2,3]
#seasons = ['JJA', 'JAS', 'ASO', 'SON']
months = ['Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov']
mmonth_numbers = [6,7,8,9,10,11]


"""
Los realtime estan operativos? si es asi el total de eventos va cambiando...
actualmente estan actualizados.
"""

# current_year = datetime.now().year
# current_month= datetime.now().month
anios_fill = [np.arange(1982, 2011), np.arange(2011, 2021)]

v = 'hgt'
print(v)
for ms in mmonth_numbers:
    print(months[ms - 6])
    r_count = 0
    for r in range(1, 25):  # loop en los miembros de ensamble
        print(r)
        count_path = 0
        for path in [dir_hc, dir_rt]:
            print(path)
            time = []
            for num in anios_fill[count_path]:
                for m in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:
                    if m < 10:
                        time.append(str(num) + '-0' + str(m) + '-01')
                    else:
                        time.append(str(num) + '-' + str(m) + '-01')
            if count_path == 0:
                dates = pd.to_datetime(time, format='%Y/%m/%d')
                anios = np.arange(1982, 2011)
                count_path = 1
            else:
                dates = pd.to_datetime(time, format='%Y/%m/%d')
                anios = np.arange(2011, 2021)

            # abre por r .. "*_r1_*"
            files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                                    dir=path,
                                    by_r=True, r=str(r))

            # ordenando por anios
            files = sorted(files, key=lambda x: x.split()[0])

            data = xr.open_mfdataset(files, decode_times=False)  # se pudre con los tiempos.... ***
            data = data.sel(X=slice(10, 330), Y=slice(-80, 20), P=200)
            if path == dir_rt:
                data = data.drop('Z')

            data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})

            data = data.drop('P')

            data = xr.decode_cf(fix_calendar(data))

            # data['time'] = dates  # reemplazando con las fechas ***


            data['L'] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            data = data.rolling(time=3, center=True).mean()
            """
            3rm, facilita la seleccion y lo hace de manera igual a DMI y N34
            """

            l_count = 0
            print('Forecast for ' + months[ms - 6])
            for l in leads:
                ic_month = ms - l  # 3rm, se toma el mes del centro del trimestre
                if ic_month == 0:
                    ic_month = 12
                elif ic_month == -1:
                    ic_month = 11
                print('issued at ' + str(ic_month))

                data_season = data.sel(time=data.time.dt.month.isin(ic_month), L=l)

                if l_count == 0:
                    data_season_f = data_season
                    l_count = 1
                else:
                    data_season_f = xr.concat([data_season_f, data_season], dim='time')

            if r_count == 0:
                data_season_r_f = data_season_f
                r_count = 1
            else:
                data_season_r_f = xr.merge([data_season_r_f, data_season_f])
                #print(data_season_r_f)

    data_season_r_f.to_netcdf(out_dir + months[ms - 6] + '_' + v + '_Leads_r_CFSv2.nc')
    print('Save' + months[ms - 6])
