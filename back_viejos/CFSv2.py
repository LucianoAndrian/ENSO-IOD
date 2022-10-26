import glob
import numpy as np
import xarray as xr

# Funciones ############################################################################################################
def SelectNMMEFiles(model_name, variable, dir, anio, in_month):
    if ((isinstance(model_name, str) == False) | (isinstance(variable,str)==False) |
            (isinstance(dir, str)==False) | (isinstance(in_month, str)==False)
            | (isinstance(anio, str)==False)):
        print('ERROR: model_name, variable, dir and in_month must be a string')
        return

    if int(in_month) < 10:
        m_in = '0'
    else:
        m_in = ''

    if in_month == '1':
        y1 = 0
        m1 = -11
        m_en = ''
    elif int(in_month) > 10:
        y1 = 1
        m1 = 1
        m_en=''
        print('Year in chagend')
        anio = str(int(anio) - 1)
        print(anio)
    else:
        y1 = 1
        m1 = 1
        m_en = '0'


    files = glob.glob(dir + variable + '_Amon_' + model_name + '_' + anio + m_in + in_month +
                      '_r*_' + anio + m_in + in_month + '-' + str(int(anio) + y1) + m_en +
                      str(int(in_month) - m1) + '.nc')
    # print(dir + variable + '_Amon_' + model_name + '_' + anio + m_in + in_month +
    #                   '_r*_' + anio + m_in + in_month + '-' + str(int(anio) + y1) + m_en +
    #                   str(int(in_month) - m1) + '.nc')

    return files
########################################################################################################################

dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
dir_rt = '/pikachu/datos/osman/nmme/monthly/real_time/'
out_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/'

variables = ['tref', 'prec']
anios = np.arange(1982, 2021)

leads = [0,1,2,3,4,5,6,7]
seasons = ['JJA', 'JAS', 'ASO', 'SON']
mmonth_seasons = [7,8,9,10]


for v in variables:
    print(v)

    for l in leads:
        print('Lead: ' + str(l))

        for m in mmonth_seasons:
            print('Forecast for ' + seasons[m - 7])
            in_month = str(m - l - 1) # VER!
            if in_month == '0':
                in_month = '12'
            elif in_month == '-1':
                in_month = '11'
            print('issued at ' + in_month)
            print('Loading years...')

            check=True
            for y in anios:
                if in_month == '12' or in_month == '11':
                    y += 1

                files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                                        dir=dir_hc, anio=str(y), in_month=in_month)

                for f in files:
                    data = xr.open_dataset(f, decode_times=False)
                    leads_data = data.L
                    r = data.M.values
                    if r <= 24:
                        s = data.S.values
                        data = data.sel(X=slice(275, 330), Y=slice(-60, 15), L=leads_data[l:l+3], S=s[0])
                        data = data.drop(['S'])
                        data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r'})
                        data = data.mean('L')
                        data = data.expand_dims({'Year':[y]})
                        if check:
                            data_final = data
                            check = False
                            print('check first file')
                        else:
                            data_final = xr.merge([data_final, data])

                files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                                        dir=dir_rt, anio=str(y), in_month=in_month)

                for f in files:
                    data = xr.open_dataset(f, decode_times=False)
                    leads_data = data.L
                    r = data.M.values
                    if r <= 24:
                        s = data.S.values
                        data = data.sel(X=slice(275, 330), Y=slice(-60, 15), L=leads_data[l:l+3], S=s[0])
                        data = data.drop(['S'])
                        data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r'})
                        data = data.mean('L')
                        data = data.expand_dims({'Year':[y]})
                        data_final = xr.merge([data_final, data])

            print('save as: ' + v + '_CFSv2_' + seasons[m - 7] + '_Lead_' + str(l) + '.nc')
            data_final.to_netcdf(out_dir + v + '_CFSv2_' + seasons[m - 7] + '_Lead_' + str(l) + '.nc' )
            del data_final
            del data