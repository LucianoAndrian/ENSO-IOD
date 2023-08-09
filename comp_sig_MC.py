import xarray as xr
import numpy as np
from multiprocessing.pool import ThreadPool
import os
import glob
import math
from datetime import datetime
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
# Functions ####################################################################
def CompositeSimple(original_data, index, mmin, mmax):
    def is_months(month, mmin, mmax):
        return (month >= mmin) & (month <= mmax)

    if len(index) != 0:
        comp_field = original_data.sel(time=original_data.time.dt.year.isin([index]))
        comp_field = comp_field.sel(
            time=is_months(month=comp_field['time.month'], mmin=mmin, mmax=mmax))
        if len(comp_field.time) != 0:
            comp_field = comp_field.mean(['time'], skipna=True)
        else:  # si sólo hay un año
            comp_field = comp_field.drop_dims(['time'])

        return comp_field
    else:
        print(' len index = 0')

def NumberPerts(data_to_concat, neutro, num = 0):
    #si las permutaciones posibles:
    # > 10000 --> permutaciones = 10000
    # < 1000 ---> permutaciones = todas las posibles
    #
    # si num != 0 permutaciones = num

    total = len(data_to_concat) + len(neutro)
    len1 = len(neutro)
    len2 = len(data_to_concat)

    total_perts = math.factorial(total) / (math.factorial(len2) * math.factorial(len1))

    if num == 0:

        if total_perts >= 10000:
            tot = 10000
            print('M = 10000')
        else:
            tot = total_perts
            print('M = ' + str(total_perts))

    else:
        tot = num

    jump = 9 #10 por .nc que guarde
    M = []
    n = 0

    while n < tot:
        aux = list(np.linspace((0 + n), (n + jump), (jump + 1)))
        M.append(aux)
        n = n + jump + 1

    return M
################################################################################
#data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'
################################################################################
seasons = ['JJA', 'SON']
min_max_months = [[6, 8],[9, 11]]

# seasons = ['JJA']
# min_max_months = [[7, 8]]
print('#######################################################################')
print(seasons[0])
periodos = [['40', '20']]


cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_un_pos', 'DMI_un_neg', 'N34_un_pos',
         'N34_un_neg']
################################################################################
for dmi_true_dipole in [True, False]:
    if dmi_true_dipole:
        nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates/'
        out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_quantiles/DMI_true_dipole/'
    else:
        nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/'
        out_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_quantiles/DMIbase/'

    for v in ['HGT200', 'HGT750', 'tcru', 'ppgpcc']:

        print('Variable: ' + v + ', dmi_true_dipole=' + str(dmi_true_dipole))

        # control variables!!! (algunas ya no estan disponibles)
        if v == 'HGT200':
            ruta = data_dir
            data = xr.open_dataset(ruta + v + '_SON_mer_d_w.nc')
        elif v == 'HGT750':
            ruta = data_dir
            data = xr.open_dataset(ruta + v + '_SON_mer_d_w.nc')
        elif v == 'vp':
            ruta = data_dir
            data = xr.open_dataset(ruta + v + '_from_w.nc')
        elif v == 'div200':
            ruta = data_dir
            data = xr.open_dataset(ruta + v + '_HS_mer_d_w.nc')

        # T y PP PROVISIONAL!
        # Para hacer JJA y SON será a mano (literal) luego solo SON
        elif v == 'tcru':
            ruta = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/'
            data = xr.open_dataset(ruta + 'tcru_w_c_d_0.25_SON.nc')
        elif v == 'ppgpcc':
            ruta = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/'
            data = xr.open_dataset(ruta + 'ppgpcc_w_c_d_1_SON.nc')

        else:
            print('que es esto <<' + v + '>> ?!!!')
            from sys import exit
            exit(1)

        if v != 'ppgpcc':
            data = data.interp(lat=np.arange(-80, 20, 0.5)[::-1],
                               lon=np.arange(0, 360, 0.5))

        for c in cases:
            print(c)
            print('----------------------------------------------------------')
            count = 0
            for s in seasons:
                print(s)

                # si la proxima tenga menos de 10000 permutaciones
                # no se sobreescribirian todas
                files = glob.glob('/pikachu/datos/luciano.andrian/observado/ncfiles/nc_comps/' + '*.nc')
                if len(files) != 0:
                    for f in files:
                        try:
                            os.remove(f)
                        except:
                            print('Error: ' + f)

                mmonth = min_max_months[count]


                def PermuDatesComposite(n, data=data, mmonth=mmonth):
                    mmin = mmonth[0]
                    mmax = mmonth[-1]
                    rn = np.random.RandomState(616)

                    for a in n:
                        dates_rn = rn.permutation(neutro_concat)
                        neutro_new = dates_rn[0:len(neutro)]
                        data_new = dates_rn[len(neutro):]

                        neutro_comp = CompositeSimple(original_data=data,
                                                      index=neutro_new,
                                                      mmin=mmin, mmax=mmax)
                        data_comp = CompositeSimple(original_data=data,
                                                    index=data_new,
                                                    mmin=mmin, mmax=mmax)

                        if (len(data_comp) != 0):
                            if a == n[0]:
                                comp = data_comp - neutro_comp
                                comp = comp.expand_dims(time=[a])
                                comp_concat = comp
                            else:
                                comp = data_comp - neutro_comp
                                comp = comp.expand_dims(time=[a])
                                try:
                                    comp_concat = xr.concat([comp_concat, comp],
                                                            dim='time')
                                except:
                                    if a != n[0]:
                                        comp_concat = comp
                        else:
                            next

                    comp_concat.to_netcdf('/pikachu/datos/luciano.andrian/observado/ncfiles/nc_comps/' + 'Comps_' +
                                          str(int(a)) + '.nc')
                    del comp
                    del data_comp
                    del neutro_new
                    del comp_concat


                # Fechas de los eventos IODS y Ninios detectados a partir de ERSSTv5 en 1920-2020
                aux = xr.open_dataset(nc_date_dir + '1920_2020_' + s + '.nc')

                neutro = aux.Neutral
                data_to_concat = aux[c]
                aux.close()

                M = NumberPerts(data_to_concat, neutro, 0)

                if (data_to_concat[0] != 0):
                    neutro_concat = np.concatenate([neutro, data_to_concat])

                    hour = datetime.now().hour
                    if (hour > 19) | (hour < 8):
                        n_thread = 30
                        pool = ThreadPool(30)
                    else:
                        n_thread = 10
                        pool = ThreadPool(10)

                    print('Threads: ', n_thread)

                    pool.map(PermuDatesComposite, [n for n in M])
                    pool.close()

                    # NO chunks! acumula RAM en cada ciclo. ~14gb en 3 ciclos...
                    aux = xr.open_mfdataset('/pikachu/datos/luciano.andrian/observado/ncfiles/nc_comps/Comps_*.nc',
                                            parallel=True
                                            ,  # chunks={'time': -1},  # 'lat':147,'lon':240},
                                            combine='nested', concat_dim="time", coords="different",
                                            compat="broadcast_equals")
                    # aux = aux['var'].astype(np.float32)
                    print('quantiles')
                    aux = aux.chunk({'time': -1})
                    qt = aux.quantile([.05, .95], dim='time', interpolation='linear')
                    qt.to_netcdf(out_dir + v + '_' + c + '1940_2020_' + s + '.nc', compute=True)

                    # if name == 'ERA20':
                    #     qt.to_netcdf('/datos/luciano.andrian/ncfiles/nc_quantiles/' + v + '_' + c + '_19' + i[0] +
                    #                  '_19' + i[1] + '_' + s + '.nc', compute=True)
                    # else:
                    #     qt.to_netcdf('/datos/luciano.andrian/ncfiles/nc_quantiles/' + v + '_' + c + '_19' + i[0] +
                    #                  '_20' + i[1] + '_' + s + '.nc', compute=True)

                    aux.close()
                    del qt
                else:
                    print('no ' + c)

            count += 1
    print('done dmi_true_dipole=' + str(dmi_true_dipole))
print('done')
########################################################################################################################
########################################################################################################################
# Para setear otras variables y periodos
########################################################################################################################
# for i in periodos:
#     # print('Período: ' + i[0] + '-' + i[1])
#     # if i[0] == '20':
#     #     ruta = data_dir1
#     #     name = 'ERA20'
#     # elif i[0] == '50':
#     #     ruta = data_dir2
#     #     name = 'ERA5'
#     # else:
#     #     print('TA MAL!')
# if v == 'hgt200':
#     ruta = data_dir
#     data = xr.open_dataset(ruta + v + '_HS_mer_d_w.nc')
#     data = data.interp(lat=np.arange(-80, 20)[::-1])
# else:
#     ruta = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/'
#     if v=='pp_prec':
#         data = xr.open_dataset(ruta + v + '_d_w_c_1950-2020_2.5.nc')
#     else:
#         data = xr.open_dataset(ruta + v + '_d_w_c_1950-2020_1.nc')
