import xarray as xr
from windspharm.xarray import VectorWind
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
################################################################################
#data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/1940_2020/'

preproc = True
compute = False
RWS = True
SF = False
# Requiere de todo el dominio...
################################################################################
if preproc:
    dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5' \
                '/downloaded/'
    out_dir = data_dir
    # Funciones ################################################################
    def Detrend(xrda, dim):
        aux = xrda.polyfit(dim=dim, deg=1)
        try:
            trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
        except:
            trend = xr.polyval(xrda[dim], aux.polyfit_coefficients[0])
        dt = xrda - trend
        return dt

    def Weights(data):
        weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                       (len(data.lon), 1)))
        data_w = data * weights
        return data_w


    #--------------------------------------------------------------------------#
    name_variables = ['u', 'v']

    for v in ['UV200', 'UV750']:
        for n_v in name_variables:
            data = xr.open_dataset(dir_files + 'ERA5_' + v + '_40-20.nc')
            if n_v == 'u':
                print('Drop v')
                data = data.drop('v')
            elif n_v == 'v':
                print('Drop u')
                data = data.drop('u')

            data = data.rename({n_v: 'var'})
            data = data.rename({'longitude': 'lon'})
            data = data.rename({'latitude': 'lat'})

            data = Weights(data)
            data = data.rolling(time=3, center=True).mean()
            # <<<<<<<<<<<<<<VERIFICAR COMO QUEDA WAF>>>>>>>>>>>>>>>>>>>>>>>>
            for mm, s_name in zip([10], ['SON']):
                aux = data.sel(time=data.time.dt.month.isin(mm))
                aux = Detrend(aux, 'time')
                print('to_netcdf...')
                aux.to_netcdf(out_dir + n_v + '_' + v + '_w_.nc')

                aux2 = aux+data.mean('time')
                aux.to_netcdf(out_dir + n_v + '_' + v + '_w_detrend.nc')
################################################################################
if compute:
    from windspharm.examples import example_data_path

    for v in ['UV200', 'UV750']:
        ds = xr.open_dataset(example_data_path('uwnd_mean.nc'))
        uwnd_aux = ds['uwnd']

        # al haber aplicado el rolling el primer y ultimo
        # mes no estan y falla Vectorwind
        uwnd = xr.open_dataset(data_dir + 'u_' + v + '_w_.nc')
        uwnd = uwnd.sel(time=slice('1940-02-01', '2020-11-01'))
        vwnd = xr.open_dataset(data_dir + 'v_' + v + '_w_.nc')
        vwnd = vwnd.sel(time=slice('1940-02-01', '2020-11-01'))

        uwnd = uwnd.interp(lat=np.linspace(uwnd_aux.latitude.values[0],
                                           uwnd_aux.latitude.values[-1], 179),
                           lon=np.linspace(uwnd_aux.longitude.values[0],
                                           uwnd_aux.longitude.values[-1], 359))
        vwnd = vwnd.interp(lat=np.linspace(uwnd_aux.latitude.values[0],
                                           uwnd_aux.latitude.values[-1], 179),
                           lon=np.linspace(uwnd_aux.longitude.values[0],
                                           uwnd_aux.longitude.values[-1], 359))

        uwnd = uwnd.to_array()
        vwnd = vwnd.to_array()

        w = VectorWind(uwnd, vwnd)

        #### DE A UNO!!! #####
        eta = w.absolutevorticity()
        eta.to_netcdf(data_dir + 'eta_' + v + '.nc')
        div = w.divergence()
        div.to_netcdf(data_dir + 'div_' + v + '.nc')
        del div
        uchi, vchi = w.irrotationalcomponent()
        uchi.to_netcdf(data_dir + 'uchi_' + v + '.nc')
        del uchi
        vchi.to_netcdf(data_dir + 'vchi_' + v + '.nc')
        del vchi
        etax, etay = w.gradient(eta)
        etax.to_netcdf(data_dir + 'etax_' + v + '.nc')
        etay.to_netcdf(data_dir + 'etay_' + v + '.nc')
        del etax, etay
        sf, vp = w.sfvp()
        sf.to_netcdf(data_dir + 'sf_from_' + v + '_w.nc')
        vp.to_netcdf(data_dir + 'vp_from_' + v + '_w.nc')
        del sf, vp

        # sf = w.streamfunction()
        # sf.to_netcdf(data_dir + 'sf_from__w_' + v + '.nc')


################################################################################

if RWS:
    for v in ['UV200', 'UV750']:
        # # Combine the components to form the Rossby wave source term.

        eta = xr.open_dataset(data_dir + 'eta_' + v + '.nc')
        eta = eta.rename({'absolute_vorticity': 'var'})
        div = xr.open_dataset(data_dir + 'div_' + v + '.nc')
        div = div.rename({'divergence': 'var'})
        uchi = xr.open_dataset(data_dir + 'uchi_' + v + '.nc')
        vchi = xr.open_dataset(data_dir + 'vchi_' + v + '.nc')
        uchi = uchi.rename({'u_chi': 'var'})
        vchi = vchi.rename({'v_chi': 'var'})
        etay = xr.open_dataset(data_dir + 'etay_' + v + '.nc')
        etax = xr.open_dataset(data_dir + 'etax_' + v + '.nc')
        etay = etay.rename({'meridional_gradient_of_absolute_vorticity': 'var'})
        etax = etax.rename({'zonal_gradient_of_absolute_vorticity': 'var'})

        S = eta * -1. * div - (uchi * etax + vchi * etay)
        S.to_netcdf(data_dir + 'RWS' + v + '.nc')
    print('END!')

################################################################################
## ah! y Streamfunction y VP ######  ###########################################
################################################################################
#------------------------------------------------------------------------------#

if SF:
    from windspharm.examples import example_data_path

    dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5' \
                '/downloaded/'
    out_dir = data_dir
    # -------------------------------------------------------------------------#
    name_variables = ['u', 'v']
    v_count = 0
    for v in ['UV200', 'UV750']:
        for n_v in name_variables:
            data = xr.open_dataset(dir_files + 'ERA5_' + v + '_40-20.nc')
            #n_v = name_variables[v_count]
            if n_v == 'u':
                print('Drop v')
                data = data.drop('v')
            elif n_v == 'v':
                print('Drop u')
                data = data.drop('u')

            data = data.rename({n_v: 'var'})
            data = data.rename({'longitude': 'lon'})
            data = data.rename({'latitude': 'lat'})

            ds = xr.open_dataset(example_data_path('uwnd_mean.nc'))
            uwnd_aux = ds['uwnd']

            data = data.interp(lat=np.linspace(uwnd_aux.latitude.values[0],
                                               uwnd_aux.latitude.values[-1],
                                               179),
                               lon=np.linspace(uwnd_aux.longitude.values[0],
                                               uwnd_aux.longitude.values[-1],
                                               359))

            data = Weights(data)
            # data = Detrend(data, 'time')

            data.to_netcdf(out_dir + n_v + '_mer_d_w_world.nc')
            v_count += 1
        del data

        uwnd = xr.open_dataset(data_dir + 'u_mer_d_w_world.nc')
        vwnd = xr.open_dataset(data_dir + 'v_mer_d_w_world.nc')

        uwnd = uwnd.to_array()
        vwnd = vwnd.to_array()

        w = VectorWind(uwnd, vwnd)
        #sf, vp = w.sfvp()
        sf = w.streamfunction()
        sf.to_netcdf(data_dir + 'sf_from_' + v + '_w.nc')
        #vp.to_netcdf(data_dir + 'vp_from_w.nc')
################################################################################



























