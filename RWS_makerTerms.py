import xarray as xr
from windspharm.xarray import VectorWind
import numpy as np
########################################################################################################################
data_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'

preproc = False
compute = False
RWS = False
SF = True
# Requiere de todo el dominio...
########################################################################################################################
if preproc:

    dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/merged/'
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
        weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180), (len(data.lon), 1)))
        data_w = data * weights
        return data_w


    #---------------------------------------------------------------------------#
    name_variables = ['u', 'v']
    v_count=0

    for v in ['UV750', 'UV750']:
        n_v = name_variables[v_count]
        data = xr.open_dataset(dir_files + 'ERA5_' + v + '_50-20_mer.nc')

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
        #data = data.sel(lon=slice(30, 340), lat=slice(20, -80))
        #data = Detrend(data, 'time')

        # ---------------------------------------------------------------------------#
        print('to_netcdf...')

        data.to_netcdf(out_dir + n_v + '_mer_d_w_world.nc')


        v_count += 1

########################################################################################################################
if compute:
    from windspharm.examples import example_data_path

    ds = xr.open_dataset(example_data_path('uwnd_mean.nc'))
    uwnd_aux = ds['uwnd']

    uwnd = xr.open_dataset(data_dir + 'u_mer_d_w_world.nc')
    vwnd = xr.open_dataset(data_dir + 'v_mer_d_w_world.nc')

    uwnd = uwnd.interp(lat=np.linspace(uwnd_aux.latitude.values[0], uwnd_aux.latitude.values[-1], 179),
                       lon=np.linspace(uwnd_aux.longitude.values[0], uwnd_aux.longitude.values[-1], 359))
    vwnd = vwnd.interp(lat=np.linspace(uwnd_aux.latitude.values[0], uwnd_aux.latitude.values[-1], 179),
                       lon=np.linspace(uwnd_aux.longitude.values[0], uwnd_aux.longitude.values[-1], 359))

    uwnd = uwnd.to_array()
    vwnd = vwnd.to_array()

    # Create a VectorWind instance to handle the computations.
    w = VectorWind(uwnd, vwnd)

    #### DE A UNO!!! #####
    # Compute components of rossby wave source: absolute vorticity, divergence,
    # irrotational (divergent) wind components, gradients of absolute vorticity.
    eta = w.absolutevorticity()
    eta.to_netcdf(data_dir + 'eta.nc')

    div = w.divergence()
    div.to_netcdf(data_dir + 'div_from_w.nc')

    uchi, vchi = w.irrotationalcomponent()
    uchi.to_netcdf(data_dir + 'uchi.nc')
    vchi.to_netcdf(data_dir + 'vchi.nc')

    etax, etay = w.gradient(eta)
    etax.to_netcdf(data_dir + 'etax.nc')
    etay.to_netcdf(data_dir + 'etay.nc')

    sf, vp = w.sfvp()
    sf.to_netcdf(data_dir + 'sf.nc')
    vp.to_netcdf(data_dir + 'vp_from_w.nc')
########################################################################################################################

if RWS:
    # # Combine the components to form the Rossby wave source term.

    eta = xr.open_dataset(data_dir + 'eta.nc')
    eta = eta.rename({'absolute_vorticity': 'var'})
    div = xr.open_dataset(data_dir + 'div_from_w.nc')
    div = div.rename({'divergence': 'var'})
    uchi = xr.open_dataset(data_dir + 'uchi.nc')
    vchi = xr.open_dataset(data_dir + 'vchi.nc')
    uchi = uchi.rename({'u_chi': 'var'})
    vchi = vchi.rename({'v_chi': 'var'})
    etay = xr.open_dataset(data_dir + 'etay.nc')
    etax = xr.open_dataset(data_dir + 'etax.nc')
    etay = etay.rename({'meridional_gradient_of_absolute_vorticity': 'var'})
    etax = etax.rename({'zonal_gradient_of_absolute_vorticity': 'var'})

    S = eta * -1. * div - (uchi * etax + vchi * etay)
    S.to_netcdf(data_dir + 'RWS.nc')
    print('END!')

########################################################################################################################
## ah! y Streamfunction y VP ######  ###################################################################################
########################################################################################################################

# Create a VectorWind instance to handle the computation of streamfunction and
# velocity potential.
if SF:
    dir_files = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/merged/'
    out_dir = data_dir
    # Funciones ################################################################
    def Detrend(xrda, dim):
        aux = xrda.polyfit(dim=dim, deg=1)
        try:
            trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
        except:
            trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
        dt = xrda - trend
        return dt

    def Weights(data):
        weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180), (len(data.lon), 1)))
        data_w = data * weights
        return data_w
    # ---------------------------------------------------------------------------#
    name_variables = ['u', 'v']
    v_count = 0

    for v in ['UV750', 'UV750']:
        n_v = name_variables[v_count]
        data = xr.open_dataset(dir_files + 'ERA5_' + v + '_50-20_mer.nc')

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
        # data = data.sel(lon=slice(30, 340), lat=slice(20, -80))
        data = Detrend(data, 'time')

        data.to_netcdf(out_dir + n_v + '_mer_d_w_world_dtr.nc')
        v_count +=1
    from windspharm.examples import example_data_path
    ds = xr.open_dataset(example_data_path('uwnd_mean.nc'))
    uwnd_aux = ds['uwnd']

    uwnd = xr.open_dataset(data_dir + 'u_mer_d_w_world.nc')
    vwnd = xr.open_dataset(data_dir + 'v_mer_d_w_world.nc')

    uwnd = uwnd.interp(lat=np.linspace(uwnd_aux.latitude.values[0], uwnd_aux.latitude.values[-1], 179),
                       lon=np.linspace(uwnd_aux.longitude.values[0], uwnd_aux.longitude.values[-1], 359))
    vwnd = vwnd.interp(lat=np.linspace(uwnd_aux.latitude.values[0], uwnd_aux.latitude.values[-1], 179),
                       lon=np.linspace(uwnd_aux.longitude.values[0], uwnd_aux.longitude.values[-1], 359))

    uwnd = uwnd.to_array()
    vwnd = vwnd.to_array()
    w = VectorWind(uwnd, vwnd)
    sf, vp = w.sfvp()
    sf.to_netcdf(data_dir + 'sf_750.nc')
    #vp.to_netcdf(data_dir + 'vp_from_w.nc')
    ########################################################################################################################



























