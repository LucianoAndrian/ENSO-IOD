"""
Calcular quintiles de los outputs de ~tesis/MonteCarlo_ANOVA.R , testeo_anova_mc.R y MonteCarlo_ACC.R
"""
########################################################################################################################
import xarray as xr
#----------------------------------------------------------------------------------------------------------------------#
dir = '/pikachu/datos/luciano.andrian/aux_nmme_quantiles/'
subdir = ['diff_Temp', 'diff_Prec']
term = 'PPt'
variables = ['temp', 'pp']
#----------------------------------------------------------------------------------------------------------------------#
for i in [0,1]:
    t = term
    path = dir + subdir[i] + '/diff_' + t + '*'
    print(path)
    aux = xr.open_mfdataset(path, combine='nested', concat_dim="time",
                            coords="different", compat="broadcast_equals",
                            parallel=True)

    aux = aux.chunk({'time': -1})
    qt = aux.quantile([.05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95], dim='time', interpolation='linear')

    qt.to_netcdf(dir + 'diff_qt_' + t + '_' + variables[i] + '.nc', compute=True)
    del aux, qt

########################################################################################################################
dir = '/pikachu/datos/luciano.andrian/ACC_aux_nmme_quantiles/'

path = dir + 'diff*'
aux = xr.open_mfdataset(path, combine='nested', concat_dim="time",
                        coords="different", compat="broadcast_equals",
                        parallel=True)
aux = aux.chunk({'time': -1})
qt = aux.quantile([.05, .1, .2, .3, .4, .5, .6, .7, .8, .9, .95], dim='time', interpolation='linear')
qt.to_netcdf(dir + 'diff_ACC_qt.nc', compute=True)