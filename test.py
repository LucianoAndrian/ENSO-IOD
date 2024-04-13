import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
data = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
data = data.sel(time=slice('1920-01-01', '2020-12-01'), lat=slice(20,-20), lon=slice(50,270) )
#data = data.groupby('time.month') - data.groupby('time.month').mean('time')

### Prueba detrend 1. xr.polyfit y polival ###
#ajusta a un polinomio de grado 1 para cada punto de reticula
aux = data.polyfit(dim='time', deg=1) # devuelve pendiente y ordenada al origen

# usando la pendiente e intercep se crea la recta de tendencia para cada punto de reticula
trend = xr.polyval(data['time'], aux.sst_polyfit_coefficients[0])
#filtra la tendencia restando al campo la recta en cada punto de reticula
detrend = data - trend
### Prueba detrend 2,. scipy signal ###
#####
plt.plot(data.mean(['lon','lat']).sst, label='original_xr')
plt.plot(detrend.mean(['lon','lat']).sst, label='detrend_xr')
plt.legend()
plt.show()
