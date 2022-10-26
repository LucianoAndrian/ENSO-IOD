import time
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
from ENSO_IOD_Funciones import OpenDatasets
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
########################################################################################################################
nc_date_dir = '/datos/luciano.andrian/ncfiles/nc_composites_dates/'
data_dir = '/datos/luciano.andrian/ncfiles/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/w_sig/'
sig_dir = '/datos/luciano.andrian/ncfiles/nc_quantiles/'

start = ('1920', '1950')
seasons = ("Full_Season", 'JJA', 'ASO', 'SON')
min_max_months = [[7,11], [6,8],[8,10],[9,11]]
variables = ['vp', 'div']

seasons = ("Full_Season", 'ASO', 'SON')
min_max_months = [[7,11], [8,10],[9,11]]


# cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
#          'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']


cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_un_pos', 'DMI_un_neg']

scales = [np.linspace(-450, 450, 21),  #hgt
          np.linspace(-3, 3, 13),  #psl
          np.linspace(-4.5e6, 4.5e6, 13),  #sf
          np.linspace(-4.5e6, 4.5e6, 13),  #vp
          np.linspace(-1, 1 ,17),  #t
          np.linspace(-1, 1 ,17),  #t
          np.linspace(-30, 30, 13)] #pp

SA = [False,False,False,False,True, True, True]
step = [1,1,1,1,1,1,1]
text = True
#
# title_case = ['DMI-ENSO simultaneous positive phase ',
#               'DMI-ENSO simultaneous negative phase ',
#               'DMI negative phase ',
#               'DMI positive phase ',
#               'DMI isolated positive phase ',
#               'DMI isolated negative phase ',
#               'ENSO positive phase ',
#               'ENSO negative phase ',
#               'ENSO isolated positive phase ',
#               'ENSO isolated negative phase ']


title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI isolated positive phase ',
              'DMI isolated negative phase ']

v_name = ['HGT 200hPa', 'PSL',
          'Stream Function', 'Potential Velocity',
          'Temperature - Cru', 'Temperature - BEIC', 'Precipitation - GPCC']
# function
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

#
def is_months(month, mmin, mmax):
    return (month >= mmin) & (month <= mmax)
########################################################################################################################
c_var=0
#for v in variables:
v = 'div200'
i = '1950'
data_dir2 = '/datos/luciano.andrian/ncfiles/ERA5/'
data = xr.open_dataset(data_dir2 + v + '.nc')

c_cases = 0
for c in cases:
    print(c)

    count = 0

    for s in seasons:
        #data_sig = xr.open_dataset(sig_dir + v + '_' + c + '_' + i + '_2020' + '_' + s + '.nc')

        aux = xr.open_dataset(nc_date_dir + 'Composite_' + i + '_2020_' + s + '.nc')
        neutro = aux.Neutral
        case = aux[c]
        aux.close()

        mmonth = min_max_months[count]

        mmin = mmonth[0]
        mmax = mmonth[-1]

        neutro_comp = CompositeSimple(original_data=data, index=neutro,
                                      mmin=mmin, mmax=mmax)

        iodw = data.sel(lat=slice(10.0, -10.0), lon=slice(50, 70))#, time=slice('195-01-01', '2020-12-31'))
        iodw = iodw.sel(time=is_months(month=iodw['time.month'], mmin=mmin, mmax=mmax))

        iode = data.sel(lat=slice(10.0, -10.0), lon=slice(90, 110))#, time=slice('1920-01-01', '2020-12-31'))
        iode = iode.sel(time=is_months(month=iode['time.month'], mmin=mmin, mmax=mmax))

        fechas = iode.time

        iodw_neutro = neutro_comp.sel(lat=slice(10.0, -10.0), lon=slice(50, 70))

        iode_neutro = neutro_comp.sel(lat=slice(0, -10.0), lon=slice(90, 110))

        w_anom = iodw - iodw_neutro
        e_anom = iode - iode_neutro

        w_anom = w_anom['var'].mean(['lon', 'lat'], skipna=True)
        e_anom = e_anom['var'].mean(['lon', 'lat'], skipna=True)

        w_anom = np.convolve(w_anom, np.ones((3,)) / 3, mode='same')
        e_anom = np.convolve(e_anom, np.ones((3,)) / 3, mode='same')

        w_anom = xr.DataArray(w_anom, coords=[fechas.time.values], dims=['time'])
        e_anom = xr.DataArray(e_anom, coords=[fechas.time.values], dims=['time'])

        aux = w_anom.sel(time=w_anom.time.dt.year.isin(case.values))
        # aux = aux.sel(
        #     time=is_months(month=aux['time.month'], mmin=mmin, mmax=mmax))

        aux2 = e_anom.sel(time=e_anom.time.dt.year.isin(case.values))
        # aux2 = aux2.sel(
        #     time=is_months(month=aux2['time.month'], mmin=mmin, mmax=mmax))


        fig, ax = plt.subplots()
        im = w_anom.plot(label='West', alpha=0.7, color='orangered')
        e_anom.plot(alpha=0.7, label='East', color='steelblue')
        plt.scatter(x=aux2.time, y=aux2.values, color='steelblue', s=40)
        plt.scatter(x=aux.time, y=aux.values, color='orangered', s=40)
        ax.grid(True)
        plt.ylim(-5e-6, 5e-6)
        plt.legend()
        plt.title(title_case[c_cases] + '\n' + s)
        plt.show()

        # data_comp = CompositeSimple(original_data=data, index=case,
        #                             mmin=mmin, mmax=mmax)
        # comp = data_comp - neutro_comp
        #
        # sig = comp.where((comp < data_sig['var'][0]) | (comp > data_sig['var'][1]))
        # try:
        #     sig = sig.where(np.isnan(sig['var']), 0)
        #     # sig = sig.where(~sig['var']!=0,1)
        # except:
        #     print('Sin valores significativos')
        #
        # print(title_case[c_cases])
        # Plot(comp=comp, comp_sig=sig, significance=True, cmap=cmap[c_var],
        #      levels=scales[c_var], number_events=case,
        #      SA=SA[c_var], step=step[c_var], text=text, dpi=200,
        #      title=v_name[c_var] + '\n' + title_case[c_cases] + '\n' + s + '  ' + i + ' - 2020',
        #      name_fig=v + '_Comp_' + c + '_' + s + '_' + i + '_2020',
        #      save=True)

        count += 1
    c_cases += 1

    # c_var += 1


from ENSO_IOD_Funciones import DMI

# indices: ----------------------------------------------------------------------------------------------------#
dmi = DMI(filter_bwa=False, start_per='1920', end_per='2020')[2]
dmi = dmi.sel(time=slice(str(1920) + '-01-01', str(2020) + '-12-01'))


def FFT(serie, minA = 0, maxA = 40):

    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.fft import fft, fftfreq

    N = np.size(serie)
    T = 1/N

    yf = fft(serie)
    xf = fftfreq(N, T)[:N // 2]

    plt.bar(xf[minA:maxA], 2/N * (np.abs(yf[minA:maxA])**2/2))
    plt.grid()
    plt.show()

