############################################ IOD ####################################

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib import colors
import matplotlib as mpl
import pandas as pd
pd.options.mode.chained_assignment = None
import statsmodels.formula.api as sm
import os
from Funciones import WaveFilter
#from Funciones import FFT2
from Funciones import DMIndex
from Funciones import MovingBasePeriodAnomaly
from Funciones import Nino34CPC
from scipy.stats import pearsonr

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'


w_dir = '/home/luciano.andrian/doc/salidas/DMI/'
ERSST = [True]#, False] # var ?¿?

western_io = slice(50,70)
start_per = '1920'
end_per = '2020'

inch_y = 6
inch_x = 6

baseline = [0,2]

filter_harmonic = True
filter_all_harmonic = True
harmonics = []

filter_bwa = True
####################################### Post processing and Plot Functions #######################################

def SYpostprocessing(dmi, full_season = False, season='JJA'):
    # only for sy2003
    # from June to November
    if full_season:
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['DMI'],
                        'Años': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['Años'],
                        'Mes': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['Mes']})
        aux = aux.groupby(['Años']).mean()
        season = 'JJASON'
    elif season == 'MJJ':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([5, 6, 7]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([5, 6, 7]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([5, 6, 7]))['Mes']})

        aux = aux.groupby(['Años']).mean()

    elif season == 'JJA':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([6, 7, 8]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([6, 7, 8]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([6, 7, 8]))['Mes']})

        aux = aux.groupby(['Años']).mean()

    elif season == 'ASO':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([8, 9, 10]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([8, 9, 10]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([8, 9, 10]))['Mes']})

        aux = aux.groupby(['Años']).mean()


    elif season == 'SON':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([9, 10, 11]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([9, 10, 11]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([9, 10, 11]))['Mes']})

        aux = aux.groupby(['Años']).mean()

    return(aux)

def PlotDMI(dmi, start_year, end_year, start_year2, end_year2,
            titulo, nombre, not_sy=True, inch_y=inch_y, inch_x = inch_x):

    positivos = str(sum(dmi['DMI'] > 0))
    negativos = str(sum(dmi['DMI'] < 0))
    total = str(len(dmi))

    anios = np.arange(int(start_year), int(end_year) + 1)
    aux = pd.DataFrame(columns=['Anios', 'DMI'], dtype=float)
    aux.Anios = anios

    if not_sy:
        dmi = dmi.groupby(['Años']).prod() #no hace nada

    # como hacer esto mas facil??
    for i in range(len(aux)):
        for j in range(len(dmi)):
            if aux.Anios[i] == dmi.index.values[j]:
                aux.DMI[i] = np.sign(dmi.DMI.values[j])

    # ----------------------------------------------------------------------------------#
    anios1 = np.arange(int(start_year), int(end_year) + 10)
    decades = np.around(len(anios1) / 10)

    aux = np.concatenate((aux.DMI, np.full((9,), 0)))
    aux = aux.reshape((int(decades), 10))
    aux = np.around(aux, 2)

    aux1 = anios1.reshape((int(decades), 10))
    aux1 = np.around(aux1, 2)

    # colors
    aux2 = np.ma.masked_where(np.isnan(aux), aux)
    aux2 = np.ma.masked_where(aux2 == 0, aux2)
    bounds = [-2, 0, 2]
    cmap = colors.ListedColormap(['dodgerblue', 'paleturquoise', 'navajowhite', 'firebrick'])
    cmap.set_bad(color="white")
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    # ----------------------------------------------------------------------------------#

    # figure
    fig, ax = plt.subplots()
    im = plt.imshow(aux2, cmap=cmap, norm=norm)

    ax.set_yticks(np.arange(int(decades) - 1) + 0.5, minor=True)
    ax.tick_params(axis='y', which="minor", length=1)
    ax.set_xticks(np.arange(10) + 0.5, minor=True)
    ax.tick_params(axis='x', which="minor", length=0)
    ax.set_xticklabels('')
    ax.set_yticklabels('')
    ax.grid(which='minor', alpha=1)
    ax.grid(which='major', alpha=0)
    plt.grid(True)

    fig.set_size_inches(inch_y, inch_x)

    plt.title(titulo + ' - ' + start_year + '-' + end_year, fontsize=12)

    plt.setp(ax.get_xticklabels(), rotation=0, ha="right",
             rotation_mode="anchor")
    for K in range(int(decades)):
        for j in range(10):
            text = ax.text(j, K, aux1[K, j],
                           ha="center", va="center", color="black")
    ax.xaxis.set_ticks_position('top')

    plt.figtext(0.5, 0.01,
                'Positivos: ' + positivos + '  ' + 'Negativos: ' + negativos + '\n Total: ' + total + '\n Climatological Baseline: ' +
                start_year2 + '-' + end_year2
                , ha='center', va='bottom', bbox={"facecolor": "green", "alpha": 0.3, "pad": 5},
                fontsize=12)

    #plt.savefig(w_dir + 'DMI_' + nombre + '.jpg', dpi=200)
    plt.show()
    plt.close()

def PlotSerieIOD(iodw, iode, timeref, dmi, start=0, end = 120*8, var=True, nombre = "fig"):
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import xarray as xr
    import pandas as pd

    dmi1 = iodw - iode
    dmi1 = xr.DataArray(dmi1, coords=[timeref.time.values], dims=['time'])

    dates = pd.DataFrame(columns=['Fechas'], dtype=float)
    for l in range(len(dmi)):
        if var:
            aux = str(int(dmi.Años[l])) + '-' + str(
                int(dmi.Mes[l])) + '-01T00:00:00.000000000'
            dates = dates.append({'Fechas': aux}, ignore_index=True)
            label = "ERRSSTv5"
        else:
            aux = str(int(dmi.Años[l])) + '-' + str(
                int(dmi.Mes[l])) + '-16T00:00:00.000000000'
            dates = dates.append({'Fechas': aux}, ignore_index=True)
            label = "HadISSTv1.1"

    prueba = pd.DataFrame(columns=['Fechas', 'DMI'])
    for i in range(len(dmi1)):
        for l in range(len(dates)):
            if dmi1.time[i] == pd.to_datetime(dates.Fechas[l]):
                prueba = prueba.append({'Fechas': dmi1.time[i].values, 'DMI': dmi1.values[i]}, ignore_index=True)


    fig, ax = plt.subplots()
    im = dmi1.plot(label=label)
    ax.grid(True)
    plt.scatter(x=prueba.Fechas, y=prueba.DMI, color='red')
    fig.set_size_inches(8, 6)
    #plt.tight_layout()
    plt.ylim((-1, 1))
    ax.xaxis_date()
    plt.axhline(y=0, color='black', linestyle='-')
    plt.title('DMI_TS_'+ label + '_' + nombre)
    datemin = np.datetime64(dmi1.time.values[start], 'Y')
    datemax = np.datetime64(dmi1.time.values[-end], 'Y') + np.timedelta64(1, 'Y')
    ax.set_xlim(datemin, datemax)
    ax.format_xdata = mdates.DateFormatter('%Y-%m')
    #plt.savefig(w_dir + 'DMI_TS_'+ label + '_' + nombre + '.jpg', dpi=200)
    plt.legend()
    #plt.show()
    plt.close()

def CompositeIOD(dmi, original_data, full_season=False, season='JJA'
                 , title='title', nombre='Fig', start='1920', end='2020', save=False):

    if full_season:
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['Mes']})
        aux = aux.groupby(['Años']).mean()
        mmin, mmax = 6, 11

        season = 'JJASON'
    elif season == 'MJJ':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([6]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([6]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([6]))['Mes']})

        aux = aux.groupby(['Años']).mean()
        mmin, mmax = 5, 7

    elif season == 'JJA':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([7]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([7]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([7]))['Mes']})

        aux = aux.groupby(['Años']).mean()
        mmin, mmax = 6, 8

    elif season == 'ASO':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([9]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([9]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([9]))['Mes']})

        aux = aux.groupby(['Años']).mean()
        mmin, mmax = 8, 10

    elif season == 'SON':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([10]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([10]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([10]))['Mes']})

        aux = aux.groupby(['Años']).mean()
        mmin, mmax = 9, 11


    dmi = aux

    def is_months(month):
        return (month >= mmin) & (month <= mmax)

    dmi_positive = dmi.index.values[np.where(dmi['DMI'] > 0)]
    dmi_negative = dmi.index.values[np.where(dmi['DMI'] < 0)]

    sst_original = original_data.sel(lat=slice(15, -15), lon=slice(30, 290),
                                     time=slice(start +'-01-01', end +'-12-31'))
    x = np.arange(int(start), int(end), 1)

    mask = np.in1d(x, dmi.index.values, invert=True)
    neutro = sst_original.sel(time=sst_original.time.dt.year.isin(x[mask]))
    neutro = neutro.sel(time=is_months(neutro['time.month']))
    neutro = neutro.sst.mean(['time'], skipna=True)

    for sign in [0, 1]:
        if (sign == 0) & (len(dmi_positive) != 0):
            dmi_sign = original_data.sel(time=original_data.time.dt.year.isin([dmi_positive]),
                                         lat=slice(15, -15), lon=slice(30, 290))

            dmi_sign = dmi_sign.sel(time=is_months(dmi_sign['time.month']))
            dmi_sign = dmi_sign.sst.mean(['time'], skipna=True)
            titulo2 = '+'
            phase = 'POS'
        elif (sign == 1) & (len(dmi_negative) != 0):
            dmi_sign = original_data.sel(time=original_data.time.dt.year.isin([dmi_negative]),
                                         lat=slice(15, -15), lon=slice(30, 290))

            dmi_sign = dmi_sign.sel(time=is_months(dmi_sign['time.month']))
            dmi_sign = dmi_sign.sst.mean(['time'], skipna=True)
            titulo2 = '-'
            phase = 'NEG'
        else:
            continue

        comp = dmi_sign - neutro


        cmap = plt.cm.get_cmap("RdBu")
        cmap.set_bad(color='gray') # esto no funca...

        fig, ax = plt.subplots(figsize=(8, 3), dpi=200)
        im = comp.plot.contourf(cbar=cmap, extend='both', levels=np.linspace(-1, 1, 11))
        anom = ax.contour(comp.lon, comp.lat, comp.values, colors='black', levels=np.linspace(-1, 1, 11))
        ax.clabel(anom, inline=1, fontsize=8, fmt='%1.1f')
        plt.title(title + '  IOD ' + titulo2)
        plt.axhline(y=0, color='gray', alpha=0.5)
        ax.add_patch(
            patches.Rectangle(xy=(50, -10), width=20, height=20, linewidth=2, color='black', fill=False),
        )

        ax.add_patch(
            patches.Rectangle(xy=(90, -10), width=20, height=10, linewidth=2, color='black', fill=False),
        )

        plt.tight_layout()
        if save:
            plt.savefig(w_dir + 'DMI_Compos_' + season + nombre + phase + '.jpg', dpi=200)
            plt.close()
        else:
            plt.show()
            plt.close()

def xrFieldTimeDetrend_sst(xrda, dim, deg=1):
    # detrend along a single dimension
    aux = xrda.polyfit(dim=dim, deg=deg)
    trend = xr.polyval(xrda[dim], aux.sst_polyfit_coefficients)
    dt = xrda - trend
    return dt

def CorrENSO_IOD(enso, iod):
    seasons = [[12, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]
    r = []
    for s in range(5):
        aux_enso = enso.where(enso.time.dt.month.isin(seasons[s]), drop=True)
        aux_iod = iod.where(iod.time.dt.month.isin(seasons[s]), drop=True)
        r.append(pearsonr(aux_enso, aux_iod)[0])

    return r

def CorrENSO_IOD_lag(enso, iod, lag_menos_uno):
    seasons = [[3,4,5],[4,5,6],[5,6,7],[6, 7, 8], [7, 8, 9], [8, 9, 10], [9, 10, 11]]

    lag = lag_menos_uno
    r = []
    for s in range((7 - lag)):
        print(s)
        print(str(seasons[s]) + str(seasons[s+lag]))
        aux_enso = enso.where(enso.time.dt.month.isin(seasons[s]), drop=True)
        aux_iod = iod.where(iod.time.dt.month.isin(seasons[s + lag]), drop=True)
        r.append(np.around(pearsonr(aux_enso, aux_iod)[0],5))

    return r

def DMI(per = 0, filter_bwa = True, filter_harmonic = True,
        filter_all_harmonic=True, harmonics = [],
        start_per='1920', end_per='2020'):

    var = True
    if per == 2:
        movinganomaly = True
        start_year = '1906'
        end_year = '2020'
        change_baseline = False
        start_year2 = '1920'
        end_year2 = '2020_30r5'
        print('30r5')

    else:
        movinganomaly = False
        start_year = '1920'
        end_year = '2020'
        change_baseline = True
        start_year2 = start_per
        end_year2 = end_per
        print('All')

    ##################################### DATA #####################################
    # ERSSTv5
    sst = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")

    dataname = 'ERSST'
    ##################################### Pre-processing #####################################
    iodw = sst.sel(lat=slice(10.0, -10.0), lon=western_io,
                       time=slice(start_year + '-01-01', end_year + '-12-31'))
    iodw = iodw.sst.mean(['lon', 'lat'], skipna=True)
    iodw2 = iodw
    if per == 2:
        iodw2 = iodw2[168:]
    # -----------------------------------------------------------------------------------#
    iode = sst.sel(lat=slice(0, -10.0), lon=slice(90, 110),
                   time=slice(start_year + '-01-01', end_year + '-12-31'))
    iode = iode.sst.mean(['lon', 'lat'], skipna=True)
    # -----------------------------------------------------------------------------------#
    bwa = sst.sel(lat=slice(20.0, -20.0), lon=slice(40, 110),
                  time=slice(start_year + '-01-01', end_year + '-12-31'))
    bwa = bwa.sst.mean(['lon', 'lat'], skipna=True)
    # ----------------------------------------------------------------------------------#

    if movinganomaly:
        iodw = MovingBasePeriodAnomaly(iodw)
        iode = MovingBasePeriodAnomaly(iode)
        bwa = MovingBasePeriodAnomaly(bwa)
    else:
        # change baseline
        if change_baseline:
            iodw = iodw.groupby('time.month') - \
                   iodw.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                       'time')

            iode = iode.groupby('time.month') - \
                   iode.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                       'time')

            bwa = bwa.groupby('time.month') - \
                  bwa.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                      'time')
            print('baseline: ' + str(start_year2) + ' - ' + str(end_year2))
        else:
            print('baseline: All period')
            iodw = iodw.groupby('time.month') - iodw.groupby('time.month').mean('time', skipna=True)
            iode = iode.groupby('time.month') - iode.groupby('time.month').mean('time', skipna=True)
            bwa = bwa.groupby('time.month') - bwa.groupby('time.month').mean('time', skipna=True)

    # ----------------------------------------------------------------------------------#
    # Detrend
    iodw_trend = np.polyfit(range(0, len(iodw)), iodw, deg=1)
    iodw = iodw - (iodw_trend[0] * range(0, len(iodw)) + iodw_trend[1])
    # ----------------------------------------------------------------------------------#
    iode_trend = np.polyfit(range(0, len(iode)), iode, deg=1)
    iode = iode - (iode_trend[0] * range(0, len(iode)) + iode_trend[1])
    # ----------------------------------------------------------------------------------#
    bwa_trend = np.polyfit(range(0, len(bwa)), bwa, deg=1)
    bwa = bwa - (bwa_trend[0] * range(0, len(bwa)) + bwa_trend[1])
    # ----------------------------------------------------------------------------------#

    # 3-Month running mean
    iodw_filtered = np.convolve(iodw, np.ones((3,)) / 3, mode='same')
    iode_filtered = np.convolve(iode, np.ones((3,)) / 3, mode='same')
    bwa_filtered = np.convolve(bwa, np.ones((3,)) / 3, mode='same')

    # Common preprocessing, for DMIs other than SY2003a
    iode_3rm = iode_filtered
    iodw_3rm = iodw_filtered

    #################################### follow SY2003a #######################################

    # power spectrum
    # aux = FFT2(iodw_3rm, maxVar=20, maxA=15).sort_values('Variance', ascending=False)
    # aux2 = FFT2(iode_3rm, maxVar=20, maxA=15).sort_values('Variance', ascending=False)

    # filtering harmonic
    if filter_harmonic:
        if filter_all_harmonic:
            for harmonic in range(15):
                iodw_filtered = WaveFilter(iodw_filtered, harmonic)
                iode_filtered = WaveFilter(iode_filtered, harmonic)
            else:
                for harmonic in harmonics:
                    iodw_filtered = WaveFilter(iodw_filtered, harmonic)
                    iode_filtered = WaveFilter(iode_filtered, harmonic)

    ## max corr. lag +3 in IODW
    ## max corr. lag +6 in IODE

    # ----------------------------------------------------------------------------------#
    # ENSO influence
    # pre processing same as before
    if filter_bwa:
        ninio3 = sst.sel(lat=slice(5.0, -5.0), lon=slice(210, 270),
                         time=slice(start_year + '-01-01', end_year + '-12-31'))
        ninio3 = ninio3.sst.mean(['lon', 'lat'], skipna=True)

        if movinganomaly:
            ninio3 = MovingBasePeriodAnomaly(ninio3)
        else:
            if change_baseline:
                ninio3 = ninio3.groupby('time.month') - \
                         ninio3.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby(
                             'time.month').mean(
                             'time')

            else:

                ninio3 = ninio3.groupby('time.month') - ninio3.groupby('time.month').mean('time', skipna=True)

            trend = np.polyfit(range(0, len(ninio3)), ninio3, deg=1)
            ninio3 = ninio3 - (trend[0] * range(0, len(ninio3)) +trend[1])

        # 3-month running mean
        ninio3_filtered = np.convolve(ninio3, np.ones((3,)) / 3, mode='same')

        # ----------------------------------------------------------------------------------#
        # removing BWA effect
        # lag de maxima corr coincide para las dos bases de datos.
        lag = 3
        x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
        result = sm.ols(formula='bwa~ninio3', data=x).fit()
        recta = result.params[1] * ninio3_filtered + result.params[0]
        iodw_f = iodw_filtered - recta

        lag = 6
        x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
        result = sm.ols(formula='bwa~ninio3', data=x).fit()
        recta = result.params[1] * ninio3_filtered + result.params[0]
        iode_f = iode_filtered - recta
        print('BWA filtrado')
    else:
        iodw_f = iodw_filtered
        iode_f = iode_filtered
        print('BWA no filtrado')
    # ----------------------------------------------------------------------------------#

    # END processing
    if movinganomaly:
        iodw_3rm = xr.DataArray(iodw_3rm, coords=[iodw.time.values], dims=['time'])
        iode_3rm = xr.DataArray(iode_3rm, coords=[iodw.time.values], dims=['time'])

        iodw_f = xr.DataArray(iodw_f, coords=[iodw.time.values], dims=['time'])
        iode_f = xr.DataArray(iode_f, coords=[iodw.time.values], dims=['time'])
        start_year = '1920'
    else:
        iodw_3rm = xr.DataArray(iodw_3rm, coords=[iodw2.time.values], dims=['time'])
        iode_3rm = xr.DataArray(iode_3rm, coords=[iodw2.time.values], dims=['time'])

        iodw_f = xr.DataArray(iodw_f, coords=[iodw2.time.values], dims=['time'])
        iode_f = xr.DataArray(iode_f, coords=[iodw2.time.values], dims=['time'])

    ####################################### compute DMI #######################################

    dmi_sy_full = DMIndex(iodw_f, iode_f, method='sy2003', limitsize=len(iodw_f) - 2)

    return dmi_sy_full, sst, iodw_f-iode_f#, iodw_f - iode_f, iodw_f, iode_f

#  return dmi_sy_full, sst, dmi_serie, iodw_f, iode_f
##########################################################################################################################

# COMPOS 15S-15N 20-290
seasons_name = ['MJJ', 'JJA', 'ASO', 'SON']
bwa = [True, False]
bwa_title = ['fil BWA', 'NO fil BWA']
bwa_name = ['BWA', 'NO_BWA']
for filter in range(2):
    dmi_sy_full, sst, aux = DMI(filter_bwa=bwa[filter], per=0)

    for s in range(4):
        CompositeIOD(dmi=dmi_sy_full[0], original_data=sst, season=seasons_name[s],
                     title=bwa_title[filter] + ' 1920_2020 ' + seasons_name[s],
                     nombre='IndianPacific' + '_' + bwa_name[filter], save=True)

##########################################################################################################################


for var in ERSST:
    for per in baseline:
        if per == 2:
            movinganomaly=True
            start_year = '1906'
            end_year = '2020'
            change_baseline = False
            start_year2 = '1920'
            end_year2 = '2020_30r5'
            print('30r5')

        else:
            movinganomaly=False
            start_year = '1920'
            end_year = '2020'
            change_baseline = True
            start_year2 = start_per
            end_year2 = end_per
            print('All')

        ##################################### DATA #####################################
        if var:
            # ERSSTv5
            sst = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")

            dataname = 'ERSST'
            ##################################### Pre-processing #####################################
            iodw = sst.sel(lat=slice(10.0, -10.0), lon=western_io,
                           time=slice(start_year + '-01-01', end_year + '-12-31'))
            iodw = iodw.sst.mean(['lon', 'lat'], skipna=True)
            iodw2 = iodw
            if per ==2:
                iodw2 = iodw2[168:]
            #-----------------------------------------------------------------------------------#
            iode = sst.sel(lat=slice(0, -10.0), lon=slice(90, 110),
                           time=slice(start_year + '-01-01', end_year + '-12-31'))
            iode = iode.sst.mean(['lon', 'lat'], skipna=True)
            #-----------------------------------------------------------------------------------#
            bwa = sst.sel(lat=slice(20.0, -20.0), lon=slice(40, 110),
                          time=slice(start_year + '-01-01', end_year + '-12-31'))
            bwa = bwa.sst.mean(['lon', 'lat'], skipna=True)

        else:
            sst = xr.open_dataset("HadISST_sst.nc")
            sst = sst.rename({'latitude': 'lat'})
            sst = sst.rename({'longitude': 'lon'})
            dataname = 'HadISST'
            western_io = slice(50.5, 70.5)
            ##################################### Pre-processing #####################################
            iodw = sst.sel(lat=slice(10.5, -10.5), lon=western_io,
                           time=slice(start_year + '-01-01', end_year + '-12-31'))
            iodw = iodw.sst.mean(['lon', 'lat'], skipna=True)
            iodw2 = iodw
            if per ==2:
                iodw2 = iodw2[168:]
            # ----------------------------------------------------------------------------------#
            iode = sst.sel(lat=slice(-0.5, -10.5), lon=slice(89.5, 110.5),
                           time=slice(start_year + '-01-01', end_year + '-12-31'))
            iode = iode.sst.mean(['lon', 'lat'], skipna=True)
            # ----------------------------------------------------------------------------------#
            bwa = sst.sel(lat=slice(20.5, -20.5), lon=slice(40.5, 110.5),
                          time=slice(start_year + '-01-01', end_year + '-12-31'))
            bwa = bwa.sst.mean(['lon', 'lat'], skipna=True)

        # ----------------------------------------------------------------------------------#

        # Detrend
        iodw_trend = np.polyfit(range(0, len(iodw)), iodw, deg=1)
        iodw = iodw - iodw_trend[0] * range(0, len(iodw))
        # ----------------------------------------------------------------------------------#
        iode_trend = np.polyfit(range(0, len(iode)), iode, deg=1)
        iode = iode - iode_trend[0] * range(0, len(iode))
        # ----------------------------------------------------------------------------------#
        bwa_trend = np.polyfit(range(0, len(bwa)), bwa, deg=1)
        bwa = bwa - bwa_trend[0] * range(0, len(bwa))
        # ----------------------------------------------------------------------------------#
        if movinganomaly:
            iodw = MovingBasePeriodAnomaly(iodw)
            iode = MovingBasePeriodAnomaly(iode)
            bwa = MovingBasePeriodAnomaly(bwa)
        else:
            # change baseline
            if change_baseline:
                iodw = iodw.groupby('time.month') - \
                       iodw.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                           'time')

                iode = iode.groupby('time.month') - \
                       iode.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                           'time')

                bwa = bwa.groupby('time.month') - \
                      bwa.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby('time.month').mean(
                          'time')
                print('baseline: ' + str(start_year2) + ' - ' + str(end_year2))
            else:
                print('baseline: All period')
                iodw = iodw.groupby('time.month') - iodw.groupby('time.month').mean('time', skipna=True)
                iode = iode.groupby('time.month') - iode.groupby('time.month').mean('time', skipna=True)
                bwa = bwa.groupby('time.month') - bwa.groupby('time.month').mean('time', skipna=True)

        # ----------------------------------------------------------------------------------#
        # 3-Month running mean
        iodw_filtered = np.convolve(iodw, np.ones((3,)) / 3, mode='same')
        iode_filtered = np.convolve(iode, np.ones((3,)) / 3, mode='same')
        bwa_filtered = np.convolve(bwa, np.ones((3,)) / 3, mode='same')

        # Common preprocessing, for DMIs other than SY2003a
        iode_3rm = iode_filtered
        iodw_3rm = iodw_filtered

        #################################### follow SY2003a #######################################

        # power spectrum
        # aux = FFT2(iodw_3rm, maxVar=20, maxA=15).sort_values('Variance', ascending=False)
        # aux2 = FFT2(iode_3rm, maxVar=20, maxA=15).sort_values('Variance', ascending=False)

        # filtering harmonic
        if filter_harmonic:
            if filter_all_harmonic:
                for harmonic in range(15):
                    iodw_filtered = WaveFilter(iodw_filtered, harmonic)
                    iode_filtered = WaveFilter(iode_filtered, harmonic)
                else:
                    for harmonic in harmonics:
                        iodw_filtered = WaveFilter(iodw_filtered, harmonic)
                        iode_filtered = WaveFilter(iode_filtered, harmonic)

        ## max corr. lag +3 in IODW
        ## max corr. lag +6 in IODE

        # ----------------------------------------------------------------------------------#
        # ENSO influence
        # pre processing same as before
        if filter_bwa:
            if var:
                ninio3 = sst.sel(lat=slice(5.0, -5.0), lon=slice(210, 270),
                                 time=slice(start_year + '-01-01', end_year + '-12-31'))
            else:
                ninio3 = sst.sel(lat=slice(5.5, -5.5), lon=slice(-150.5, -90.5),
                                 time=slice(start_year + '-01-01', end_year + '-12-31'))

            ninio3 = ninio3.sst.mean(['lon', 'lat'], skipna=True)

            if movinganomaly:
                ninio3 = MovingBasePeriodAnomaly(ninio3)
            else:
                trend = np.polyfit(range(0, len(ninio3)), ninio3, deg=1)
                ninio3 = ninio3 - trend[0] * range(0, len(ninio3))
                if change_baseline:
                    ninio3 = ninio3.groupby('time.month') - \
                             ninio3.sel(time=slice(start_year2 + '-01-01', end_year2 + '-12-31')).groupby(
                                 'time.month').mean(
                                 'time')

                else:

                    ninio3 = ninio3.groupby('time.month') - ninio3.groupby('time.month').mean('time', skipna=True)

            # 3-month running mean
            ninio3_filtered = np.convolve(ninio3, np.ones((3,)) / 3, mode='same')
            # ----------------------------------------------------------------------------------#

            # removing BWA effect
            # lag de maxima corr coincide para las dos bases de datos.

            lag = 3
            x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
            result = sm.ols(formula='bwa~ninio3', data=x).fit()
            recta = result.params[1] * ninio3_filtered + result.params[0]
            iodw_f = iodw_filtered - recta

            lag = 6
            x = pd.DataFrame({'bwa': bwa_filtered[lag:], 'ninio3': ninio3_filtered[:-lag]})
            result = sm.ols(formula='bwa~ninio3', data=x).fit()
            recta = result.params[1] * ninio3_filtered + result.params[0]
            iode_f = iode_filtered - recta
            print('BWA filtrado')
        else:
            iodw_f = iodw_filtered
            iode_f = iode_filtered
            print('BWA no filtrado')
        # ----------------------------------------------------------------------------------#

        # END processing
        if movinganomaly:
            iodw_3rm = xr.DataArray(iodw_3rm, coords=[iodw.time.values], dims=['time'])
            iode_3rm = xr.DataArray(iode_3rm, coords=[iodw.time.values], dims=['time'])

            iodw_f = xr.DataArray(iodw_f, coords=[iodw.time.values], dims=['time'])
            iode_f = xr.DataArray(iode_f, coords=[iodw.time.values], dims=['time'])
            start_year = '1920'
        else:

            iodw_3rm = xr.DataArray(iodw_3rm, coords=[iodw2.time.values], dims=['time'])
            iode_3rm = xr.DataArray(iode_3rm, coords=[iodw2.time.values], dims=['time'])

            iodw_f = xr.DataArray(iodw_f, coords=[iodw2.time.values], dims=['time'])
            iode_f = xr.DataArray(iode_f, coords=[iodw2.time.values], dims=['time'])

        ####################################### compute DMI #######################################

        # Simple
        #dmi_sy_simple = DMIndex(iodw_3rm, iode_3rm, method='sy2003', limitsize=len(iodw_f) - 2)
        # dmi_sy_simple_notsdanom = DMIndex(iodw_3rm, iode_3rm, method='sy2003',
        #                                   sst_anom_sd=False, limitsize=len(iodw_f) - 2)

        # Full
        dmi_sy_full = DMIndex(iodw_f, iode_f, method='sy2003', limitsize=len(iodw_f) - 2)
        # dmi_sy_full_notsdanom = DMIndex(iodw_f, iode_f, method='sy2003',
        #                                 sst_anom_sd=False, limitsize=len(iodw_f) - 2)
        ####################################### ########## #######################################

        # COR
        dmi = iodw_f - iode_f
        aux = Nino34CPC(sst)[0]

        CorrENSO_IOD(enso=aux, iod=dmi)

        CorrENSO_IOD_lag(enso=aux, iod=dmi, lag_menos_uno=lag)


        ###################################### Plot #######################################
        seasons = ['MJJ', 'JJA', 'ASO', 'SON']
        for s in seasons:
            PlotDMI(dmi=SYpostprocessing(dmi_sy_full[0], full_season=False, season=s), end_year=end_year,
                    start_year=start_year,
                    end_year2=end_year2, start_year2=start_year2,
                    titulo='DMI - SY2003 full - ' + s + ' ' + dataname,
                    nombre='_SY_full_J-N' + end_year + '-' + end_year2 + dataname,
                    not_sy=False,
                    inch_y=inch_x, inch_x=inch_y)

##########################################################################################################################


# drafts
#
# sy_events = [1958,1960,1961,1963, 1964, 1967, 1971,1972,1974,1975,1977,1982,1983,
#              1989, 1992, 1993, 1994, 1996, 1997,1998]
#
# def SYearPlot(original_data, sy_years, mmin, mmax, dataname):
#
#     sst_original = original_data.sel(lat=slice(20, -30), lon=slice(30, 120),
#                                      time=slice('1958-01-01', '1999-12-31'))
#
#     sst_original = xrFieldTimeDetrend_sst(xrda=sst_original, dim='time')
#
#     x = np.arange(1958, 1999, 1)
#
#     mask = np.in1d(x, sy_years, invert=True)
#     neutro = sst_original.sel(time=sst_original.time.dt.year.isin(x[mask]))
#     neutro = neutro.sst.mean(['time'], skipna=True)
#
#     # esto no funca...
#     cmap = plt.cm.get_cmap("RdBu")
#     cmap.set_bad(color='gray')
#
#     def is_months(month):
#         return (month >= mmin) & (month <= mmax)
#
#     for l in range(len(sy_years)):
#         dmi_sign = sst_original.sel(time=sst_original.time.dt.year.isin(sy_years[l]),
#                                      lat=slice(20, -30), lon=slice(30, 120))
#
#         dmi_sign = dmi_sign.sel(time=is_months(dmi_sign['time.month']))
#         dmi_sign = dmi_sign.sst.mean(['time'], skipna=True)
#
#         comp = dmi_sign - neutro
#
#         fig, ax = plt.subplots(figsize=(6, 3), dpi=100)
#         im = comp.plot.contourf(cbar=cmap, extend='both', levels=np.linspace(-3, 3, 31))
#         anom = ax.contour(comp.lon, comp.lat, comp.values, colors='black', levels=np.linspace(-3, 3, 31))
#         ax.clabel(anom, inline=1, fontsize=8, fmt='%1.1f')
#         plt.title(str(sy_years[l]))
#         plt.axhline(y=0, color='gray', alpha=0.5)
#         ax.add_patch(
#             patches.Rectangle(xy=(50, -10), width=20, height=20, linewidth=2, color='black', fill=False),
#         )
#
#         ax.add_patch(
#             patches.Rectangle(xy=(90, -10), width=20, height=10, linewidth=2, color='black', fill=False),
#         )
#
#         plt.tight_layout()
#         #plt.savefig(w_dir + 'SY_' + dataname + '_' + str(sy_years[l]) + '.jpg', dpi=100)
#         plt.show()
#         plt.close()
#
# sy_events=[1998]
# SYearPlot(original_data=sst, sy_years=sy_events,mmin=6,mmax=11, dataname=dataname)
