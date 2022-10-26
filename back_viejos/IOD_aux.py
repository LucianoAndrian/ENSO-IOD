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

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'


w_dir = '/home/luciano.andrian/doc/salidas/DMI/'
ERSST = [True, False]

western_io = slice(50,70)
start_per = ['1920','1981']
end_per = ['2020', '2010']

inch_y = 6
inch_x = 6

baseline = [0,1,2]

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

    plt.savefig(w_dir + 'DMI_' + nombre + '.jpg', dpi=200)
    #plt.show()
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
    plt.savefig(w_dir + 'DMI_TS_'+ label + '_' + nombre + '.jpg', dpi=200)
    plt.legend()
    #plt.show()
    plt.close()

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

        else:
            movinganomaly=False
            start_year = '1920'
            end_year = '2020'
            change_baseline = True
            start_year2 = start_per[per]
            end_year2 = end_per[per]

        ##################################### DATA #####################################
        if var:
            # ERSSTv5
            sst = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
            dataname = 'ERSST'
            western_io = slice(60, 80)
            ##################################### Pre-processing #####################################
            iodw = sst.sel(lat=slice(10.0, -10.0), lon=western_io,
                           time=slice(start_year + '-01-01', end_year + '-12-31'))
            iodw = iodw.sst.mean(['lon', 'lat'], skipna=True)
            iodw2 = iodw
            if per ==2:
                iodw2 = iodw2[168:]
            # ----------------------------------------------------------------------------------#
            iode = sst.sel(lat=slice(0, -10.0), lon=slice(90, 110),
                           time=slice(start_year + '-01-01', end_year + '-12-31'))
            iode = iode.sst.mean(['lon', 'lat'], skipna=True)
            # ----------------------------------------------------------------------------------#
            bwa = sst.sel(lat=slice(20.0, -20.0), lon=slice(40, 110),
                          time=slice(start_year + '-01-01', end_year + '-12-31'))
            bwa = bwa.sst.mean(['lon', 'lat'], skipna=True)

        else:
            sst = xr.open_dataset("HadISST_sst.nc")
            sst = sst.rename({'latitude': 'lat'})
            sst = sst.rename({'longitude': 'lon'})
            dataname = 'HadISST'
            western_io = slice(59.5, 80.5)
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

            # detrend

            # anomalies
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

        else:
            iodw_f = iodw_filtered
            iode_f = iode_filtered
            print('test')
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

        dmi_sy_simple = DMIndex(iodw_3rm, iode_3rm, method='sy2003', limitsize=len(iodw_f) - 2)
        dmi_sy_full = DMIndex(iodw_f, iode_f, method='sy2003', limitsize=len(iodw_f) - 2)

        dmi_sy_simple_notsdanom = DMIndex(iodw_3rm, iode_3rm, method='sy2003',
                                          sst_anom_sd=False, limitsize=len(iodw_f) - 2)
        dmi_sy_full_notsdanom = DMIndex(iodw_f, iode_f, method='sy2003',
                                        sst_anom_sd=False, limitsize=len(iodw_f) - 2)

        ####################################### Plot #######################################
        # seasons = [str('MJJ'), str('JJA'), str('ASO'), str('SON')]
        # for s in range(4):
        #     PlotDMI(dmi=SYpostprocessing(dmi_sy_full, full_season=False, season=seasons[s]), end_year=end_year,
        #             start_year=start_year,
        #             end_year2=end_year2, start_year2=start_year2,
        #             titulo='DMI - SY2003 full ' + seasons[s],
        #             nombre='SY_full_' + seasons[s] + end_year + '-' + end_year2 + dataname,
        #             not_sy=False,
        #             inch_y=inch_x, inch_x=inch_y)
        #
        #     PlotDMI(dmi=SYpostprocessing(dmi_sy_full_notsdanom, full_season=False, season=seasons[s]),
        #             end_year=end_year,
        #             start_year=start_year,
        #             end_year2=end_year2, start_year2=start_year2,
        #             titulo='DMI - SY2003 full ' + seasons[s] + ' no-SD',
        #             nombre='SY_full_' + seasons[s] + '_notsd' + end_year + '-' + end_year2 + dataname,
        #             not_sy=False,
        #             inch_y=inch_x, inch_x=inch_y)
        #
        #     PlotDMI(dmi=SYpostprocessing(dmi_sy_simple, full_season=False, season=seasons[s]), end_year=end_year,
        #             start_year=start_year,
        #             end_year2=end_year2, start_year2=start_year2,
        #             titulo='DMI - SY2003 simple ' + seasons[s],
        #             nombre='SY_simple_' + seasons[s] + end_year + '-' + end_year2 + dataname,
        #             not_sy=False,
        #             inch_y=inch_x, inch_x=inch_y)
        #
        #     PlotDMI(dmi=SYpostprocessing(dmi_sy_simple_notsdanom, full_season=False, season=seasons[s]),
        #             end_year=end_year,
        #             start_year=start_year,
        #             end_year2=end_year2, start_year2=start_year2,
        #             titulo='DMI - SY2003 simple ' + seasons[s] + ' no-SD',
        #             nombre='SY_simple_' + seasons[s] + '_notsd' + end_year + '-' + end_year2 + dataname,
        #             not_sy=False,
        #             inch_y=inch_x, inch_x=inch_y)
        #
        # PlotDMI(dmi=SYpostprocessing(dmi_sy_full, full_season=True), end_year=end_year,
        #         start_year=start_year,
        #         end_year2=end_year2, start_year2=start_year2,
        #         titulo='DMI - SY2003 full ',
        #         nombre='_SY_full_J-N' + end_year + '-' + end_year2 + dataname,
        #         not_sy=False,
        #         inch_y=inch_x, inch_x=inch_y)
        #
        # PlotDMI(dmi=SYpostprocessing(dmi_sy_full_notsdanom, full_season=True), end_year=end_year,
        #         start_year=start_year,
        #         end_year2=end_year2, start_year2=start_year2,
        #         titulo='DMI - SY2003 full no SD',
        #         nombre='_SY_full_J-N_notsd' + end_year + '-' + end_year2 + dataname,
        #         not_sy=False,
        #         inch_y=inch_x, inch_x=inch_y)
        #
        # PlotDMI(dmi=SYpostprocessing(dmi_sy_simple, full_season=True), end_year=end_year,
        #         start_year=start_year,
        #         end_year2=end_year2, start_year2=start_year2,
        #         titulo='DMI - SY2003 simple ',
        #         nombre='_SY_simple_J-N' + end_year + '-' + end_year2 + dataname,
        #         not_sy=False,
        #         inch_y=inch_x, inch_x=inch_y)
        #
        # PlotDMI(dmi=SYpostprocessing(dmi_sy_simple_notsdanom, full_season=True), end_year=end_year,
        #         start_year=start_year,
        #         end_year2=end_year2, start_year2=start_year2,
        #         titulo='DMI - SY2003 simple ' + ' no-SD',
        #         nombre='_SY_simple_J-N_notsd' + end_year + '-' + end_year2 + dataname,
        #         not_sy=False,
        #         inch_y=inch_x, inch_x=inch_y)
#
# #-----------------------------------------------------------------------------------------------------------------------
#         PlotSerieIOD(iodw=iodw_3rm, iode=iode_3rm, dmi=dmi_sy_simple_notsdanom,
#                      timeref=iodw2, start=120 * 0, end=120 * 2, nombre="SY_simple_noSD_1" + dataname + end_year2, var=var)
#
#         # PlotSerieIOD(iodw=iodw_3rm, iode=iode_3rm, dmi=dmi_sy_simple_notsdanom,
#         #              timeref=iodw2, start=120 * 4, end=120 * 2, nombre="SY_simple_noSD_2" + dataname + end_year2, var=var)
#
#         PlotSerieIOD(iodw=iodw_3rm, iode=iode_3rm, dmi=dmi_sy_simple_notsdanom,
#                      timeref=iodw2, start=120 * 8, end=1, nombre="SY_simple_noSD_3" + dataname + end_year2, var=var)
# #-----------------------------------------------------------------------------------------------------------------------
#         PlotSerieIOD(iodw=iodw_3rm, iode=iode_3rm, dmi=dmi_sy_simple,
#                      timeref=iodw2, start=120 * 0, end=120 * 2, nombre="SY_simple_1" + dataname + end_year2, var=var)
#
#         # PlotSerieIOD(iodw=iodw_3rm, iode=iode_3rm, dmi=dmi_sy_simple,
#         #              timeref=iodw2, start=120 * 4, end=120 * 2, nombre="SY_simple_2" + dataname + end_year2, var=var)
#
#         PlotSerieIOD(iodw=iodw_3rm, iode=iode_3rm, dmi=dmi_sy_simple,
#                      timeref=iodw2, start=120 * 8, end=1, nombre="SY_simple_3" + dataname + end_year2, var=var)
#-----------------------------------------------------------------------------------------------------------------------
#         PlotSerieIOD(iodw=iodw_f, iode=iode_f, dmi=dmi_sy_full_notsdanom,
#                      timeref=iodw2, start=120 * 0, end=120 * 2, nombre="SY_full_noSD_1" + dataname + end_year2, var=var)
#
#         # PlotSerieIOD(iodw=iodw_f, iode=iode_f, dmi=dmi_sy_simple,
#         #              timeref=iodw2, start=120 * 4, end=120 * 2, nombre="SY_full_noSD_2" + dataname + end_year2, var=var)
#
#         PlotSerieIOD(iodw=iodw_f, iode=iode_f, dmi=dmi_sy_full_notsdanom,
#                      timeref=iodw2, start=120 * 8, end=1, nombre="SY_full_noSD_3" + dataname + end_year2, var=var)
# # -----------------------------------------------------------------------------------------------------------------------
#         PlotSerieIOD(iodw=iodw_f, iode=iode_f, dmi=dmi_sy_full,
#                      timeref=iodw2, start=120 * 0, end=120 * 2, nombre="SY_full_1" + dataname + end_year2, var=var)
#         #
#         # PlotSerieIOD(iodw=iodw_f, iode=iode_f, dmi=dmi_sy_simple,
#         #              timeref=iodw2, start=120 * 4, end=120 * 2, nombre="SY_full_2" + dataname + end_year2, var=var)
#
#         PlotSerieIOD(iodw=iodw_f, iode=iode_f, dmi=dmi_sy_full,
#                      timeref=iodw2, start=120 * 8, end=1, nombre="SY_full_3" + dataname + end_year2, var=var)






def CompositeIOD(dmi, original_data=sst, full_season = False, season='JJA', title='title'):
    # only for sy2003
    # from June to November
    if full_season:
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['DMI'],
                        'Años': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['Años'],
                        'Mes': dmi.where(dmi.Mes.isin([6, 7, 8, 9, 10, 11]))['Mes']})
        aux = aux.groupby(['Años']).mean()
    elif season == 'MJJ':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([6]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([6]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([6]))['Mes']})

        aux = aux.groupby(['Años']).mean()

    elif season == 'JJA':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([7]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([7]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([7]))['Mes']})

        aux = aux.groupby(['Años']).mean()

    elif season == 'ASO':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([9]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([9]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([9]))['Mes']})

        aux = aux.groupby(['Años']).mean()


    elif season == 'SON':
        aux = pd.DataFrame({'DMI': dmi.where(dmi.Mes.isin([10]))['DMI'],
                            'Años': dmi.where(dmi.Mes.isin([10]))['Años'],
                            'Mes': dmi.where(dmi.Mes.isin([10]))['Mes']})

        aux = aux.groupby(['Años']).mean()


    dmi = aux

    dmi_positive = dmi.index.values[np.where(dmi['DMI'] > 0)]
    dmi_negative = dmi.index.values[np.where(dmi['DMI'] < 0)]

    sst_original = original_data.sel(lat=slice(20, -30), lon=slice(30, 120), time=slice('1920-01-01', '2020-12-31'))
    x = np.arange(1920, 2021, 1)

    mask = np.in1d(x, dmi.index.values, invert=True)
    neutro = sst_original.sel(time=sst_original.time.dt.year.isin(x[mask]))
    neutro = neutro.sst.mean(['time'], skipna=True)


    for sign in [0,1]:
        if sign==0:
            dmi_sign = original_data.sel(time=original_data.time.dt.year.isin([dmi_positive]),lat=slice(20, -30), lon=slice(30, 120))
            dmi_sign = dmi_sign.sst.mean(['time'], skipna=True)
            titulo2 = '+'
        else:
            dmi_sign = original_data.sel(time=original_data.time.dt.year.isin([dmi_negative]),lat=slice(20, -30), lon=slice(30, 120))
            dmi_sign = dmi_sign.sst.mean(['time'], skipna=True)
            titulo2 = '-'
#  ESTA RESTA ES PARA CADA EVENTO POR SEPARADO? O EL PROMEDIO DE POSITIVOS Y NEGATIVOS
        comp = dmi_sign - neutro

        fig, ax = plt.subplots(figsize=(6, 3), dpi=300)
        im = comp.plot.contourf(cbar=plt.get_cmap('RdBu'), extend='both', levels=np.linspace(-1.5, 1.5, 31))
        anom = ax.contour(comp.lon, comp.lat, comp.values, colors='black', levels=np.linspace(-1.5, 1.5, 31))
        ax.clabel(anom, inline=1, fontsize=8, fmt='%1.1f')
        plt.title(title + '  IOD ' + titulo2)
        plt.grid(True)
        ax.add_patch(
            patches.Rectangle( xy=(50, -10), width=20, height=20, linewidth=2, color='red', fill=False),
        )

        ax.add_patch(
            patches.Rectangle( xy=(90, -10), width=20, height=10, linewidth=2, color='red', fill=False),
        )


        plt.tight_layout()
        plt.show()

CompositeIOD(dmi=dmi_sy_full, original_data=sst, season='SON', title='SON')


