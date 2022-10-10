################################## FFT frequency ##################################
def FFT(serie, minA = 0, maxA = 20):

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

################################## FFT % variance ##################################
# more slow (can be very slow...)
def FFT2(serie, minA=0, maxA=20, maxVar=100, printstep=False):

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd

    N = np.size(serie)
    K = int(N/2)

    VAR = np.var(serie)

    CA = np.zeros(K)
    NARM = np.zeros(K)
    A = np.zeros(K)
    B = np.zeros(K)
    AM = np.zeros(K)
    C = np.zeros(K)
    CA = np.zeros(K)

    for i in range(K):
        NARM[i] = i
        SUM = 0
        SAM = 0
        if printstep:
            print(str(i))

        for j in range(N):
            SUM = SUM + serie[j] * np.sin(i * 2 * np.pi * (j) / N)
            SAM = SAM + (serie[j] * np.cos(i * 2 * np.pi * (j) / N))


        A[i] = 2*SUM/N
        B[i] = 2*SAM/N
        AM[i] = (A[i]**2 + B[i]**2)**0.5
        C[i] = (((AM[i]**2)/2)/VAR)*100
        CA[i] = np.sum(C)

    # last harmonic
    SUM = 0
    for j in range(N):
        SUM = SUM + serie[j]*np.cos(int(K)*2*np.pi*j/N)

    K = int(K) -1
    B[K] = SUM/N
    AM[K] = B[K]
    C[K] = ((AM[K]**2) / VAR)*100
    CA[K] = CA[K - 1] + C[K]
    A[K] = 0
    NARM[K]= K


    C = np.round(C,2)

    plt.bar(NARM[minA:maxA], C[minA:maxA])
    plt.ylim(0,maxVar)
    plt.xticks(np.arange(1,maxA,step=1))
    plt.ylabel('% of Variance')
    plt.xlabel('# harmonic')
    plt.grid()
    plt.show()

    results = pd.DataFrame({'Harmonic':NARM[minA:maxA], 'Variance': C[minA:maxA]})

    return(results)


################################## WaveFilter ##################################
# Harmonic Filter
def WaveFilter(serie, harmonic):

    import numpy as np

    sum = 0
    sam = 0
    N = np.size(serie)

    sum = 0
    sam = 0

    for j in range(N):
        sum = sum + serie[j] * np.sin(harmonic * 2 * np.pi * j / N)
        sam = sam + serie[j] * np.cos(harmonic * 2 * np.pi * j / N)

    A = 2*sum/N
    B = 2*sam/N

    xs = np.zeros(N)

    for j in range(N):
        xs[j] = A * np.sin(2 * np.pi * harmonic * j / N) + B * np.cos(2 * np.pi * harmonic * j / N)

    fil = serie - xs
    return(fil)


################################## DMIndex ##################################
# method: rao2002, sy2003, vf2005, v2018

def DMIndex(iodw, iode, sst_anom_sd=True, limitsize=1210, xsd=0.5, method='sy2003'):

    import numpy as np
    from itertools import groupby
    import pandas as pd

    # dipole mode index
    dmi = iodw - iode

    if method == 'sy2003':

        # criteria
        western_sign = np.sign(iodw)
        eastern_sign = np.sign(iode)
        opposite_signs = western_sign != eastern_sign

        sd = np.std(dmi) * xsd
        print(str(sd))
        sdw = np.std(iodw.values) * xsd
        sde = np.std(iode.values) * xsd

        results = []
        for k, g in groupby(enumerate(opposite_signs.values), key=lambda x: x[1]):
            if k:
                g = list(g)
                results.append([g[0][0], len(g)])


        iods = pd.DataFrame(columns=['DMI', 'Años', 'Mes'], dtype=float)
        dmi_raw = []
        for m in range(0, len(results)):
            # True values
            len_true = results[m][1]

            # True values for at least 3 consecutive seasons
            if len_true >= 3:

                for l in range(0, len_true):

                    if l < (len_true - 2):

                        main_month_num = results[m][0] + 1 + l
                        if main_month_num != limitsize:
                            main_month_name = dmi[main_month_num]['time.month'].values  # "name" 1 2 3 4 5

                            main_season = dmi[main_month_num]
                            b_season = dmi[main_month_num - 1]
                            a_season = dmi[main_month_num + 1]

                            # abs(dmi) > sd....(0.5*sd)
                            aux = (abs(main_season.values) > sd) & \
                                  (abs(b_season) > sd) & \
                                  (abs(a_season) > sd)

                            if sst_anom_sd:
                                if aux:
                                    sstw_main = iodw[main_month_num]
                                    sstw_b = iodw[main_month_num - 1]
                                    sstw_a = iodw[main_month_num + 1]
                                    #
                                    aux2 = (abs(sstw_main) > sdw) & \
                                           (abs(sstw_b) > sdw) & \
                                           (abs(sstw_a) > sdw)
                                    #
                                    sste_main = iode[main_month_num]
                                    sste_b = iode[main_month_num - 1]
                                    sste_a = iode[main_month_num + 1]

                                    aux3 = (abs(sste_main) > sde) & \
                                           (abs(sste_b) > sde) & \
                                           (abs(sste_a) > sde)

                                    if aux3 & aux2:
                                        iods = iods.append({'DMI': np.around(dmi[main_month_num].values, 2),
                                                            'Años': np.around(dmi[main_month_num]['time.year'].values),
                                                            'Mes': np.around(dmi[main_month_num]['time.month'].values)},
                                                           ignore_index=True)


                                        a = results[m][0]
                                        dmi_raw.append([np.arange(a, a + results[m][1]),
                                                        dmi[np.arange(a, a + results[m][1])].values])


                            else:
                                if aux:
                                    iods = iods.append({'DMI': np.around(dmi[main_month_num].values, 2),
                                                        'Años': np.around(dmi[main_month_num]['time.year'].values),
                                                        'Mes': np.around(dmi[main_month_num]['time.month'].values)},
                                                       ignore_index=True)



    elif method == 'rao2002':

        sd = np.std(dmi) * 1
        dmi_anual = dmi.groupby('time.year').mean()
        iods = pd.DataFrame(columns=['DMI', 'Años'], dtype=float)

        # hace esto mas rapido y facil?
        for i in range(len(dmi_anual)):
            if (np.abs(dmi_anual.values[i]) > sd):
                iods = iods.append({
                    'DMI':np.around(dmi_anual[i].values,2),
                    'Años': np.around(dmi_anual[i].year.values)},
                    ignore_index=True)


    elif method == 'vf2005':

        sd = np.std(dmi) * 0.5
        # average from june to August
        dmi = dmi.sel(time=dmi.time.dt.month.isin([6,7,8])).groupby('time.year').mean()

        iods = pd.DataFrame(columns=['DMI', 'Años'], dtype=float)
        for i in range(len(dmi)):
            if (np.abs(dmi.values[i]) > sd):
                iods = iods.append({
                    'DMI':np.around(dmi[i].values,2),
                    'Años': np.around(dmi[i]['year'].values)},
                    ignore_index=True)

    elif method == 'v2018':

        sd = np.std(dmi)*1
        # seems to averaged from june to november (unclear)
        dmi = dmi.sel(time=dmi.time.dt.month.isin([6,7,8,9,10,11])).groupby('time.year').mean()
        iods = pd.DataFrame(columns=['DMI', 'Años'], dtype=float)
        for i in range(len(dmi)):
            if (np.abs(dmi.values[i]) > sd):
                iods = iods.append({
                    'DMI':np.around(dmi[i].values,2),
                    'Años': np.around(dmi[i]['year'].values)},
                    ignore_index=True)


    return iods, dmi_raw


################################## MovingBasePeriodAnomaly ##################################

def MovingBasePeriodAnomaly(data, start='1920'):

    import xarray as xr
    # first five years
    start_num = int(start)

    initial = data.sel(time=slice(start + '-01-01', str(start_num + 5) + '-12-31')).groupby('time.month') - \
              data.sel(time=slice(str(start_num - 14) + '-01-01', str(start_num + 5 + 10) + '-12-31')).groupby(
                  'time.month').mean('time')


    start_num = start_num + 6
    result = initial

    while (start_num != 2016):

        aux = data.sel(time=slice(str(start_num) + '-01-01', str(start_num + 4) + '-12-31')).groupby('time.month') - \
              data.sel(time=slice(str(start_num - 15) + '-01-01', str(start_num + 4 + 10) + '-12-31')).groupby(
                  'time.month').mean('time')

        start_num = start_num + 5

        result = xr.concat([result, aux], dim='time')

    # 2016-2020 use base period 1991-2020
    aux = data.sel(time=slice(str(start_num) + '-01-01', str(start_num + 4) + '-12-31')).groupby('time.month') - \
          data.sel(time=slice('1991-01-01', '2020-12-31')).groupby('time.month').mean('time')

    result = xr.concat([result, aux], dim='time')

    return (result)

######################################### Nino34CPC #########################################

def Nino34CPC(data):

    import xarray as xr
    import numpy as np
    from itertools import groupby
    from Funciones import MovingBasePeriodAnomaly

    sst = data
    # N34
    ninio34 = sst.sel(lat=slice(4.0, -4.0), lon=slice(190, 240)
                      , time=slice('1906-01-01', '2020-12-31'))
    ninio34 = ninio34.sst.mean(['lon', 'lat'], skipna=True)

    # compute monthly anomalies
    ninio34 = MovingBasePeriodAnomaly(ninio34)

    # compute 5-month running mean
    ninio34_filtered = np.convolve(ninio34, np.ones((3,)) / 3, mode='same')  #
    ninio34_f = xr.DataArray(ninio34_filtered, coords=[ninio34.time.values], dims=['time'])

    aux = abs(ninio34_f) > 0.5
    results = []
    for k, g in groupby(enumerate(aux.values), key=lambda x: x[1]):
        if k:
            g = list(g)
            results.append([g[0][0], len(g)])

    n34 = []
    for m in range(0, len(results)):
        # True values
        len_true = results[m][1]

        # True values for at least 5 consecutive seasons
        if len_true >= 5:
            a = results[m][0]
            n34.append([np.arange(a, a + results[m][1]), ninio34_f[np.arange(a, a + results[m][1])].values])

    return ninio34_f, n34

