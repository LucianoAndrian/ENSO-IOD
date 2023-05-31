"""
Anomalía de PP (y T cuando ande) en regiones de SA según la magnitud de los índices
"""
########################################################################################################################
import xarray as xr
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/' # índices por estaciones
fields_dir = '/pikachu/datos/luciano.andrian/cases_fields/' # campos de las variables PP  ( y T cuando ande)
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/bins/'
save = True
# Funciones ############################################################################################################

def SelectBins(serie, serie_sel, bmin, bmax):
    x = serie.where(serie_sel >= bmin)
    x = x.where(x < bmax)
    return x


def SelectDatesBins(bin, bin_data, r):
    aux = bin.sel(r=r)
    bin_data = bin_data.sel(r=r)
    dates_bins = aux.time[np.where(~np.isnan(aux.sst))]
    bin_data_f = bin_data.sel(time=bin_data.time.isin(dates_bins))
    return bin_data_f


def PlotBars(x, bin_n, bin_n_err, bin_n_len,
             bin_d, bin_d_err, bin_d_len,
             title='', name_fig='fig', out_dir=out_dir, save=False,
             ymin=-80, ymax=45):
    import matplotlib.pyplot as plt

    fig = plt.figure(1, figsize=(7, 7), dpi=300)
    ax = fig.add_subplot(111)
    plt.hlines(y=0, xmin=-4, xmax=4, color='k')
    ax.bar(x, bin_n, color='royalblue', alpha=1, width=0.15, label='Niño3.4')
    ax.errorbar(x, bin_n, yerr=bin_n_err, capsize=4, fmt='o', alpha=1,
                elinewidth=0.9, ecolor='navy', mfc='w', mec='navy', markersize=5)
    ax2 = ax.twinx()
    ax2.bar(x + 0.075, bin_n_len, color='royalblue', alpha=0.5, width=0.15)

    # ax.bar(x - 0.075, np.nan_to_num(bin_d), color='indianred', alpha=1, width=0.15, label='DMI')
    # ax.errorbar(x - 0.075, bin_d, yerr=bin_d_err, capsize=4, fmt='o', alpha=1,
    #             elinewidth=0.9, ecolor='firebrick', mec='firebrick', mfc='w', markersize=5)
    # ax2.bar(x - 0.075, bin_d_len, color='firebrick', alpha=0.5, width=0.15)

    ax.legend(loc='upper left')
    ax.set_ylim(ymin, ymax)
    ax2.set_ylim(0, 3000)

    ax.set_ylabel('prec anomaly [mm/month]', fontsize=10)
    ax2.set_ylabel('number of samples', fontsize=10)

    ax.set_xlabel('SST index (of std)', fontsize=15)
    ax.xaxis.set_tick_params(labelsize=12)
    ax.grid(True)

    plt.title(title, fontsize=15)
    plt.xlim(-3.5, 3.5)

    if save:
        plt.savefig(out_dir + name_fig)
    else:
        plt.show()

    plt.close('all')


def OpenDataSet(name, interp=False, lat_interp=None, lon_interp=None):

    if name == 'pp_gpcc':
        # GPCC2018
        aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_viejo/' + 'pp_gpcc.nc')
        pp_gpcc = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
        if interp:
            pp_gpcc = aux.interp(lon=lon_interp, lat=lat_interp)
        pp_gpcc = pp_gpcc.rename({'precip': 'var'})
        pp_gpcc = pp_gpcc.sel(time='1982-01-01')
        pp_gpcc = xr.where(np.isnan(pp_gpcc), np.nan, 1)

        return pp_gpcc

mask = OpenDataSet('pp_gpcc', interp=True,
                   lat_interp=np.linspace(-60,15,76),
                   lon_interp=np.linspace(275,330,56))

def Compute(box, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6, dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5,
            dmi_b_6,
            n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5,
            n34_b_6):
    r_count = 0
    for r in range(1, 25):
        if r_count == 0:
            d_b0 = SelectDatesBins(dmi_b0, box, r)

            d_b1 = SelectDatesBins(dmi_b1, box, r)
            d_b2 = SelectDatesBins(dmi_b2, box, r)
            d_b3 = SelectDatesBins(dmi_b3, box, r)
            d_b4 = SelectDatesBins(dmi_b4, box, r)
            d_b5 = SelectDatesBins(dmi_b5, box, r)
            d_b6 = SelectDatesBins(dmi_b6, box, r)

            d_b_1 = SelectDatesBins(dmi_b_1, box, r)
            d_b_2 = SelectDatesBins(dmi_b_2, box, r)
            d_b_3 = SelectDatesBins(dmi_b_3, box, r)
            d_b_4 = SelectDatesBins(dmi_b_4, box, r)
            d_b_5 = SelectDatesBins(dmi_b_5, box, r)
            d_b_6 = SelectDatesBins(dmi_b_6, box, r)

            n_b0 = SelectDatesBins(n34_b0, box, r)

            n_b1 = SelectDatesBins(n34_b1, box, r)
            n_b2 = SelectDatesBins(n34_b2, box, r)
            n_b3 = SelectDatesBins(n34_b3, box, r)
            n_b4 = SelectDatesBins(n34_b4, box, r)
            n_b5 = SelectDatesBins(n34_b5, box, r)
            n_b6 = SelectDatesBins(n34_b6, box, r)

            n_b_1 = SelectDatesBins(n34_b_1, box, r)
            n_b_2 = SelectDatesBins(n34_b_2, box, r)
            n_b_3 = SelectDatesBins(n34_b_3, box, r)
            n_b_4 = SelectDatesBins(n34_b_4, box, r)
            n_b_5 = SelectDatesBins(n34_b_5, box, r)
            n_b_6 = SelectDatesBins(n34_b_6, box, r)
            r_count = 1

        else:
            d_b0 = xr.concat([d_b0, SelectDatesBins(dmi_b0, box, r)], dim='time')

            d_b1 = xr.concat([d_b1, SelectDatesBins(dmi_b1, box, r)], dim='time')
            d_b2 = xr.concat([d_b2, SelectDatesBins(dmi_b2, box, r)], dim='time')
            d_b3 = xr.concat([d_b3, SelectDatesBins(dmi_b3, box, r)], dim='time')
            d_b4 = xr.concat([d_b4, SelectDatesBins(dmi_b4, box, r)], dim='time')
            d_b5 = xr.concat([d_b5, SelectDatesBins(dmi_b4, box, r)], dim='time')
            d_b6 = xr.concat([d_b6, SelectDatesBins(dmi_b4, box, r)], dim='time')

            d_b_1 = xr.concat([d_b_1, SelectDatesBins(dmi_b_1, box, r)], dim='time')
            d_b_2 = xr.concat([d_b_2, SelectDatesBins(dmi_b_2, box, r)], dim='time')
            d_b_3 = xr.concat([d_b_3, SelectDatesBins(dmi_b_3, box, r)], dim='time')
            d_b_4 = xr.concat([d_b_4, SelectDatesBins(dmi_b_4, box, r)], dim='time')
            d_b_5 = xr.concat([d_b_5, SelectDatesBins(dmi_b_3, box, r)], dim='time')
            d_b_6 = xr.concat([d_b_6, SelectDatesBins(dmi_b_4, box, r)], dim='time')

            d_b0_mean = d_b0.mean(skipna=True)

            d_b1_mean = d_b1.mean(skipna=True)
            d_b2_mean = d_b2.mean(skipna=True)
            d_b3_mean = d_b3.mean(skipna=True)
            d_b4_mean = d_b4.mean(skipna=True)
            d_b5_mean = d_b5.mean(skipna=True)
            d_b6_mean = d_b6.mean(skipna=True)

            d_b_1_mean = d_b_1.mean(skipna=True)
            d_b_2_mean = d_b_2.mean(skipna=True)
            d_b_3_mean = d_b_3.mean(skipna=True)
            d_b_4_mean = d_b_4.mean(skipna=True)
            d_b_5_mean = d_b_5.mean(skipna=True)
            d_b_6_mean = d_b_6.mean(skipna=True)

            d_b0_std = d_b0.std(skipna=True)

            d_b1_std = d_b1.std(skipna=True)
            d_b2_std = d_b2.std(skipna=True)
            d_b3_std = d_b3.std(skipna=True)
            d_b4_std = d_b4.std(skipna=True)
            d_b5_std = d_b5.std(skipna=True)
            d_b6_std = d_b6.std(skipna=True)

            d_b_1_std = d_b_1.std(skipna=True)
            d_b_2_std = d_b_2.std(skipna=True)
            d_b_3_std = d_b_3.std(skipna=True)
            d_b_4_std = d_b_4.std(skipna=True)
            d_b_5_std = d_b_5.std(skipna=True)
            d_b_6_std = d_b_6.std(skipna=True)

            n_b0 = xr.concat([n_b0, SelectDatesBins(n34_b0, box, r)], dim='time')

            n_b1 = xr.concat([n_b1, SelectDatesBins(n34_b1, box, r)], dim='time')
            n_b2 = xr.concat([n_b2, SelectDatesBins(n34_b2, box, r)], dim='time')
            n_b3 = xr.concat([n_b3, SelectDatesBins(n34_b3, box, r)], dim='time')
            n_b4 = xr.concat([n_b4, SelectDatesBins(n34_b4, box, r)], dim='time')
            n_b5 = xr.concat([n_b5, SelectDatesBins(n34_b4, box, r)], dim='time')
            n_b6 = xr.concat([n_b6, SelectDatesBins(n34_b4, box, r)], dim='time')

            n_b_1 = xr.concat([n_b_1, SelectDatesBins(n34_b_1, box, r)], dim='time')
            n_b_2 = xr.concat([n_b_2, SelectDatesBins(n34_b_2, box, r)], dim='time')
            n_b_3 = xr.concat([n_b_3, SelectDatesBins(n34_b_3, box, r)], dim='time')
            n_b_4 = xr.concat([n_b_4, SelectDatesBins(n34_b_4, box, r)], dim='time')
            n_b_5 = xr.concat([n_b_5, SelectDatesBins(n34_b_3, box, r)], dim='time')
            n_b_6 = xr.concat([n_b_6, SelectDatesBins(n34_b_4, box, r)], dim='time')

            n_b0_mean = n_b0.mean(skipna=True)

            n_b1_mean = n_b1.mean(skipna=True)
            n_b2_mean = n_b2.mean(skipna=True)
            n_b3_mean = n_b3.mean(skipna=True)
            n_b4_mean = n_b4.mean(skipna=True)
            n_b5_mean = n_b5.mean(skipna=True)
            n_b6_mean = n_b6.mean(skipna=True)

            n_b_1_mean = n_b_1.mean(skipna=True)
            n_b_2_mean = n_b_2.mean(skipna=True)
            n_b_3_mean = n_b_3.mean(skipna=True)
            n_b_4_mean = n_b_4.mean(skipna=True)
            n_b_5_mean = n_b_5.mean(skipna=True)
            n_b_6_mean = n_b_6.mean(skipna=True)

            n_b0_std = n_b0.std(skipna=True)

            n_b1_std = n_b1.std(skipna=True)
            n_b2_std = n_b2.std(skipna=True)
            n_b3_std = n_b3.std(skipna=True)
            n_b4_std = n_b4.std(skipna=True)
            n_b5_std = n_b5.std(skipna=True)
            n_b6_std = n_b6.std(skipna=True)

            n_b_1_std = n_b_1.std(skipna=True)
            n_b_2_std = n_b_2.std(skipna=True)
            n_b_3_std = n_b_3.std(skipna=True)
            n_b_4_std = n_b_4.std(skipna=True)
            n_b_5_std = n_b_5.std(skipna=True)
            n_b_6_std = n_b_6.std(skipna=True)

    bin_n = [n_b_6_mean.prec.values + 0, n_b_5_mean.prec.values + 0, n_b_4_mean.prec.values + 0,
             n_b_3_mean.prec.values + 0, n_b_2_mean.prec.values + 0, n_b_1_mean.prec.values + 0,
             n_b0_mean.prec.values + 0, n_b1_mean.prec.values + 0, n_b2_mean.prec.values + 0,
             n_b3_mean.prec.values + 0, n_b4_mean.prec.values + 0, n_b5_mean.prec.values + 0,
             n_b6_mean.prec.values + 0]

    bin_n_err = [n_b_6_std.prec.values + 0, n_b_5_std.prec.values + 0, n_b_4_std.prec.values + 0,
                 n_b_3_std.prec.values + 0, n_b_2_std.prec.values + 0, n_b_1_std.prec.values + 0,
                 n_b0_std.prec.values + 0, n_b1_std.prec.values + 0, n_b2_std.prec.values + 0,
                 n_b3_std.prec.values + 0, n_b4_std.prec.values + 0, n_b5_std.prec.values + 0,
                 n_b6_std.prec.values + 0]

    bin_d = [d_b_6_mean.prec.values + 0, d_b_5_mean.prec.values + 0, d_b_4_mean.prec.values + 0,
             d_b_3_mean.prec.values + 0, d_b_2_mean.prec.values + 0, d_b_1_mean.prec.values + 0,
             d_b0_mean.prec.values + 0, d_b1_mean.prec.values + 0, d_b2_mean.prec.values + 0,
             d_b3_mean.prec.values + 0, d_b4_mean.prec.values + 0, d_b5_mean.prec.values + 0,
             d_b6_mean.prec.values + 0]

    bin_d_err = [d_b_6_std.prec.values + 0, d_b_5_std.prec.values + 0, d_b_4_std.prec.values + 0,
                 d_b_3_std.prec.values + 0, d_b_2_std.prec.values + 0, d_b_1_std.prec.values + 0,
                 d_b0_std.prec.values + 0, d_b1_std.prec.values + 0, d_b2_std.prec.values + 0,
                 d_b3_std.prec.values + 0, d_b4_std.prec.values + 0, d_b5_std.prec.values + 0,
                 d_b6_std.prec.values + 0]

    bin_d_len = [len(d_b_6.time.values) + 0, len(d_b_5.time.values) + 0, len(d_b_4.time.values) + 0,
                 len(d_b_3.time.values) + 0, len(d_b_2.time.values) + 0, len(d_b_1.time.values) + 0,
                 len(d_b0.time.values) + 0, len(d_b1.time.values) + 0, len(d_b2.time.values) + 0,
                 len(d_b3.time.values) + 0, len(d_b4.time.values) + 0, len(d_b5.time.values) + 0,
                 len(d_b6.time.values) + 0]

    bin_n_len = [len(n_b_6.time.values) + 0, len(n_b_5.time.values) + 0, len(n_b_4.time.values) + 0,
                 len(n_b_3.time.values) + 0, len(n_b_2.time.values) + 0, len(n_b_1.time.values) + 0,
                 len(n_b0.time.values) + 0, len(n_b1.time.values) + 0, len(n_b2.time.values) + 0,
                 len(n_b3.time.values) + 0, len(n_b4.time.values) + 0, len(n_b5.time.values) + 0,
                 len(n_b6.time.values) + 0]

    return bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len
########################################################################################################################

seasons = ['JJA', 'JAS', 'ASO', 'SON']
for s in seasons:

    # indices
    data_dmi = xr.open_dataset(dates_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc')
    data_n34 = xr.open_dataset(dates_dir + 'N34_' + s + '_Leads_r_CFSv2.nc')
    # normalizando por sd
    data_dmi /= data_dmi.std(['time', 'r'])
    data_n34 /= data_n34.std(['time', 'r'])

    # prec
    data_prec = xr.open_dataset(fields_dir + 'prec_' + s.lower() + '.nc')
    data_prec *= 30
    # anom
    data_prec = data_prec - data_prec.mean(['time', 'r'])

    # mask
    mask_xr = np.ones((data_prec.dims['lat'], data_prec.dims['lon'])) * np.isfinite(data_prec.prec)
    mask_xr = xr.where(mask_xr == 0, np.nan, 1)

    for i in range(0, 156):  # q elegancia la de francia...
        for l in range(0, 24):
            mask_xr.values[i, l, :, :] = xr.where(~np.isnan(mask['var']), 1, np.nan)
    data_prec *= mask_xr

    # delimitando cada columna
    dmi_b0 = SelectBins(data_dmi, data_dmi, -0.25, 0.25)

    dmi_b1 = SelectBins(data_dmi, data_dmi, 0.25, 0.75)
    dmi_b2 = SelectBins(data_dmi, data_dmi, 0.75, 1.25)
    dmi_b3 = SelectBins(data_dmi, data_dmi, 1.25, 1.75)
    dmi_b4 = SelectBins(data_dmi, data_dmi, 1.75, 2.25)
    dmi_b5 = SelectBins(data_dmi, data_dmi, 2.25, 2.75)
    dmi_b6 = SelectBins(data_dmi, data_dmi, 2.75, 3.25)

    dmi_b_1 = SelectBins(data_dmi, data_dmi, -0.75, -0.25)
    dmi_b_2 = SelectBins(data_dmi, data_dmi, -1.25, -0.75)
    dmi_b_3 = SelectBins(data_dmi, data_dmi, -1.75, -1.25)
    dmi_b_4 = SelectBins(data_dmi, data_dmi, -2.25, -1.75)
    dmi_b_5 = SelectBins(data_dmi, data_dmi, -2.75, -2.25)
    dmi_b_6 = SelectBins(data_dmi, data_dmi, -2.25, -2.75)

    n34_b0 = SelectBins(data_n34, data_n34, -0.25, 0.25)

    n34_b1 = SelectBins(data_n34, data_n34, 0.25, 0.75)
    n34_b2 = SelectBins(data_n34, data_n34, 0.75, 1.25)
    n34_b3 = SelectBins(data_n34, data_n34, 1.25, 1.75)
    n34_b4 = SelectBins(data_n34, data_n34, 1.75, 2.25)
    n34_b5 = SelectBins(data_n34, data_n34, 2.25, 2.75)
    n34_b6 = SelectBins(data_n34, data_n34, 2.75, 3.25)

    n34_b_1 = SelectBins(data_n34, data_n34, -0.75, -0.25)
    n34_b_2 = SelectBins(data_n34, data_n34, -1.25, -0.75)
    n34_b_3 = SelectBins(data_n34, data_n34, -1.75, -1.25)
    n34_b_4 = SelectBins(data_n34, data_n34, -2.25, -1.75)
    n34_b_5 = SelectBins(data_n34, data_n34, -2.75, -2.25)
    n34_b_6 = SelectBins(data_n34, data_n34, -3.25, -2.75)

    x = np.linspace(-3, 3, 13)

    # Seleccion caja, computo y plot ----------------------------------------------------------------------------------#
    SESA = data_prec.sel(lat=slice(-39, -17), lon=slice(296, 315))
    SESA = SESA.mean(['lon', 'lat'])

    bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
        Compute(SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6,
                n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

    PlotBars(x, bin_n, bin_n_err, bin_n_len,
             bin_d, bin_d_err, bin_d_len,
             title='SESA - ' + s + ' - ' + 'Prec. Anomaly',
             name_fig='PREC_SESA_' + s + '_N34.jpg', save=save,
             ymax=60)

    S_SESA = data_prec.sel(lat=slice(-39, -29), lon=slice(296, 315))
    S_SESA = S_SESA.mean(['lon', 'lat'])

    bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
        Compute(S_SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6,
                n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

    PlotBars(x, bin_n, bin_n_err, bin_n_len,
             bin_d, bin_d_err, bin_d_len,
             title='S-SESA - ' + s + ' - ' + 'Prec. Anomaly',
             name_fig='PREC_S-SESA_' + s + '_N34.jpg', save=save,
             ymax=60)

    N_SESA = data_prec.sel(lat=slice(-29, -17), lon=slice(300, 320))
    N_SESA = N_SESA.mean(['lon', 'lat'])

    bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
        Compute(N_SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
                , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

    PlotBars(x, bin_n, bin_n_err, bin_n_len,
             bin_d, bin_d_err, bin_d_len,
             title='N-SESA - ' + s + ' - ' + 'Prec. Anomaly',
             name_fig='PREC_N-SESA_' + s + '_N34.jpg', save=save,
             ymax=60)



    #
    # C_BRAZIL = data_prec.sel(lat=slice(-20, -10), lon=slice(300, 325))
    # C_BRAZIL = C_BRAZIL.mean(['lon', 'lat'])
    #
    # bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
    #     Compute(C_BRAZIL, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
    #             , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
    #             n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
    #             , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)
    #
    # PlotBars(x, bin_n, bin_n_err, bin_n_len,
    #          bin_d, bin_d_err, bin_d_len,
    #          title='C-BRAZIL - ' + s + ' - ' + 'Prec. Anomaly',
    #          name_fig='PREC_C-BRAZIL_' + s + '.jpg', save=save,
    #          ymax=60)
    #
    #
    # ANDES_CHILE = data_prec.sel(lat=slice(-40, -30), lon=slice(285, 290))
    # ANDES_CHILE = ANDES_CHILE.mean(['lon', 'lat'])
    #
    # bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
    #     Compute(ANDES_CHILE, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
    #             , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
    #             n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
    #             , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)
    #
    # PlotBars(x, bin_n, bin_n_err, bin_n_len,
    #          bin_d, bin_d_err, bin_d_len,
    #          title='ANDES-CHILE - ' + s + ' - ' + 'Prec. Anomaly',
    #          name_fig='PREC_And-Chi_' + s + '.jpg', save=save,
    #          ymax=60)
    #
    # ANDES_CHILE_Sur = data_prec.sel(lat=slice(-60, -43), lon=slice(285, 288))
    # ANDES_CHILE = ANDES_CHILE_Sur.mean(['lon', 'lat'])
    #
    # bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
    #     Compute(ANDES_CHILE_Sur, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
    #             , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
    #             n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
    #             , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)
    #
    # PlotBars(x, bin_n, bin_n_err, bin_n_len,
    #          bin_d, bin_d_err, bin_d_len,
    #          title='ANDES-CHILE_Sur - ' + s + ' - ' + 'Prec. Anomaly',
    #          name_fig='PREC_And-ChiS_' + s + '.jpg', save=save,
    #          ymax=60)
    #
    # SESA = data_prec.sel(lat=slice(-35, -20), lon=slice(300, 320))
    # SESA = SESA.mean(['lon', 'lat'])
    #
    # bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
    #     Compute(SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
    #             , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
    #             n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
    #             , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)
    #
    # PlotBars(x, bin_n, bin_n_err, bin_n_len,
    #          bin_d, bin_d_err, bin_d_len,
    #          title='SESA* - ' + s + ' - ' + 'Prec. Anomaly',
    #          name_fig='PREC_SESA_' + s + '.jpg', save=save,
    #          ymax=60)
########################################################################################################################