########################################################################################################################
import xarray as xr
import numpy as np

########################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/'
dir_leads = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/bins/dc/'
save=True
########################################################################################################################

def SelectBins(serie, serie_sel, bmin, bmax):
    x = serie.where(serie_sel >= bmin)
    x = x.where(x < bmax)
    return x


def SelectDatesBins(bin, bin_data, r):
    aux = bin.sel(r=r)
    bin_data = bin_data.sel(r=r)
    dates_bins = aux.time[np.where(~np.isnan(aux.index))]
    bin_data_f = bin_data.sel(time=bin_data.time.isin(dates_bins))
    return bin_data_f


def PlotBars(x, bin_n, bin_n_err, bin_n_len,
             bin_d, bin_d_err, bin_d_len,
             title='', name_fig='fig', out_dir=out_dir, save=False,
             ymin=-80, ymax=45):
    import matplotlib.pyplot as plt

    fig = plt.figure(1, figsize=(7, 7), dpi=200)
    ax = fig.add_subplot(111)
    plt.hlines(y=0, xmin=-4, xmax=4, color='k')
    ax.bar(x + 0.075, bin_n, color='royalblue', alpha=1, width=0.15, label='Niño3.4')
    ax.errorbar(x + 0.075, bin_n, yerr=bin_n_err, capsize=4, fmt='o', alpha=1,
                elinewidth=0.9, ecolor='navy', mfc='w', mec='navy', markersize=5)
    ax2 = ax.twinx()
    ax2.bar(x + 0.075, bin_n_len, color='royalblue', alpha=0.5, width=0.15)

    ax.bar(x - 0.075, np.nan_to_num(bin_d), color='indianred', alpha=1, width=0.15, label='DMI')
    ax.errorbar(x - 0.075, bin_d, yerr=bin_d_err, capsize=4, fmt='o', alpha=1,
                elinewidth=0.9, ecolor='firebrick', mec='firebrick', mfc='w', markersize=5)
    ax2.bar(x - 0.075, bin_d_len, color='firebrick', alpha=0.5, width=0.15)

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
        aux = xr.open_dataset('/datos/luciano.andrian/ncfiles/' + 'pp_gpcc.nc')
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

    bin_n = [n_b_6_mean.variable.values + 0, n_b_5_mean.variable.values + 0, n_b_4_mean.variable.values + 0,
             n_b_3_mean.variable.values + 0, n_b_2_mean.variable.values + 0, n_b_1_mean.variable.values + 0,
             n_b0_mean.variable.values + 0, n_b1_mean.variable.values + 0, n_b2_mean.variable.values + 0,
             n_b3_mean.variable.values + 0, n_b4_mean.variable.values + 0, n_b5_mean.variable.values + 0,
             n_b6_mean.variable.values + 0]

    bin_n_err = [n_b_6_std.variable.values + 0, n_b_5_std.variable.values + 0, n_b_4_std.variable.values + 0,
                 n_b_3_std.variable.values + 0, n_b_2_std.variable.values + 0, n_b_1_std.variable.values + 0,
                 n_b0_std.variable.values + 0, n_b1_std.variable.values + 0, n_b2_std.variable.values + 0,
                 n_b3_std.variable.values + 0, n_b4_std.variable.values + 0, n_b5_std.variable.values + 0,
                 n_b6_std.variable.values + 0]

    bin_d = [d_b_6_mean.variable.values + 0, d_b_5_mean.variable.values + 0, d_b_4_mean.variable.values + 0,
             d_b_3_mean.variable.values + 0, d_b_2_mean.variable.values + 0, d_b_1_mean.variable.values + 0,
             d_b0_mean.variable.values + 0, d_b1_mean.variable.values + 0, d_b2_mean.variable.values + 0,
             d_b3_mean.variable.values + 0, d_b4_mean.variable.values + 0, d_b5_mean.variable.values + 0,
             d_b6_mean.variable.values + 0]

    bin_d_err = [d_b_6_std.variable.values + 0, d_b_5_std.variable.values + 0, d_b_4_std.variable.values + 0,
                 d_b_3_std.variable.values + 0, d_b_2_std.variable.values + 0, d_b_1_std.variable.values + 0,
                 d_b0_std.variable.values + 0, d_b1_std.variable.values + 0, d_b2_std.variable.values + 0,
                 d_b3_std.variable.values + 0, d_b4_std.variable.values + 0, d_b5_std.variable.values + 0,
                 d_b6_std.variable.values + 0]

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
mmonth_seasons_names = [0, 1, 2, 3]  # son las de abajo...
seasons = ['JJA', 'JAS', 'ASO', 'SON']
mmonth_seasons = [7, 8, 9, 10]
sets = [[0, 1, 2, 3]]

variables = ['prec', 'tref']

for m_name in mmonth_seasons_names:
    print(seasons[m_name])
    data_dmi_s = xr.open_dataset(dates_dir + seasons[m_name] + '_DMI_Leads_r_CFSv2.nc')
    data_n34_s = xr.open_dataset(dates_dir + seasons[m_name] + '_N34_Leads_r_CFSv2.nc')
    data_tref_s = xr.open_dataset(dir_leads + seasons[m_name] + '_tref_Leads_r_CFSv2.nc')
    data_prec_s = xr.open_dataset(dir_leads + seasons[m_name] + '_prec_Leads_r_CFSv2.nc')

    for s in sets:
        print('Set: ' + str(s))
        l = np.arange(len(s))
        ms = mmonth_seasons[m_name]
        data_dmi = data_dmi_s.sel(time=data_dmi_s.time.dt.month.isin(ms - l))
        data_n34 = data_n34_s.sel(time=data_n34_s.time.dt.month.isin(ms - l))

        dmi_sd = data_dmi.std(['time', 'r'])
        n34_sd = data_n34.std(['time', 'r'])
        # data_dmi = data_dmi.__mul__(1 / dmi_sd)
        # data_n34 = data_n34.__mul__(1 / n34_sd)

        data_prec = data_prec_s.sel(time=data_prec_s.time.dt.month.isin(ms - l))
        data_tref = data_tref_s.sel(time=data_tref_s.time.dt.month.isin(ms - l))

        data_prec *= 30
        data_prec = data_prec - data_prec.mean(['r', 'time'])

        data_tref = data_tref - data_tref.mean(['r', 'time'])

        data_prec = data_prec.rename({'prec': 'variable'})
        data_tref = data_tref.rename({'tref': 'variable'})

        mask_xr = np.ones((data_prec.dims['lat'], data_prec.dims['lon'])) * np.isfinite(data_prec.variable)
        mask_xr = xr.where(mask_xr == 0, np.nan, 1)
        data_prec *= mask_xr
        data_tref *= mask_xr

        # data_dmi = data_dmi_s.where(np.abs(data_dmi) > 0.75*data_dmi.std(['time','r']))
        #
        #
        # data_n34 = data_n34_s.where(np.abs(data_n34) > data_n34.std(['time','r']))

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

        S_SESA = data_prec.sel(lat=slice(-35, -29), lon=slice(300, 310))
        S_SESA = S_SESA.mean(['lon', 'lat'])

        S_SESA_exp = data_prec.sel(lat=slice(-35, -25), lon=slice(295, 315))
        S_SESA_exp = S_SESA_exp.mean(['lon', 'lat'])

        N_SESA = data_prec.sel(lat=slice(-29, -20), lon=slice(300, 320))
        N_SESA = N_SESA.mean(['lon', 'lat'])

        C_BRAZIL = data_prec.sel(lat=slice(-20, -10), lon=slice(300, 325))
        C_BRAZIL = C_BRAZIL.mean(['lon', 'lat'])

        ANDES_CHILE = data_prec.sel(lat=slice(-40, -30), lon=slice(285, 290))
        ANDES_CHILE = ANDES_CHILE.mean(['lon', 'lat'])

        # E_PATAGONIA = data_prec.sel(lat=slice(-50, -30), lon=slice(285, 290))
        # E_PATAGONIA = E_PATAGONIA.mean(['lon', 'lat'])

        # # S_SESA= data_prec.sel(lat=slice(-15,2), lon=slice(311, 325))
        # S_SESA = S_SESA.mean(['lon', 'lat'])

        S_SESA_t = data_tref.sel(lat=slice(-35, -29), lon=slice(300, 310))
        S_SESA_t = S_SESA_t.mean(['lon', 'lat'])

        S_SESA_exp_t = data_tref.sel(lat=slice(-35, -25), lon=slice(295, 315))
        S_SESA_exp_t = S_SESA_exp_t.mean(['lon', 'lat'])

        N_SESA_t = data_tref.sel(lat=slice(-29, -20), lon=slice(300, 320))
        N_SESA_t = N_SESA_t.mean(['lon', 'lat'])

        C_BRAZIL_t = data_tref.sel(lat=slice(-20, -10), lon=slice(300, 325))
        C_BRAZIL_t = C_BRAZIL_t.mean(['lon', 'lat'])

        ANDES_CHILE_t = data_tref.sel(lat=slice(-40, -30), lon=slice(285, 290))
        ANDES_CHILE_t = ANDES_CHILE_t.mean(['lon', 'lat'])

        # E_PATAGONIA_t = data_tref.sel(lat=slice(-50, -30), lon=slice(285, 290))
        # E_PATAGONIA_t = E_PATAGONIA_t.mean(['lon', 'lat'])


        x = np.linspace(-3, 3, 13)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(S_SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                    , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                    n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6,
                    n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='S-SESA - ' + seasons[m_name] + ' - ' + 'Prec. Anomaly',
                 name_fig='PREC_S-SESA_' + seasons[m_name] + '.jpg', save=save)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(S_SESA_exp, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                    , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                    n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6,
                    n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='S-SESA_exp - ' + seasons[m_name] + ' - ' + 'Prec. Anomaly',
                 name_fig='PREC_S-SESA_exp_' + seasons[m_name] + '.jpg', save=save)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(N_SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                    , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                    n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
                    , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='N-SESA - ' + seasons[m_name] + ' - ' + 'Prec. Anomaly',
                 name_fig='PREC_N-SESA_' + seasons[m_name] + '.jpg', save=save)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(ANDES_CHILE, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                    , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                    n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
                    , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='ANDES-CHILE - ' + seasons[m_name] + ' - ' + 'Prec. Anomaly',
                 name_fig='PREC_And-Chi_' + seasons[m_name] + '.jpg', save=save)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(C_BRAZIL, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                    , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                    n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
                    , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='C-BRAZIL - ' + seasons[m_name] + ' - ' + 'Prec. Anomaly',
                 name_fig='PREC_C-BRAZIL_' + seasons[m_name] + '.jpg', save=save)


        # bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
        #     Compute(E_PATAGONIA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
        #             , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
        #             n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
        #             , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)
        #
        # PlotBars(x, bin_n, bin_n_err, bin_n_len,
        #          bin_d, bin_d_err, bin_d_len,
        #          title='E-PATAGONIA - ' + seasons[m_name] + ' - ' + 'Prec. Anomaly',
        #          name_fig='PREC_E_PATAGONIA_' + seasons[m_name] + '.jpg', save=save)


        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(S_SESA_t, dmi_b0,dmi_b1,dmi_b2,dmi_b3,dmi_b4,dmi_b5,dmi_b6
                    ,dmi_b_1,dmi_b_2,dmi_b_3,dmi_b_4,dmi_b_5,dmi_b_6,
            n34_b0,n34_b1,n34_b2,n34_b3,n34_b4,n34_b5,n34_b6,n34_b_1
                    ,n34_b_2,n34_b_3,n34_b_4,n34_b_5,n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='S-SESA - ' + seasons[m_name] + ' - ' + 'T. Anomaly', ymin=-3, ymax=2,
                 name_fig='TREF_S-SESA_' + seasons[m_name] + '.jpg', save=save)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(S_SESA_exp, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
                    , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                    n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6,
                    n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='S-SESA_exp - ' + seasons[m_name] + ' - ' + 'T. Anomaly',  ymin=-3, ymax=2,
                 name_fig='Tref_S-SESA_exp_' + seasons[m_name] + '.jpg', save=save)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(N_SESA_t, dmi_b0,dmi_b1,dmi_b2,dmi_b3,dmi_b4,dmi_b5,dmi_b6
                    ,dmi_b_1,dmi_b_2,dmi_b_3,dmi_b_4,dmi_b_5,dmi_b_6,
            n34_b0,n34_b1,n34_b2,n34_b3,n34_b4,n34_b5,n34_b6,n34_b_1
                    ,n34_b_2,n34_b_3,n34_b_4,n34_b_5,n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='N-SESA - ' + seasons[m_name] + ' - ' + 'T. Anomaly', ymin=-3, ymax=2,
                 name_fig='TREF_N-SESA_' + seasons[m_name] + '.jpg', save=save)

        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(ANDES_CHILE_t, dmi_b0,dmi_b1,dmi_b2,dmi_b3,dmi_b4,dmi_b5,dmi_b6
                    ,dmi_b_1,dmi_b_2,dmi_b_3,dmi_b_4,dmi_b_5,dmi_b_6,
            n34_b0,n34_b1,n34_b2,n34_b3,n34_b4,n34_b5,n34_b6,n34_b_1
                    ,n34_b_2,n34_b_3,n34_b_4,n34_b_5,n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='ANDES-CHILE - ' + seasons[m_name] + ' - ' + 'T. Anomaly', ymin=-3, ymax=2,
                 name_fig='TREF_And-Chi_' + seasons[m_name] + '.jpg', save=save)


        bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
            Compute(C_BRAZIL_t, dmi_b0,dmi_b1,dmi_b2,dmi_b3,dmi_b4,dmi_b5,dmi_b6
                    ,dmi_b_1,dmi_b_2,dmi_b_3,dmi_b_4,dmi_b_5,dmi_b_6,
            n34_b0,n34_b1,n34_b2,n34_b3,n34_b4,n34_b5,n34_b6,n34_b_1
                    ,n34_b_2,n34_b_3,n34_b_4,n34_b_5,n34_b_6)

        PlotBars(x, bin_n, bin_n_err, bin_n_len,
                 bin_d, bin_d_err, bin_d_len,
                 title='C-BRAZIL- ' + seasons[m_name] + ' - ' + 'T. Anomaly', ymin=-3, ymax=2,
                 name_fig='TREF_C-BRAZIL_' + seasons[m_name] + '.jpg', save=save)

        # bin_n, bin_n_err, bin_n_len, bin_d, bin_d_err, bin_d_len = \
        #     Compute(E_PATAGONIA_t, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6
        #             , dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
        #             n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6, n34_b_1
        #             , n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)
        #
        # PlotBars(x, bin_n, bin_n_err, bin_n_len,
        #          bin_d, bin_d_err, bin_d_len,
        #          title='E-PATAGONIA - ' + seasons[m_name] + ' - ' + 'T. Anomaly',  ymin=-3, ymax=2,
        #          name_fig='Tref_E_PATAGONIA_' + seasons[m_name] + '.jpg', save=save)


#
# ########################################################################################################################
#         import matplotlib.pyplot as plt
#         x = np.linspace(-3,3,13)
#         y = np.zeros_like(x)
#         z = np.zeros_like(x)
#         dx = np.full_like(x, fill_value=0.5)
#         dy = dx
#         dz = bin_n
#
#         y2 = np.linspace(-3,3,13)
#         x2 = np.zeros_like(y2)
#         z2 = np.zeros_like(y2)
#         dy2 = np.full_like(y2, fill_value=0.5)
#         dx2 = dx
#         dz2 = bin_d
#
#         X, Y = np.meshgrid(np.linspace(-3,3,13), np.linspace(-3,3,13))
#         X, Z0 = np.meshgrid(np.linspace(-3,3,13), np.zeros_like(x))
#
#         Xv, Yv = np.meshgrid(bin_n, bin_d)
#         Z = (Xv+Yv)
#         import matplotlib.cm as cm
#
#         xpos = X.flatten() / 2.
#         ypos = Y.flatten() / 2.
#         zpos = np.zeros_like(xpos)
#
#         dx3 = X[1] - X[0]
#         dy3 = Y[1] - Y[0]
#         dz3 = Z.flatten()
#
#         cmap = cm.get_cmap('BrBG')  # Get desired colormap - you can change this!
#         max_height = np.max(dz3)  # get range of colorbars so we can normalize
#         min_height = np.min(dz3)
#         # scale each z to [0,1], and get their rgb values
#         rgba = [cmap((k - min_height) / (max_height-min_height)) for k in dz3]
#         cmap = cm.get_cmap('RdBu_r')  # Get desired colormap - you can change this!
#
#
#         fig = plt.figure(dpi=400)
#         ax = fig.add_subplot(projection='3d',  computed_zorder=False)
#         ax.bar3d(X.ravel(), Y.ravel(), Z0.ravel(), 0.5, 0.5, Z.ravel(), color=rgba, alpha=1, zorder=4)
#         #ax.bar3d(x, y, z, dx, dy, dz, zsort='average', color='blue', alpha=0.5, zorder=4.5, label='Niño3.4')
#         # ax.bar3d(x2, y2, z2, dx2, dy2, dz2, zsort='average', color='firebrick', alpha=0.5, zorder=4.4, label='DMI')
#
#         x_r = np.linspace(-3.5, 3.5, 13)
#         y_r = np.zeros_like(x_r)
#         ax.plot(x_r, y_r+3.5, 0, zdir='z', c='k', zorder=0)
#         ax.plot(y_r+3.5, x_r, 0, zdir='z', c='k', zorder=0)
#         ax.plot(x_r, y_r, -40.5, zdir='z', c='k', zorder=0)
#         ax.plot(y_r, x_r, -40.5, zdir='z', c='k', zorder=0)
#         ax.view_init(20, 30+180)
#         ax.set_zlim(-40,40)
#         ax.set_ylim(-3.5,3.5)
#         ax.set_xlim(-3.5,3.5)
#         ax.set_xlabel('Niño3.4 - SST index (of std)', fontsize=10)
#         ax.set_ylabel('DMI - SST index (of std)', fontsize=10)
#         #plt.legend(labels=['Niño3.4', 'DMI'])
#         plt.title(seasons[m_name])
#         plt.tight_layout()
#         plt.show()
# # ########################################################################################################################
#
#
#
#
