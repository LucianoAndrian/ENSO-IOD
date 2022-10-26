"""
Anomalía de PP (y T cuando ande) en regiones de SA según la magnitud de los índices
"""
########################################################################################################################
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/' # índices por estaciones
fields_dir = '/pikachu/datos/luciano.andrian/cases_fields/' # campos de las variables PP  ( y T cuando ande)
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/bins/2D/'
save = True
dpi=300
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

def Compute2(box, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6, dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5,
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

            d_b0_mean = d_b0

            d_b1_mean = d_b1
            d_b2_mean = d_b2
            d_b3_mean = d_b3
            d_b4_mean = d_b4
            d_b5_mean = d_b5
            d_b6_mean = d_b6

            d_b_1_mean = d_b_1
            d_b_2_mean = d_b_2
            d_b_3_mean = d_b_3
            d_b_4_mean = d_b_4
            d_b_5_mean = d_b_5
            d_b_6_mean = d_b_6

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

            n_b0_mean = n_b0

            n_b1_mean = n_b1
            n_b2_mean = n_b2
            n_b3_mean = n_b3
            n_b4_mean = n_b4
            n_b5_mean = n_b5
            n_b6_mean = n_b6

            n_b_1_mean = n_b_1
            n_b_2_mean = n_b_2
            n_b_3_mean = n_b_3
            n_b_4_mean = n_b_4
            n_b_5_mean = n_b_5
            n_b_6_mean = n_b_6

    bin_n = [n_b_6_mean.prec.values + 0, n_b_5_mean.prec.values + 0, n_b_4_mean.prec.values + 0,
             n_b_3_mean.prec.values + 0, n_b_2_mean.prec.values + 0, n_b_1_mean.prec.values + 0,
             n_b0_mean.prec.values + 0, n_b1_mean.prec.values + 0, n_b2_mean.prec.values + 0,
             n_b3_mean.prec.values + 0, n_b4_mean.prec.values + 0, n_b5_mean.prec.values + 0,
             n_b6_mean.prec.values + 0]


    bin_d = [d_b_6_mean.prec.values + 0, d_b_5_mean.prec.values + 0, d_b_4_mean.prec.values + 0,
             d_b_3_mean.prec.values + 0, d_b_2_mean.prec.values + 0, d_b_1_mean.prec.values + 0,
             d_b0_mean.prec.values + 0, d_b1_mean.prec.values + 0, d_b2_mean.prec.values + 0,
             d_b3_mean.prec.values + 0, d_b4_mean.prec.values + 0, d_b5_mean.prec.values + 0,
             d_b6_mean.prec.values + 0]


    aux_d_pos = []
    aux_d_neg = []
    aux_n_pos = []
    aux_n_neg = []

    for i in range(0, 13):
        aux_d_pos.append(bin_d[i][np.where(bin_d[i] > 0)])
        aux_d_neg.append(bin_d[i][np.where(bin_d[i] < 0)])

        aux_n_pos.append(bin_n[i][np.where(bin_n[i] > 0)])
        aux_n_neg.append(bin_n[i][np.where(bin_n[i] < 0)])

    return aux_d_pos, aux_d_neg, aux_n_pos, aux_n_neg, bin_d, bin_n


def PlotCuad(bin_d, bin_n, dpi=100, title='title', name_fig='fig', cmap='BrBG',
             levels=np.linspace(-1,1,11), save=save, num_max=800, cmap_num='Reds'):
    import matplotlib.pyplot as plt

    x = np.linspace(-3, 3, 13)
    X, Z0 = np.meshgrid(np.linspace(-3, 3, 13), np.zeros_like(x))
    X, Z_num = np.meshgrid(np.linspace(-3, 3, 13), np.zeros_like(x))

    Xv, Yv = np.meshgrid(bin_n, bin_d)
    Yv = Yv[::-1] # para que queden los valores negativos del DMI abajo en el grafico
    for i in range(0,13):
        for j in range(0,13):
            Z0[i,j] = np.nanmean(np.concatenate((Xv[i,j], Yv[i,j])))
            Z_num[i,j] = len(np.concatenate((Xv[i,j], Yv[i,j])))


    # Z = (Xv + Yv[::-1])
    Z = Z0

    vmin = min(levels)
    vmax = max(levels)

    fig = plt.figure(dpi=dpi, figsize=(7,7))
    x_r = np.linspace(0, 12, 13)
    y_r = [6,6,6,6,6,6,6,6,6,6,6,6,6]

    ax = fig.add_subplot(111)
    im=ax.imshow(Z, vmin=vmin, vmax=vmax,cmap=cmap)

    # ax.plot(x_r, y_r, c='k', alpha =.1)
    # ax.plot(y_r, x_r, c='k', alpha=.1)

    for i in range(0,13):
        for j in range(0,13):
            ax.text(j, i, Z_num[i,j].astype(np.int64), ha='center', va='center', color='k')

    ax.set_xticks(np.linspace(0,12,13), np.linspace(-3,3,13))
    ax.set_yticks(np.linspace(0, 12, 13), np.linspace(3, -3, 13))

    ax.set_xlabel('Niño3.4 - SST index (of std)', fontsize=11)
    ax.set_ylabel('DMI - SST index (of std)', fontsize=11)
    fig.suptitle(title, size=12)
    ax.set_xticks(np.arange(0, 12, 0.5), minor=True)
    ax.set_yticks(np.arange(0, 12, 0.5), minor=True)
    ax.grid(which='minor', alpha=0.5, color='k')
    plt.colorbar(im,ticks=levels,# boundaries=levels,
                 fraction=0.046, pad=0.04)

    plt.tight_layout()
    if save:
        plt.savefig(out_dir + name_fig)
        plt.close('all')
    else:
        plt.show()

def PlotCuadNumSamples(bin_d, bin_n, dpi=100, title='title', name_fig='fig', cmap='BrBG',
             levels=np.linspace(-1,1,11), save=save):
    import matplotlib.pyplot as plt

    x = np.linspace(-3, 3, 13)
    X, Z0 = np.meshgrid(np.linspace(-3, 3, 13), np.zeros_like(x))
    X, Z_num = np.meshgrid(np.linspace(-3, 3, 13), np.zeros_like(x))

    Xv, Yv = np.meshgrid(bin_n, bin_d)
    Yv = Yv[::-1] # para que queden los valores negativos del DMI abajo en el grafico
    for i in range(0,13):
        for j in range(0,13):
            Z0[i,j] = np.nanmean(np.concatenate((Xv[i,j], Yv[i,j])))
            Z_num[i,j] = len(np.concatenate((Xv[i,j], Yv[i,j])))


    # Z = (Xv + Yv[::-1])
    Z = Z0

    vmin = min(levels)
    vmax = max(levels)

    fig = plt.figure(dpi=dpi, figsize=(7,7))
    x_r = np.linspace(0, 12, 13)
    y_r = [6,6,6,6,6,6,6,6,6,6,6,6,6]

    ax = fig.add_subplot(111)
    im=ax.imshow(Z_num, vmin=vmin, vmax=vmax,cmap=cmap)

    # ax.plot(x_r, y_r, c='k', alpha =.1)
    # ax.plot(y_r, x_r, c='k', alpha=.1)

    for i in range(0,13):
        for j in range(0,13):
            ax.text(j, i, Z_num[i,j].astype(np.int64), ha='center', va='center', color='k')

    ax.set_xticks(np.linspace(0,12,13), np.linspace(-3,3,13))
    ax.set_yticks(np.linspace(0, 12, 13), np.linspace(3, -3, 13))

    ax.set_xlabel('Niño3.4 - SST index (of std)', fontsize=11)
    ax.set_ylabel('DMI - SST index (of std)', fontsize=11)
    fig.suptitle('Number of Samples - ' + title, size=12)
    ax.set_xticks(np.arange(0, 12, 0.5), minor=True)
    ax.set_yticks(np.arange(0, 12, 0.5), minor=True)
    ax.grid(which='minor', alpha=0.5, color='k')
    plt.colorbar(im,ticks=levels,# boundaries=levels,
                 fraction=0.046, pad=0.04)

    plt.tight_layout()
    if save:
        plt.savefig(out_dir + name_fig)
        plt.close('all')
    else:
        plt.show()

def PlotSeries(bin_d, bin_n, dpi=100, ymin=0, ymax=50,
                   title='title', name_fig='title', save=save):

        bin_d_mean = [np.nanmean(i) for i in bin_d]
        bin_d_std = [np.nanstd(i) for i in bin_d]

        bin_n_mean = [np.nanmean(i) for i in bin_n]
        bin_n_std = [np.nanstd(i) for i in bin_n]

        fig = plt.figure(dpi=dpi, figsize=(7, 5))
        ax = fig.add_subplot(111)
        ax.plot(bin_n_mean, label='N34', color='royalblue', linewidth=2.5)
        ax.fill_between(np.linspace(0, 12, 13),
                        [a + 0.5 * b for a, b in zip(bin_n_mean, bin_n_std)],
                        [a - 0.5 * b for a, b in zip(bin_n_mean, bin_n_std)],
                        alpha=0.2, color='royalblue')

        ax.plot(bin_d_mean, label='DMI', color='indianred', linewidth=2.5)
        ax.fill_between(np.linspace(0, 12, 13),
                        [a + 0.5 * b for a, b in zip(bin_d_mean, bin_d_std)],
                        [a - 0.5 * b for a, b in zip(bin_d_mean, bin_d_std)],
                        alpha=0.2, color='indianred')

        ax2 = ax.twinx()
        ax2.bar(np.linspace(0, 12, 13) - 0.10, [len(i) for i in bin_d], color='indianred', alpha=0.5, width=0.2)
        ax2.bar(np.linspace(0, 12, 13) + 0.10, [len(i) for i in bin_n], color='royalblue', alpha=0.5, width=0.2)

        ax.set_ylim(ymin, ymax)
        ax2.set_ylim(0, 1200)
        ax2.set_yticks(np.linspace(0,800,9))

        ax2.set_ylabel('number of samples', fontsize=10)

        ax.legend(loc='upper left')
        #plt.grid()
        fig.suptitle(title)

        ax.set_xticks(np.linspace(0, 12, 13), np.linspace(-3, 3, 13))
        #ax.set_xticks(np.arange(0, 12, 0.5), minor=True)
        ax.grid(which='minor', alpha=0.5, color='k')
        ax.grid(which='major', alpha=0.5, color='k')
        ax2.grid(which='major', alpha=0.6, color='gray', linestyle='--')
        ax.set_xlabel('std')
        ax.set_ylabel('[mm]')

        if save:
            plt.savefig(out_dir + name_fig)
            plt.close('all')
        else:
            plt.show()

from matplotlib import colors
color = ['white', '#F9FD7C', plt.cm.tab20c(7), plt.cm.tab20c(6), plt.cm.tab20c(5), plt.cm.tab20c(4),
             plt.cm.tab20b(15), plt.cm.tab20b(14), plt.cm.tab20b(13), plt.cm.tab20b(12)]
cmap_pos = colors.ListedColormap(color)
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
    S_SESA = data_prec.sel(lat=slice(-35, -29), lon=slice(300, 310))
    S_SESA = S_SESA.mean(['lon', 'lat'])

    bin_d_pos, bin_d_neg, bin_n_pos, bin_n_neg, bin_d, bin_n = \
        Compute2(S_SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6,
                 dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                 n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6,
                 n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

    PlotCuad(bin_d_pos, bin_n_pos, dpi=dpi,
             title='S-SESA - ' + s + ' -  PP>0 [mm]',
             name_fig='Cuad_S-SESA_' + s + 'pos.jpg',
             cmap='PuBuGn', levels=np.linspace(20, 50, 4), save=save)
    PlotCuadNumSamples(bin_d_pos, bin_n_pos, dpi=dpi,
             title='S-SESA - ' + s + ' -  PP>0 [mm]',
             name_fig='num_Cuad_S-SESA_' + s + 'pos.jpg',
             cmap=cmap_pos, levels=np.linspace(0, 800, 11), save=save)
    PlotSeries(bin_d_pos, bin_n_pos,dpi=dpi,
               title='S-SESA - ' + s + ' -  PP>0 [mm]',
               name_fig='Serie_S-SESA_' + s + 'pos.jpg',
               ymin=-10, ymax=70, save=save)


    PlotCuad(bin_d_neg, bin_n_neg,dpi=dpi,
             title='S-SESA - ' + s + ' -  PP<0 [mm]',
             name_fig='Cuad_S-SESA_' + s + 'neg.jpg',
             cmap='afmhot', levels=np.linspace(-40, -20, 6), save=save)
    PlotCuadNumSamples(bin_d_neg, bin_n_neg,dpi=dpi,
             title='S-SESA - ' + s + ' -  PP<0 [mm]',
             name_fig='num_Cuad_S-SESA_' + s + 'neg.jpg',
             cmap=cmap_pos, levels=np.linspace(0, 1600, 11), save=save)
    PlotSeries(bin_d_neg, bin_n_neg,dpi=dpi,
               title='S-SESA - ' + s + ' -  PP<0 [mm]',
               name_fig='Serie_S-SESA_' + s + 'neg.jpg',
               ymin=-70, ymax=0, save=save)


    PlotCuad(bin_d_pos, bin_n_neg,dpi=dpi,
             title='S-SESA - ' + s + ' -  N34 PP<0 - DMI PP>0 [mm]',
             name_fig='Cuad_S-SESA_' + s + 'N34neg_DMIpos.jpg',
             cmap='BrBG', levels=np.linspace(-40, 40, 5), save=save)
    PlotCuadNumSamples(bin_d_pos, bin_n_neg,dpi=dpi,
             title='S-SESA - ' + s + ' -  N34 PP<0 - DMI PP>0 [mm]',
             name_fig='num_Cuad_S-SESA_' + s + 'N34neg_DMIpos.jpg',
             cmap=cmap_pos, levels=np.linspace(0, 800, 11), save=save)
    PlotSeries(bin_d_pos, bin_n_neg, dpi=dpi,
               title='S-SESA - ' + s + ' -  N34 PP<0 - DMI PP>0 [mm]',
               name_fig='Serie_S-SESA_' + s + 'N34neg_DMIpos.jpg',
               ymin=-80, ymax=60, save=save)


    PlotCuad(bin_d_neg, bin_n_pos,dpi=dpi,
             title='S-SESA - ' + s + ' -  N34 PP>0 - DMI PP<0 [mm]',
             name_fig='Cuad_S-SESA_' + s + 'N34pos_DMIneg.jpg',
             cmap='BrBG', levels=np.linspace(-40, 40, 5), save=save)
    PlotCuadNumSamples(bin_d_neg, bin_n_pos,dpi=dpi,
             title='S-SESA - ' + s + ' -  N34 PP>0 - DMI PP<0 [mm]',
             name_fig='num_Cuad_S-SESA_' + s + 'N34pos_DMIneg.jpg',
             cmap=cmap_pos, levels=np.linspace(0, 800, 11), save=save)
    PlotSeries(bin_d_neg, bin_n_pos,dpi=dpi,
               title='S-SESA - ' + s + ' -  N34 PP>0 - DMI PP<0 [mm]',
               name_fig='Serie_S-SESA_' + s + 'N34pos_DMIneg.jpg',
               ymin=-80, ymax=60, save=save)


    PlotCuad(bin_d, bin_n, dpi=dpi,
             title='S-SESA - ' + s + ' -  N34 DMI [mm]',
             name_fig='Cuad_S-SESA_' + s + 'DMIN34.jpg',
             cmap='BrBG', levels=np.linspace(-20, 20, 5), save=save)
    PlotCuadNumSamples(bin_d, bin_n, dpi=dpi,
             title='S-SESA - ' + s + ' -  N34 DMI [mm]',
             name_fig='num_Cuad_S-SESA_' + s + 'DMIN34.jpg',
             cmap=cmap_pos, levels=np.linspace(0, 1600, 11), save=save)
    PlotSeries(bin_d, bin_n, dpi=dpi,
               title='S-SESA - ' + s + ' -  N34 DMI [mm]',
               name_fig='Serie_S-SESA_' + s + 'DMIN34.jpg',
               ymin=-80, ymax=60, save=save)


    N_SESA = data_prec.sel(lat=slice(-29, -20), lon=slice(300, 320))
    N_SESA = N_SESA.mean(['lon', 'lat'])

    bin_d_pos, bin_d_neg, bin_n_pos, bin_n_neg, bin_d, bin_n = \
        Compute2(N_SESA, dmi_b0, dmi_b1, dmi_b2, dmi_b3, dmi_b4, dmi_b5, dmi_b6,
                 dmi_b_1, dmi_b_2, dmi_b_3, dmi_b_4, dmi_b_5, dmi_b_6,
                 n34_b0, n34_b1, n34_b2, n34_b3, n34_b4, n34_b5, n34_b6,
                 n34_b_1, n34_b_2, n34_b_3, n34_b_4, n34_b_5, n34_b_6)

    PlotCuad(bin_d_pos, bin_n_pos, dpi=dpi,
             title='N-SESA - ' + s + ' -  PP>0 [mm]',
             name_fig='Cuad_N-SESA_' + s + 'pos.jpg',
             cmap='PuBuGn', levels=np.linspace(20, 50, 4), save=save)
    PlotCuadNumSamples(bin_d_pos, bin_n_pos, dpi=dpi,
             title='N-SESA - ' + s + ' -  PP>0 [mm]',
             name_fig='num_Cuad_N-SESA_' + s + 'pos.jpg',
             cmap=cmap_pos, levels=np.linspace(0, 800, 11), save=save)
    PlotSeries(bin_d_pos, bin_n_pos,dpi=dpi,
               title='N-SESA - ' + s + ' -  PP>0 [mm]',
               name_fig='Serie_N-SESA_' + s + 'pos.jpg',
               ymin=0, ymax=70, save=save)


    PlotCuad(bin_d_neg, bin_n_neg,dpi=dpi,
             title='N-SESA - ' + s + ' -  PP<0 [mm]',
             name_fig='Cuad_N-SESA_' + s + 'neg.jpg',
             cmap='afmhot', levels=np.linspace(-40, -20, 6), save=save)
    PlotSeries(bin_d_neg, bin_n_neg,dpi=dpi,
               title='N-SESA - ' + s + ' -  PP<0 [mm]',
               name_fig='Serie_N-SESA_' + s + 'neg.jpg',
               ymin=-70, ymax=0, save=save)


    PlotCuad(bin_d_pos, bin_n_neg,dpi=dpi,
             title='N-SESA - ' + s + ' -  N34 PP<0 - DMI PP>0 [mm]',
             name_fig='Cuad_N-SESA_' + s + 'N34neg_DMIpos.jpg',
             cmap='BrBG', levels=np.linspace(-40, 40, 5), save=save)
    PlotSeries(bin_d_pos, bin_n_neg, dpi=dpi,
               title='N-SESA - ' + s + ' -  N34 PP<0 - DMI PP>0 [mm]',
               name_fig='Serie_N-SESA_' + s + 'N34neg_DMIpos.jpg',
               ymin=-60, ymax=60, save=save)


    PlotCuad(bin_d_neg, bin_n_pos,dpi=dpi,
             title='N-SESA - ' + s + ' -  N34 PP>0 - DMI PP<0 [mm]',
             name_fig='Cuad_N-SESA_' + s + 'N34pos_DMIneg.jpg',
             cmap='BrBG', levels=np.linspace(-40, 40, 5), save=save)
    PlotSeries(bin_d_neg, bin_n_pos,dpi=dpi,
               title='N-SESA - ' + s + ' -  N34 PP>0 - DMI PP<0 [mm]',
               name_fig='Serie_N-SESA_' + s + 'N34pos_DMIneg.jpg',
               ymin=-60, ymax=60, save=save)


    PlotCuad(bin_d, bin_n, dpi=dpi,
             title='N-SESA - ' + s + ' -  N34 DMI [mm]',
             name_fig='Cuad_N-SESA_' + s + 'DMIN34.jpg',
             cmap='BrBG', levels=np.linspace(-20, 20, 5), save=save)
    PlotSeries(bin_d, bin_n, dpi=dpi,
               title='N-SESA - ' + s + ' -  N34 DMI [mm]',
               name_fig='Serie_N-SESA_' + s + 'DMIN34.jpg',
               ymin=-60, ymax=60, save=save)
########################################################################################################################