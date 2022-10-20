"""
PDF de las anomalías de PP (y T cuando ande...) en distans regiones de SA para cada caso
de IOD, ENSO y ENSO-IOD
"""
########################################################################################################################
import xarray as xr
from matplotlib import colors
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
warnings.filterwarnings('ignore')
import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
from threadpoolctl import threadpool_limits
from sklearn.cluster import KMeans
########################################################################################################################
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/PDFs/'

save = True
dpi = 200
# Funciones ############################################################################################################
def best_fit_distribution(data, size, start, end):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    #y, x = np.histogram(data, bins=bins, density=True, weights=np.ones(len(data))/len(data)*100)
    y, x = np.histogram(data)
    x = (x + np.roll(x, -1))[:-1] / 2.0

    # Best holders
    best_distributions = []

    # Estimate distribution parameters from data
    # for ii, distribution_name in enumerate([d for d in _distn_names if not d in ['arcsine', 'beta', 'laplace', 'rdist', 'reciprocal' ,
    #                                                                         'halfcauchy' ,'bradford', 'genpareto' ,'gennorm',
    #                                                                         'genexpon', 'dgamma', 'foldcauchy', 'expon' ,'dweibull',
    #                                                                         'cauchy' ,'levy_stable', 'studentized_range', 'wrapcauchy',
    #                                                                         'wald', 'uniform', 'tukeylambda', 'truncnorm', 'truncexpon',
    #                                                                         'triang', 'semicircular', 'reciprocal', 'rdist', 'powerlaw',
    #                                                                         'pareto', 'lomax', 'lognorm', 'loglaplace', 'levy_l','levy',
    #                                                                              'gausshyper', 'foldnorm']]):

                                                                            # algunas distribuciones que no se tienen en cuenta
    for distribution_name in ['norm']:
        distribution = getattr(st, distribution_name)

        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                params = distribution.fit(data)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution.pdf(x, loc=loc, scale=scale, *arg)
                sse = np.sum(np.power(y - pdf, 2.0))

                # identify if this distribution is better
                best_distributions.append((distribution, params, sse, distribution_name))

        except Exception:
            pass

    best_dist = sorted(best_distributions, key=lambda x: x[2])[0]
    dist = best_dist[0]
    params =  best_dist[1]

    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf, best_dist[-1]

def EPDF(data):
    hist, bins = np.histogram(np.nan_to_num(data), bins=20, density=True)
    bin_center = 0.5 * (bins[1:] + bins[:-1])
    return bin_center, hist


def KM(data, n_clusters, centroids=False, reshape_dim=[56,76]):
    with threadpool_limits(limits=1):#con 1 mas rapido q con 40...
        kmeans = KMeans(n_clusters=n_clusters, n_init=100, random_state=666).fit(data)

    y_pred = kmeans.predict(data)
    if centroids:
        return y_pred.reshape(reshape_dim[0], reshape_dim[1]).T, kmeans.cluster_centers_
    else:
        return y_pred.reshape(reshape_dim[0], reshape_dim[1]).T

def MakeMask(DataArray, dataname='mask'):
    import regionmask
    mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(DataArray)
    mask = xr.where(np.isnan(mask), mask, 1)
    mask = mask.to_dataset(name=dataname)
    return mask

def CFSv2Data(v):
    out_data_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
    dir_hc = '/pikachu/datos/osman/nmme/monthly/hindcast/'
    try:
        data_cfsv2 = xr.open_dataset(out_data_dir + v + '_data_cfsv2_noRoll.nc')
    except:
        def fix_calendar(ds, timevar='time'):
            """
            agrega los dias a los archivos nc de NMME
            """
            if ds[timevar].attrs['calendar'] == '360':
                ds[timevar].attrs['calendar'] = '360_day'
            return ds

        from ENSO_IOD_Funciones import SelectNMMEFiles
        files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable=v,
                                dir=dir_hc, All=True)
        files = sorted(files, key=lambda x: x.split()[0])

        # abriendo todos los archivos
        data_cfsv2 = xr.open_mfdataset(files, decode_times=False).sel(
            L=[0.5, 1.5, 2.5, 3.5])  # xr no entiende la codificacion de Leads, r y las fechas
        data_cfsv2 = data_cfsv2.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
        data_cfsv2['L'] = [0, 1, 2, 3]
        data_cfsv2 = xr.decode_cf(fix_calendar(data_cfsv2))  # corrigiendo fechas
        data_cfsv2 = data_cfsv2.sel(lon=slice(275, 330), lat=slice(-60, 15))
        data_cfsv2 = data_cfsv2.sel(r=slice(1, 24))
        data_cfsv2.to_netcdf(out_data_dir + v + '_data_cfsv2_noRoll.nc')
        data_cfsv2 = data_cfsv2.compute()
    return data_cfsv2

def KmMaks(data_cfsv2, v, l, pb):
    data2 = data_cfsv2.sel(L=l)
    # data2 = data2.sel(lon=slice(290, 330), lat=slice(-60, 0))
    data_clim = data2.groupby('time.month').mean()
    data_clim = data_clim.mean(['r'])
    data_clim = data_clim.rename({'month': 'time'})
    data_clim.load()
    data_clim *= 30
    # data_clim *= mask
    # if l == 0:
    #     pb = 0
    # else:
    #     pb = 2

    pp2 = data_clim.sel(lon=slice(pre_box_lon[pb][0], pre_box_lon[pb][1]),
                        lat=slice(pre_box_lat[pb][0], pre_box_lat[pb][1]))
    #### k-means #####
    X = pp2.stack(new=['lon', 'lat'])[v].T
    X[np.where(np.isnan(X))] = -99

    for n_c in [7]:
        reshape_lat_dim = len(range(pre_box_lat[pb][0], pre_box_lat[pb][1])) + 1
        reshape_lon_dim = len(range(pre_box_lon[pb][0], pre_box_lon[pb][1])) + 1

        clusters, centroids = KM(data=X, n_clusters=n_c, centroids=True,
                                 reshape_dim=[reshape_lon_dim, reshape_lat_dim])

        xr_cluster_label = MakeMask(pp2, dataname='cluster')
        xr_cluster_label['cluster'].values = clusters

        xr_cluster_label *= MakeMask(xr_cluster_label, 'cluster')

    return xr_cluster_label
########################################################################################################################
# Region seleccionada con k-means a paritr del CFSv2 ------------------------------------------------------------------#
pre_box = [0, 1, 2]
pre_box_lat = [[-60, 0], [-60, -20], [-60, -20]]
pre_box_lon = [[290, 330], [290, 330], [280, 290]]

seasons = ['JJA', 'JAS', 'ASO', 'SON']
cases = ['dmi_puros_pos', 'dmi_puros_neg',
        'n34_puros_pos', 'n34_puros_neg',
        'sim_pos', 'sim_neg',
        'dmi_neg_n34_pos', 'dmi_pos_n34_neg',
        'neutros']
#cases2 = ['dmi_puros_pos', 'dmi_puros_neg', 'n34_puros_pos', 'n34_puros_neg', 'sim_pos', 'sim_neg']

positive_cases = ['dmi_puros_pos', 'n34_puros_pos', 'sim_pos']
negative_cases = ['dmi_puros_neg', 'n34_puros_neg', 'sim_neg']

title_case = ['DMI pure - positive',
              'DMI pure - negative',
              'El Niño pure', 'La Niña pure',
              'DMI positive - El Niño',
              'DMI negative - La Niña']
# cajas ---------------------------------------------------------------------------------------------------------------#
box_name = ['S_SESA', 'SESA', 'N_SESA', 'C_Brazil', 'Chile', 'Chile_sur']
box_color =['Spectral_r', 'afmhot_r', 'BrBG', 'BrBG_r', 'RdBu', 'RdBu_r', 'Accent']
box_lats = [[-35,-29],[-35,-20],[-29,-20],[-20,0],[-40,-30], [-60,-43]]
box_lons = [[300, 310],[300,320],[300,320],[300,325],[285,290], [285,288]]
#---------------------------------------------------------------------------------------------------------------------#
v='prec'

data_cfsv2 = CFSv2Data(v)
# mask = OpenDataSet('pp_gpcc', interp=True,
#                    lat_interp=np.linspace(-60,15,76),
#                    lon_interp=np.linspace(275,330,56))


fix_factor = [30,1]
startend = [80,4]
cluster_set = [0,1]
cluster_name = ['SESA', 'outSESA']
units = ['[mm]', '[ºC]']

ylim = [[0.020,0.06], [0.7,0.7]]


xr_cluster_label = KmMaks(data_cfsv2, v=v, l=0, pb=0)
mask_cluster = xr.where(xr_cluster_label == 3, 1, np.nan).sel(lat=slice(-40, -20), lon=slice(290, 320))
mask_cluster = mask_cluster.rename({'cluster':'var'})
v_count=0
for v in ['prec','tref']:
    for s in seasons:
        for cs in cluster_set:
            if cs == 0:
                mask_cluster = xr.where(xr_cluster_label == 3, 1, np.nan).sel(lat=slice(-40, -20), lon=slice(290, 320))
            else:
                mask_cluster = xr.where((xr_cluster_label == 2) | (xr_cluster_label == 0), 1, np.nan) \
                    .sel(lat=slice(-40, -10), lon=slice(290, 320))

            mask_cluster = mask_cluster.rename({'cluster': 'var'})

            # Climatologia
            neutro = xr.open_dataset(cases_dir + v + '_neutros_' + s + '.nc').rename({v: 'var'})
            neutro *= fix_factor[v_count]
            neutro *= mask_cluster
            # neutro = neutro.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
            # neutro = neutro.drop(['r', 'L'])

            c_count = 0
            for c in cases:
                case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s + '.nc').rename({v: 'var'})
                case *= mask_cluster
                case *= fix_factor[v_count]
                # case = case.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
                # case = case.drop(['L', 'r'])

                # creando la climatologia con TODOS los eventos (esto es mas rapido que hacer otro SELECT_EVENTS
                if c_count == 0:
                    clim_full = xr.concat([neutro, case], dim='time')
                    c_count = 1
                else:
                    clim_full = xr.concat([clim_full, case], dim='time')

            clim_full = clim_full.mean(['lon', 'lat'])

            fig, ax = plt.subplots(1, 2, figsize=(11, 4.5), dpi=dpi)
            ax[0].set_xlim(-startend[v_count], startend[v_count])
            ax[1].set_xlim(-startend[v_count], startend[v_count])
            ax[0].set_ylim(0, ylim[v_count][cs])
            ax[1].set_ylim(0, ylim[v_count][cs])
            ax[0].set_xlabel(v + '- anomaly ' + units[v_count])
            ax[1].set_xlabel(v + '- anomaly ' + units[v_count])

            ax[1].set_title('Clim. - Negative Phases', color='royalblue')
            ax[0].set_title('Clim - Positive Phases', color='firebrick')

            cases_color = ['darkorange', 'blue', 'firebrick', 'cyan', 'forestgreen', 'purple']
            positive_cases_colors = ['red', '#F47B00', 'forestgreen']
            negative_cases_colors = ['blue', '#55F400', 'magenta']

            aux_clim_full = clim_full - clim_full.mean('time')
            aux_clim_full = aux_clim_full['var'].values
            # epdf_clim_full_x, epdf_clim_full = EPDF(aux_clim_full)
            pdf_clim_full, name_clim_full_dist = best_fit_distribution(np.nan_to_num(aux_clim_full),
                                                                       len(aux_clim_full),
                                                                       -1 * startend[v_count], startend[v_count])
            # ax[0].hist(aux_clim_full, bins=30, density=True, color='k', alpha=0.3)
            ax[0].plot(pdf_clim_full, lw=1.5, color='k', label='climatology')

            # ax[1].hist(aux_clim_full, bins=30, density=True, color='k', alpha=0.3)
            ax[1].plot(pdf_clim_full, lw=1.5, color='k', label='climatology')

            c_count = 0
            for c in positive_cases:
                case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s + '.nc').rename({v: 'var'})
                case *= mask_cluster
                case *= fix_factor[v_count]
                case = case.mean(['lon', 'lat'])

                # signal
                case_anom = case - clim_full.mean('time')
                case_anom = case_anom['var'].values
                # epdf_case_anom_x, epdf_case_anom = EPDF(case_anom)
                pdf_case, name_case_dist = best_fit_distribution(np.nan_to_num(case_anom),
                                                                 len(case_anom),
                                                                 -1 * startend[v_count], startend[v_count])

                ax[0].hist(case_anom, bins=30, density=True, color=positive_cases_colors[c_count], alpha=0.3)
                ax[0].plot(pdf_case, lw=2.5, color=positive_cases_colors[c_count], label=c)
                c_count += 1

            c_count = 0
            for c in negative_cases:
                case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s + '.nc').rename({v: 'var'})
                case *= mask_cluster
                case *= fix_factor[v_count]
                case = case.mean(['lon', 'lat'])

                # signal
                case_anom = case - clim_full.mean('time')
                case_anom = case_anom['var'].values
                # epdf_case_anom_x, epdf_case_anom = EPDF(case_anom)
                pdf_case, name_case_dist = best_fit_distribution(np.nan_to_num(case_anom),
                                                                 len(case_anom),
                                                                 -1 * startend[v_count], startend[v_count])

                ax[1].hist(case_anom, bins=30, density=True, color=negative_cases_colors[c_count], alpha=0.3)
                ax[1].plot(pdf_case, lw=2.5, color=negative_cases_colors[c_count], label=c)
                c_count += 1

            ax[0].grid(alpha=0.5)
            ax[1].grid(alpha=0.5)
            ax[0].legend(loc='best')
            ax[1].legend(loc='best')
            # ax[1][0].legend()
            # ax[1][1].legend()

            fig.suptitle(v + ' - ' + cluster_name[cs] + ' - ' + s)
            # fig.legend(bbox_to_anchor=(1, 0.5), loc="center right")  # , bbox_transform=fig.transFigure, ncol=len(cases2))
            # plt.subplots_adjust(right=0.85)
            plt.yticks(size=10)
            plt.xticks(size=10)
            # plt.tight_layout(rect=[0, 0, 0.7, 0])
            if save:
                plt.savefig(out_dir + v + '_PDF_' + cluster_name[cs] + '_' + s + '.jpg', box_inches="tight")
            else:
                plt.show()

    v_count += 1

########################################################################################################################