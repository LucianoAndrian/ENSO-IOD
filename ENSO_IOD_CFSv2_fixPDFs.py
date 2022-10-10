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
#warnings.filterwarnings( "ignore", module = "matplotlib\..*" )
warnings.filterwarnings('ignore')
import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
import matplotlib.pyplot as plt
########################################################################################################################
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/PDFs/'
save=True
# Funciones ############################################################################################################

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
    hist, bins = np.histogram(np.nan_to_num(data), bins=30, density=True)
    bin_center = 0.5 * (bins[1:] + bins[:-1])
    return bin_center, hist
########################################################################################################################

seasons = ['JJA', 'JAS', 'ASO', 'SON']
cases = ['dmi_puros_pos', 'dmi_puros_neg', 'n34_puros_pos', 'n34_puros_neg', 'sim_pos', 'sim_neg']

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



mask = OpenDataSet('pp_gpcc', interp=True,
                   lat_interp=np.linspace(-60,15,76),
                   lon_interp=np.linspace(275,330,56))

#for v in variables
v = 'prec'
fix_factor = 30
startend = 80

# s = 'SON'
# b_count = 0
# b = box_name[0]
#
# c_count=0
# c = cases[0]

for s in seasons:
    b_count = 0
    for b in box_name:

        min_lat, max_lat = box_lats[b_count]
        min_lon, max_lon = box_lons[b_count]

        # Climatologia
        neutro = xr.open_dataset(cases_dir + v + '_neutros_' + s + '.nc').rename({v:'var'})
        neutro *= fix_factor
        neutro *= mask
        neutro = neutro.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
        neutro = neutro.drop(['r', 'L'])

        c_count=0
        for c in cases:
            case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s + '.nc').rename({v:'var'})
            case *= fix_factor
            case = case.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
            case = case.drop(['L', 'r'])

            # creando la climatologia con TODOS los eventos (esto es mas rapido que hacer otro SELECT_EVENTS
            if c_count == 0:
                clim_full = xr.concat([neutro, case], dim='time')
                c_count = 1
            else:
                clim_full = xr.concat([clim_full, case], dim='time')

        clim_full = clim_full.mean(['lon', 'lat'])


        fig, ax = plt.subplots(1, 2, figsize=(11, 4), dpi=300)
        ax[0].set_xlim(-startend, startend)
        ax[1].set_xlim(-startend, startend)

        ax[0].set_title('Clim. - Empirical')
        ax[1].set_title('Clim - Theoretical')

        cases_color = ['orange', 'blue', 'firebrick', 'cyan', 'forestgreen', 'purple']

        c_count = 0
        for c in cases:
            case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s + '.nc').rename({v:'var'})
            case *= fix_factor
            case = case.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
            case = case.drop(['L', 'r'])
            case = case.mean(['lon', 'lat'])

            # signal
            aux_clim_full = clim_full - clim_full.mean('time')
            case_anom = case - clim_full.mean('time')

            aux_clim_full = aux_clim_full['var'].values

            epdf_clim_full_x, epdf_clim_full = EPDF(aux_clim_full)
            pdf_clim_full, name_clim_full_dist = best_fit_distribution(np.nan_to_num(aux_clim_full),
                                                                       len(aux_clim_full),
                                                                       -1 * startend, startend)

            case_anom = case_anom['var'].values
            epdf_case_anom_x, epdf_case_anom = EPDF(case_anom)
            pdf_case, name_case_dist = best_fit_distribution(np.nan_to_num(case_anom),
                                                                       len(case_anom),
                                                                       -1 * startend, startend)

            ax[0].plot(epdf_clim_full_x, epdf_clim_full, lw=2, color='k')
            ax[0].plot(epdf_case_anom_x, epdf_case_anom, lw=2, label=c, color=cases_color[c_count])

            ax[1].plot(pdf_clim_full, lw=2, color='k')
            ax[1].plot(pdf_case, lw=2,
                       color=cases_color[c_count])
            c_count += 1

        ax[0].grid(alpha=0.5)
        ax[1].grid(alpha=0.5)
        # ax[1][0].legend()
        # ax[1][1].legend()

        fig.suptitle(v + ' - ' + b + ' - ' + s)
        fig.legend(bbox_to_anchor=(1, 0.5), loc="center right")  # , bbox_transform=fig.transFigure, ncol=len(cases2))
        plt.subplots_adjust(right=0.85)
        plt.yticks(size=10)
        plt.xticks(size=10)
        plt.tight_layout(rect=[0, 0, 0.7, 0])
        if save:
            plt.savefig(out_dir + v + '_PDF_' + b + '_' + s + '.jpg', box_inches="tight")
        else:
            plt.show()

        b_count += 1
########################################################################################################################