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
cases_dir = '/pikachu/datos/luciano.andrian/cases/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/PDFs/'
save=True
dpi=200
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
    hist, bins = np.histogram(np.nan_to_num(data), bins=20, density=True)
    bin_center = 0.5 * (bins[1:] + bins[:-1])
    return bin_center, hist

########################################################################################################################
variables = ('prec', 'tref')
seasons = ('JJA', 'JAS', 'ASO', 'SON')
set_name= ('3')

cases = ['sim_pos', 'sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']#, 'sim_DMIneg_N34pos', 'sim_DMIpos_N34neg']

cases2 = ['sim_pos', 'sim_neg', 'DMI_un_pos', 'DMI_un_neg',
         'N34_un_pos', 'N34_un_neg']#, 'sim_DMIneg_N34pos', 'sim_DMIpos_N34neg']

title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI negative phase ',
              'DMI positive phase ',
              'DMI isolated positive phase ',
              'DMI isolated negative phase ',
              'ENSO positive phase ',
              'ENSO negative phase ',
              'ENSO isolated positive phase ',
              'ENSO isolated negative phase ',
              'DMI negative and ENSO positive',
              'DMI positive and ENSO negative']

box_name = ['S_SESA', 'S_SESA_exp', 'N_SESA', 'C_Brazil', 'Chile']
box_color =['Spectral_r', 'Spectral', 'BrBG', 'BrBG_r', 'RdBu', 'RdBu_r']
box_lats = [[-35,-29],[-35,-25],[-29,-20],[-20,0],[-40,-30]]
box_lons = [[300, 310],[295,315],[300,320],[300,325],[285,290]]


maxmins=[60,2]

mask = OpenDataSet('pp_gpcc', interp=True,
                   lat_interp=np.linspace(-60,15,76),
                   lon_interp=np.linspace(275,330,56))


v = 'prec'
fix_factor=30
s='SON'
c='N34_un_neg'
s_n= '3'
v_count = 0
for v in variables:
    if v == 'prec':
        fix_factor = 30
        startend=60
    else:
        fix_factor=1
        startend = 3
    for s in seasons:
        print(s)
        for s_n in set_name:
            print(s_n)

            b_count=0
            for b_name in box_name:

                min_lat, max_lat = box_lats[b_count]
                min_lon, max_lon = box_lons[b_count]

                # Climatologia
                data_neutral = xr.open_dataset(cases_dir + v + '_' + s + '_Set' + s_n + '_NEUTRO.nc').drop(['L', 'r'])
                data_neutral = data_neutral.rename({v: 'var'})
                data_neutral *= fix_factor
                data_neutral *= mask
                data_neutral = data_neutral.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
                c2_count = 0
                for c2 in cases2:
                    print(c2)
                    data_case = xr.open_dataset(cases_dir + v + '_' + s + '_Set' + s_n + '_' + c2 + '.nc').drop(
                        ['L', 'r'])
                    data_case = data_case.rename({v: 'var'})
                    data_case *= fix_factor
                    data_case = data_case.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))

                    # creando la climatologia con TODOS los eventos (esto es mas rapido que hacer otro SELECT_EVENTS
                    if c2_count == 0:
                        data_clim_full = xr.concat([data_neutral, data_case], dim='time')
                        c2_count = 1
                    else:
                        data_clim_full = xr.concat([data_clim_full, data_case], dim='time')

                data_clim_full = data_clim_full.mean(['lon', 'lat'])
                data_neutral = data_neutral.mean(['lon', 'lat'])

                c_count = 0
                fig, ax = plt.subplots(1, 2, figsize=(11, 4), dpi=300)
                ax[0].set_xlim(-startend, startend)
                ax[1].set_xlim(-startend, startend)
                # ax[1][0].set_xlim(-startend, startend)
                # ax[1][1].set_xlim(-startend, startend)

                ax[0].set_title('Clim. - Empirical')
                ax[1].set_title('Clim - Theorical')
                # ax[1][0].set_title('Clim. Neutral - Theorical')
                # ax[1][1].set_title('Clim. Full - Theorical')

                cases_color = ['gold', 'dodgerblue', 'firebrick', 'cyan', 'forestgreen', 'purple']
                for c in cases2:
                    print(c)
                    data_case = xr.open_dataset(cases_dir + v + '_' + s + '_Set' + s_n + '_' + c + '.nc').drop(
                        ['L', 'r'])
                    data_case = data_case.rename({v: 'var'})
                    data_case *= fix_factor
                    data_case = data_case.sel(lat=slice(min_lat, max_lat), lon=slice(min_lon, max_lon))
                    data_case = data_case.mean(['lon', 'lat'])

                    # signal
                    clim_neutral = data_neutral - data_neutral.mean('time')
                    clim_full = data_clim_full - data_clim_full.mean('time')
                    case_neutral_anom = data_case - data_neutral.mean('time')
                    case_full_anom = data_case - data_clim_full.mean('time')


                    # PDFs
                    # clim_neutral = clim_neutral['var'].values
                    # epdf_clim_neutral_x, epdf_clim_neutral = EPDF(clim_neutral)
                    # pdf_clim_neutral, name_clim_neutral_dist = best_fit_distribution(np.nan_to_num(clim_neutral),
                    #                                                                  len(clim_neutral),
                    #                                                                  -1*startend, startend)

                    clim_full = clim_full['var'].values
                    epdf_clim_full_x, epdf_clim_full = EPDF(clim_full)
                    pdf_clim_full, name_clim_full_dist = best_fit_distribution(np.nan_to_num(clim_full),
                                                                               len(clim_full),
                                                                               -1*startend, startend)

                    # case_neutral_anom = case_neutral_anom['var'].values
                    # epdf_case_neutral_anom_x, epdf_case_neutral_anom = EPDF(case_neutral_anom)
                    # pdf_case_neutral, name_case_neutral_dist = best_fit_distribution(np.nan_to_num(case_neutral_anom),
                    #                                                                  len(case_neutral_anom),
                    #                                                                  -1*startend, startend)

                    case_full_anom = case_full_anom['var'].values
                    epdf_case_full_anom_x, epdf_case_full_anom = EPDF(case_full_anom)
                    pdf_case_full, name_case_full_dist = best_fit_distribution(np.nan_to_num(case_full_anom),
                                                                               len(case_full_anom),
                                                                               -1*startend, startend)

                    # ax[0][0].plot(epdf_clim_neutral_x, epdf_clim_neutral, lw=2, color='k')
                    # ax[0][0].plot(epdf_case_neutral_anom_x, epdf_case_neutral_anom, lw=2, label=c, color=cases_color[c_count])

                    ax[0].plot(epdf_clim_full_x, epdf_clim_full, lw=2, color='k')
                    ax[0].plot(epdf_case_full_anom_x, epdf_case_full_anom, lw=2, label=c, color=cases_color[c_count])

                    # ax[1][0].plot(pdf_clim_neutral, lw=2, color='k')
                    # ax[1][0].plot(pdf_case_neutral, lw=2, label=c,
                    #            color=cases_color[c_count])

                    ax[1].plot(pdf_clim_full, lw=2, color='k')
                    ax[1].plot(pdf_case_full, lw=2,
                               color=cases_color[c_count])
                    c_count +=1


                ax[0].grid(alpha=0.5)
                ax[1].grid(alpha=0.5)
                # ax[1][0].legend()
                # ax[1][1].legend()

                fig.suptitle(v + ' - ' + b_name + ' - ' + s)
                fig.legend(bbox_to_anchor=(1, 0.5), loc="center right")#, bbox_transform=fig.transFigure, ncol=len(cases2))
                plt.subplots_adjust(right=0.85)
                plt.yticks(size=10)
                plt.xticks(size=10)
                plt.tight_layout(rect=[0, 0, 0.7, 0])
                if save:
                    plt.savefig(out_dir + v + '_PDF_' + b_name + '_' + s + '.jpg',box_inches = "tight")
                else:
                    plt.show()

                b_count += 1
