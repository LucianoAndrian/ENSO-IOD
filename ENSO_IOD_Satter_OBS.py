"""
Scatter obs
"""
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from ENSO_IOD_Funciones import DMI2, Nino34CPC
###############################################################################
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/paper1/1940_2020/'

save = True
dpi = 300

################################################################################
def auxScatter(n34, n34_3, dmi, dmi_3, s):
    dmi_todos = dmi_3.sel(time=dmi_3.time.dt.month.isin([s]))
    dmi_criteria_y = dmi.where((dmi.Mes == s)).Años.dropna().values

    n34_todos = n34.sel(time=n34.time.dt.month.isin([s]))
    n34_criteria_y = n34_3.where((n34_3.Mes == s)).Años.dropna().values

    sim_y = np.intersect1d(n34_criteria_y, dmi_criteria_y)

    dmi_sim = dmi_todos.sel(time=dmi_todos.time.dt.year.isin(sim_y))
    n34_sim = n34_todos.sel(time=n34_todos.time.dt.year.isin(sim_y))

    dmi_sim_pos = dmi_sim.where(dmi_sim > 0)
    n34_sim_pos = n34_sim.where(n34_sim > 0)

    dmi_sim_neg = dmi_sim.where(dmi_sim < 0)
    n34_sim_neg = n34_sim.where(n34_sim < 0)

    dmi_pos_n34_neg = dmi_sim_pos.where(~np.isnan(n34_sim_neg.values))
    dmi_neg_n34_pos = dmi_sim_neg.where(~np.isnan(n34_sim_pos.values))

    dmi_dates_ref = dmi_todos.time.dt.year
    mask = np.in1d(dmi_dates_ref, dmi_criteria_y)
    aux_dmi = dmi_todos.sel(time=dmi_todos.time.dt.year.isin(dmi_dates_ref[mask]))

    n34_dates_ref = n34_todos.time.dt.year
    mask = np.in1d(n34_dates_ref, n34_criteria_y)
    aux_n34 = n34_todos.sel(time=n34_todos.time.dt.year.isin(n34_dates_ref[mask]))

    aux_dates_ref = aux_dmi.time.dt.year
    mask = np.in1d(aux_dates_ref, sim_y, invert=True)
    dmi_un = aux_dmi.sel(time=aux_dmi.time.dt.year.isin(aux_dates_ref[mask]))

    dmi_un_pos = dmi_un.where(dmi_un > 0)
    dmi_un_pos_n34_values = n34_todos.sel(time=n34_todos.time.isin(dmi_un_pos.time))

    dmi_un_neg = dmi_un.where(dmi_un < 0)
    dmi_un_neg_n34_values = n34_todos.sel(time=n34_todos.time.isin(dmi_un_neg.time))

    aux_dates_ref = aux.time.dt.year
    mask = np.in1d(aux_dates_ref, sim_y, invert=True)
    n34_un = aux_n34.sel(time=aux_n34.time.dt.year.isin(aux_dates_ref[mask]))

    n34_un_pos = n34_un.where(n34_un > 0)
    n34_un_pos_dmi_values = dmi_todos.sel(time=dmi_todos.time.isin(n34_un_pos.time))
    n34_un_neg = n34_un.where(n34_un < 0)
    n34_un_neg_dmi_values = dmi_todos.sel(time=dmi_todos.time.isin(n34_un_neg.time))

    return dmi_un_pos, dmi_un_pos_n34_values, dmi_un_neg, dmi_un_neg_n34_values, \
           n34_un_pos, n34_un_pos_dmi_values, n34_un_neg, n34_un_neg_dmi_values, \
           dmi_sim_pos, n34_sim_pos, dmi_sim_neg, n34_sim_neg, dmi_todos, \
           n34_todos, dmi_pos_n34_neg, dmi_neg_n34_pos
################################################################################
i = 1920
end = 2020
seasons = ['SON']
seasons_n = [10]
################################################################################
# indices: ---------------------------------------------------------------------
dmi, dmi_2, dmi_3 = DMI2(filter_bwa=False, start_per=str(i), end_per=str(end),
                         sst_anom_sd=False, opposite_signs_criteria=False)

dmi_td, dmi_2, dmi_3_td = DMI2(filter_bwa=False, start_per=str(i), end_per=str(end),
                         sst_anom_sd=True, opposite_signs_criteria=True)

aux = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
n34, n34_2, n34_3 = Nino34CPC(aux, start=i)

dmi_3 = dmi_3.sel(time=slice('1940-01-01', '2020-12-01'))
dmi = dmi.where(dmi.Años>=1940).dropna()
dmi_3_td = dmi_3_td.sel(time=slice('1940-01-01', '2020-12-01'))
dmi_td = dmi_td.where(dmi_td.Años>=1940).dropna()

n34 = n34.sel(time=slice('1940-01-01', '2020-12-01'))
n34_3 = n34_3.where(n34_3.Años>=1940).dropna()

dmi_3 = dmi_3 / dmi_3.std('time')
dmi_3_td = dmi_3_td / dmi_3_td.std('time')
n34 = n34 / n34.std('time')

for s in seasons_n:
    dmi_un_pos, dmi_un_pos_n34_values, dmi_un_neg, dmi_un_neg_n34_values, \
    n34_un_pos, n34_un_pos_dmi_values, n34_un_neg, n34_un_neg_dmi_values, \
    dmi_sim_pos, n34_sim_pos, dmi_sim_neg, n34_sim_neg, dmi_todos, n34_todos, \
    dmi_pos_n34_neg, dmi_neg_n34_pos = auxScatter(n34, n34_3, dmi, dmi_3, s)

    dmi_un_pos_td, dmi_un_pos_n34_values_td, dmi_un_neg_td, dmi_un_neg_n34_values_td, \
    n34_un_pos_td, n34_un_pos_dmi_values_td, n34_un_neg_td, n34_un_neg_dmi_values_td, \
    dmi_sim_pos_td, n34_sim_pos_td, dmi_sim_neg_td, n34_sim_neg_td, dmi_todos, n34_todos, \
    dmi_pos_n34_neg_td, dmi_neg_n34_pos_td = auxScatter(n34, n34_3, dmi_td, dmi_3_td, s)
    #---------------------------------------------------------------------------

    fig, ax = plt.subplots(dpi=dpi)
    # todos
    plt.scatter(x=dmi_todos, y=n34_todos, marker='.', label='Niño3.4 vs DMI',
                s=20, edgecolor='k', color='dimgray', alpha=1)
    # dmi puros
    plt.scatter(x=dmi_un_pos.values, y=dmi_un_pos_n34_values.values, marker='>',
                s=70, edgecolor='firebrick', facecolor='firebrick', alpha=1,
                label='IOD+, NO ENSO')
    plt.scatter(x=dmi_un_neg.values, y=dmi_un_neg_n34_values.values, marker='<',
                s=70, facecolor='limegreen', edgecolor='limegreen', alpha=1,
                label='IOD-, NO ENSO')
    # n34 puros
    plt.scatter(y=n34_un_pos.values, x=n34_un_pos_dmi_values.values, marker='^',
                s=70, edgecolors='navy', facecolor='navy', alpha=1,
                label='Niña, NO IOD')
    plt.scatter(y=n34_un_neg.values, x=n34_un_neg_dmi_values.values, marker='v',
                s=70, edgecolors='deeppink', facecolor='deeppink', alpha=1,
                label='Niña, NO IOD')
    # sim
    plt.scatter(x=dmi_sim_pos.values, y=n34_sim_pos.values, marker='s', s=50,
                edgecolor='red',color='red', alpha=1, label='Niño & IOD+')
    plt.scatter(x=dmi_sim_neg.values, y=n34_sim_neg.values, marker='s', s=50,
                edgecolor='deepskyblue', color='deepskyblue', alpha=1,
                label='Niña & IOD-')
    # sim opp. sing
    plt.scatter(x=dmi_pos_n34_neg.values, y=n34_sim_neg.values, marker='s', s=50,
                edgecolor='orange',color='orange', alpha=1, label='Niña & IOD+')
    plt.scatter(x=dmi_neg_n34_pos.values, y=n34_sim_neg.values, marker='s', s=50,
                edgecolor='gold', color='gold', alpha=1, label='Niño & IOD-')

    plt.scatter(x=dmi_sim_pos_td.values, y=n34_sim_pos_td.values, marker='+',
                s=70,
                color='k', alpha=1)
    plt.scatter(x=dmi_sim_neg_td.values, y=n34_sim_neg_td.values, marker='+',
                s=70,
                color='k', alpha=1)
    plt.scatter(x=dmi_pos_n34_neg_td.values, y=n34_sim_neg_td.values, marker='+',
                s=70,
                color='k', alpha=1)
    plt.scatter(x=dmi_neg_n34_pos_td.values, y=n34_sim_neg_td.values, marker='+',
                s=70,
                color='k', alpha=1)
    plt.scatter(x=dmi_un_pos_td.values, y=dmi_un_pos_n34_values_td.values,
                marker='+',
                s=70, color='k', alpha=1)
    plt.scatter(x=dmi_un_neg_td.values, y=dmi_un_neg_n34_values_td.values,
                marker='+',
                s=70, color='k', alpha=1)

    plt.legend(loc=(.01,.55))

    plt.ylim((-5, 5))
    plt.xlim((-5, 5))
    plt.axhspan(-.5, .5, alpha=0.2, color='black', zorder=0)
    plt.axvspan(-.5, .5, alpha=0.2, color='black', zorder=0)
    # ax.grid(True)
    fig.set_size_inches(6, 6)
    plt.xlabel('IOD', size=15)
    plt.ylabel('Niño 3.4', size=15)

    plt.text(-4.8, 4.6, 'EN/IOD-', dict(size=15))
    plt.text(-.3, 4.6, 'EN', dict(size=15))
    plt.text(+3.1, 4.6, 'EN/IOD+', dict(size=15))
    plt.text(+3.8, -.1, 'IOD+', dict(size=15))
    plt.text(+3.1, -4.9, 'LN/IOD+', dict(size=15))
    plt.text(-.3, -4.9, 'LN', dict(size=15))
    plt.text(-4.8, -4.9, ' LN/IOD-', dict(size=15))
    plt.text(-4.8, -.1, 'IOD-', dict(size=15))
    plt.title('SON' + ' - ' + 'OBS')
    plt.tight_layout()
    if save:
        plt.savefig(out_dir + 'ENSO_IOD_Scatter_comparision_dmi_SON_OBS_oneTrue.jpg')
    else:
        plt.show()
########################################################################################################################