"""
Anomalía de PP (y T cuando ande) en regiones de SA según la magnitud de los índices
3er intento
"""
########################################################################################################################
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os
from matplotlib import colors
from sklearn.cluster import KMeans
from threadpoolctl import threadpool_limits
from ENSO_IOD_Funciones import SelectNMMEFiles
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
dates_dir = '/pikachu/datos/luciano.andrian/DMI_N34_Leads_r/'
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/' # campos de las variables PP  ( y T cuando ande)
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/RegionsAnoms/'
save = True
dpi = 300
# Funciones ############################################################################################################
def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds

def SelectBins(data, min, max, sd=1):
    # sd opcional en caso de no estar escalado
    # En este caso los >= y <=0 son necesarios siempre ya que se elecciona en funcion de los cases,
    # es decir en los dmi entre -0.5 y 0.5 de N34, TODOS son dmi_puros_positivos

    if np.abs(min) > np.abs(max):
        return (data >= min*sd) & (data <= max*sd)
    elif np.abs(min) < np.abs(max):
        return (data >= min*sd) & (data <= max*sd)
    elif np.abs(min) == np.abs(max):
        return (data > min*sd) & (data < max*sd)

def SelectVariables(dates, data):
    t_count=0
    for t in dates.index:
        r_t = t.r.values
        L_t = t.L.values
        t_t = t.values
        try: #q elegancia la de francia...
            t_t*1
            t_t = t.time.values
        except:
            pass

        if t_count == 0:
            aux = data.where(data.L == L_t).sel(r=r_t, time=t_t)
            t_count += 1
        else:
            aux = xr.concat([aux,
                             data.where(data.L == L_t).sel(r=r_t, time=t_t)],
                            dim='time')
    return aux

def BinsByCases(v, v_name, fix_factor, s, mm, c, c_count, box_lat, box_lon,
                bin_limits, bins_by_cases_dmi, bins_by_cases_n34,
                mask=False):

    # 1. se abren los archivos de los índices (completos y se pesan por su SD)
    # tambien los archivos de la variable en cuestion pero para cada "case" = c

    data_dates_dmi_or = xr.open_dataset(dates_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc')
    data_dates_dmi_or /= data_dates_dmi_or.mean('r').std() #corregido!

    data_dates_n34_or = xr.open_dataset(dates_dir + 'N34_' + s + '_Leads_r_CFSv2.nc')
    data_dates_n34_or /= data_dates_n34_or.mean('r').std()

    # 1.1 Climatología y case
    if v == 'tref':
        end_nc_file = '_nodetrend.nc'
    else:
        end_nc_file = '.nc'

    clim = xr.open_dataset(cases_dir + v + '_' + s.lower() + end_nc_file).rename({v_name: 'var'}) * fix_factor
    case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s.upper() + end_nc_file).rename({v_name: 'var'}) * fix_factor

    # Anomalía
    for l in [0, 1, 2, 3]:
        if l == 0:
            anom = case.sel(time=case.time.dt.month.isin(mm - l)) - \
                   clim.sel(time=clim.time.dt.month.isin(mm - l)).mean(['r', 'time'])
        else:
            anom2 = case.sel(time=case.time.dt.month.isin(mm - l)) - \
                    clim.sel(time=clim.time.dt.month.isin(mm - l)).mean(['r', 'time'])
            anom = xr.concat([anom, anom2], dim='time')

    # 1.2 Se selecciona la región
    if mask:
        anom = anom * mask
    anom = anom.sel(lat=slice(box_lat[0], box_lat[1]),
                    lon=slice(box_lon[0], box_lon[1]))
    BOX = anom.mean(['lon', 'lat'])
    BOX = BOX.sortby(BOX.time.dt.month)
    # if mask:
    #     # PROBANDO!!!
    #     BOX = anom*data_mask
    #     #BOX = BOX.sel(lat=slice(min(box_lat), max(box_lat)), lon=slice(min(box_lon), max(box_lon)))
    #     BOX = BOX.mean(['lon', 'lat'])
    #     BOX = BOX.sortby(BOX.time.dt.month)
    # else:
    #     BOX = anom.sel(lat=slice(min(box_lat), max(box_lat)), lon=slice(min(box_lon), max(box_lon)))
    #     BOX = BOX.mean(['lon', 'lat'])
    #     BOX = BOX.sortby(BOX.time.dt.month)

    # 2. Vinculo fechas case -> índices DMI y N34 para poder clasificarlos
    # las fechas entre el case variable y el case indices COINCIDEN,
    # DE ESA FORMA SE ELIGIERON LOS CASES VARIABLE
    # pero diferen en orden. Para evitar complicar la selección usando r y L
    # con .sortby(..time.dt.month) en cada caso se simplifica el problema
    # y coinciden todos los eventos en fecha, r y L

    cases_date_dir = '/pikachu/datos/luciano.andrian/cases_dates/'
    try:
        aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
            .rename({'__xarray_dataarray_variable__': 'index'})
    except:
        aux_cases = xr.open_dataset(cases_date_dir + c + '_f_' + s + '.nc') \
            .rename({'sst': 'index'})

    case_sel_dmi = SelectVariables(aux_cases, data_dates_dmi_or)
    case_sel_dmi = case_sel_dmi.sortby(case_sel_dmi.time.dt.month)

    case_sel_dmi_n34 = SelectVariables(aux_cases, data_dates_n34_or)
    case_sel_dmi_n34 = case_sel_dmi_n34.sortby(case_sel_dmi_n34.time.dt.month)
    # 2.1 uniendo var, dmi y n34
    data_merged = xr.Dataset(
        data_vars=dict(
            var=(['time'], BOX['var'].values),
            dmi=(['time'], case_sel_dmi.sst.values),
            n34=(['time'], case_sel_dmi_n34.sst.values),
        ),
        coords=dict(
            time=BOX.time.values
        )
    )

    bins_aux_dmi = bins_by_cases_dmi[c_count]
    bins_aux_n34 = bins_by_cases_n34[c_count]

    # 3. Seleccion en cada bin
    anom_bin_main = list()
    num_bin_main = list()
    for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
        bins_aux = data_merged.where(SelectBins(data_merged.dmi,
                                                bin_limits[bins_aux_dmi[ba_dmi]][0],
                                                bin_limits[bins_aux_dmi[ba_dmi]][1]))
        anom_bin = list()
        num_bin = list()
        for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
            bin_f = bins_aux.where(SelectBins(bins_aux.n34,
                                              bin_limits[bins_aux_n34[ba_n34]][0],
                                              bin_limits[bins_aux_n34[ba_n34]][1]))
            anom_bin.append(bin_f['var'].mean('time'))
            num_bin.append(len(np.where(~np.isnan(bin_f.dmi))[0])) # en algunos casos hay un CAMPO faltante...

        anom_bin_main.append(anom_bin)
        num_bin_main.append(num_bin)

    return anom_bin_main, num_bin_main

def Plot2D(aux, aux_num, cmap, vmin, vmax, levels, dpi, title, name_fig, save, color_thr=20):
    fig = plt.figure(dpi=dpi, figsize=(8,7))
    ax = fig.add_subplot(111)
    im = ax.imshow(aux, cmap=cmap, vmin=vmin, vmax=vmax)

    for i in range(0, 16):
        for j in range(0, 16):
            if np.abs(aux[i,j])>color_thr:
                color_num = 'white'
            else:
                color_num = 'k'
            if ~np.isnan(aux_num[i, j]):
                ax.text(j, i, aux_num[i, j].astype(np.int64), ha='center', va='center', color=color_num)

    ax.set_xticks(np.arange(-0.5, 16, 1), np.linspace(-4, 4, 17))
    ax.set_yticks(np.arange(-0.5, 16, 1), np.linspace(-4, 4, 17))
    ax.set_xlim([-.5, 15])
    ax.set_ylim([-.5, 15])

    ax.set_ylabel('Niño3.4 - SST index (of std)', fontsize=11)
    ax.set_xlabel('DMI - SST index (of std)', fontsize=11)
    fig.suptitle(title, size=12)

    plt.axhline(y=8.5, color='k', linestyle='-', linewidth=2)
    plt.axhline(y=6.5, color='k', linestyle='-', linewidth=2)
    plt.axvline(x=8.5, color='k', linestyle='-', linewidth=2)
    plt.axvline(x=6.5, color='k', linestyle='-', linewidth=2)

    ax.set_xticks(np.arange(0, 16, 0.5), minor=True)
    ax.set_yticks(np.arange(0, 16, 0.5), minor=True)
    ax.margins(0)
    ax.grid(which='major', alpha=0.5, color='k')
    plt.colorbar(im, ticks=levels,
                 fraction=0.046, pad=0.04,
                 boundaries=levels)
    plt.tight_layout()
    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

def Plot2D_1_1(aux, aux_num, cmap, vmin, vmax, levels, dpi, title, name_fig, save, color_thr=20):
    fig = plt.figure(dpi=dpi, figsize=(8, 7))
    ax = fig.add_subplot(111)
    im = ax.imshow(aux, cmap=cmap, vmin=vmin, vmax=vmax)

    for i in range(0, 9):
        for j in range(0, 9):
            if np.abs(aux[i, j]) > color_thr:
                color_num = 'white'
            else:
                color_num = 'k'
            if ~np.isnan(aux_num[i, j]):
                ax.text(j, i, aux_num[i, j].astype(np.int64), ha='center', va='center', color=color_num)

    ax.set_xticks(np.arange(-0.5, 9, 1), np.linspace(-4.5, 4.5, 10))
    ax.set_yticks(np.arange(-0.5, 9, 1), np.linspace(-4.5, 4.5, 10))
    ax.set_xlim([-.5, 8.5])
    ax.set_ylim([-.5, 8.5])

    ax.set_ylabel('Niño3.4 - SST index (of std)', fontsize=11)
    ax.set_xlabel('DMI - SST index (of std)', fontsize=11)
    fig.suptitle(title, size=12)

    plt.axhline(y=4.5, color='k', linestyle='-', linewidth=2)
    plt.axhline(y=3.5, color='k', linestyle='-', linewidth=2)
    plt.axvline(x=4.5, color='k', linestyle='-', linewidth=2)
    plt.axvline(x=3.5, color='k', linestyle='-', linewidth=2)

    ax.set_xticks(np.arange(0, 9, 1), minor=True)
    ax.set_yticks(np.arange(0, 9, 1), minor=True)
    ax.margins(0)
    ax.grid(which='major', alpha=0.5, color='k')
    plt.colorbar(im, ticks=levels,
                 fraction=0.046, pad=0.04,
                 boundaries=levels)
    plt.tight_layout()
    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

########################################################################################################################
# Region seleccionada con k-means a paritr del CFSv2 ------------------------------------------------------------------#
pre_box = [0, 1, 2]
pre_box_lat = [[-60, 0], [-60, -20], [-60, -20]]
pre_box_lon = [[290, 330], [290, 330], [280, 290]]
#----------------------------------------------------------------------------------------------------------------------#
cases = ['dmi_puros_pos', 'dmi_puros_neg',
        'n34_puros_pos', 'n34_puros_neg',
        'sim_pos', 'sim_neg',
        'dmi_neg_n34_pos', 'dmi_pos_n34_neg',
        'neutros']

# Limitesde cada "bin"
bin_limits = [[-4, -3.5], [-3.5, -3], [-3, -2.5], [-2.5, -2], [-2, -1.5], [-1.5, -1], [-1, -.5], # 6
              [-.5, 0], [0, .5], # 8
              [.5, 1], [1, 1.5], [1.5, 2], [2, 2.5], [2.5, 3], [3, 3.5], [3.5, 4]] # 15

# limites que se van a usar en cada caso.
# en un determinado rango de valores. ej. los dmi_puros_positivos, son dmi > 0.5 con -0.5>n34<0.5

bins_by_cases_dmi = [[9, 10, 11, 12, 13, 14, 15], [0, 1, 2, 3, 4, 5, 6],
                     [7, 8], [7, 8],
                     [9, 10, 11, 12, 13, 14, 15], [0, 1, 2, 3, 4, 5, 6],
                     [0, 1, 2, 3, 4, 5, 6],[9, 10, 11, 12, 13, 14, 15],
                     [7, 8]]

bins_by_cases_n34 = [[7, 8], [7, 8],
                     [9, 10, 11, 12, 13, 14, 15], [0, 1, 2, 3, 4, 5, 6],
                     [9, 10, 11, 12, 13, 14, 15], [0, 1, 2, 3, 4, 5, 6],
                     [9, 10, 11, 12, 13, 14, 15], [0, 1, 2, 3, 4, 5, 6],
                     [7, 8]]

# Limitesde cada "bin"
bin_limits_1 = [[-4.5, -3.5], [-3.5, -2.5], [-2.5,-1.5],
              [-1.5, -0.5], [-0.5, 0.5], [0.5, 1.5],
              [1.5, 2.5], [2.5, 3.5], [3.5, 4.5]]

# limites que se van a usar en cada caso.
# en un determinado rango de valores. ej. los dmi_puros_positivos, son dmi > 0.5 con -0.5>n34<0.5

bins_by_cases_dmi_1 = [[5,6,7,8], [0,1,2,3],
                     [4], [4],
                     [5,6,7,8], [0, 1, 2, 3],
                     [0, 1, 2, 3],[5,6,7,8],
                     [4]]

bins_by_cases_n34_1 = [[4], [4],
                     [5,6,7,8], [0,1,2,3],
                     [5,6,7,8], [0, 1, 2, 3],
                     [5,6,7,8], [0, 1, 2, 3],
                     [4]]

cbar_Prueba = colors.ListedColormap(['#002A3D','#074D4F', '#1E6D5A' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#DCBC75', '#995D13','#6A3D07','#543005','#3F2404'][::-1])
cbar_Prueba.set_under('#3F2404')
cbar_Prueba.set_over('#002A3D')
cbar_Prueba.set_bad(color='white')

cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')
########################################################################################################################

# fix_factor = [30, 1]
# variables = ['prec', 'tref']
# variables = ['prec']
# color_thr = [20, 0.8]
# vmin = [-44, -1]
# vmax = [44, 1]
# colorbars = [cbar_Prueba, cbar_t]
# levels = [np.arange(-44, 44, 8), np.linspace(-1,1,12)]

fix_factor = 30
variables = ['prec']
color_thr = 20
vmin = -44
vmax = 44
colorbars = cbar_Prueba
levels = np.arange(-44, 44, 8)
seasons = ['JJA', 'SON']
mmonth = [7,10]

seasons = ['SON']
mmonth = [10]

from ENSO_IOD_Funciones import MakeMask

# cajas ------------------------------------------------------------------------
box_name = ['S-SESA', 'N-SESA', 'NeB', 'Chile-Cuyo']# 'Patagonia']
box_lats = [[-39,-25], [-29,-17], [-15,2], [-40,-30]]
box_lons = [[296, 306], [305, 315], [311,325], [285,293]]#, [288,300]]
#------------------------------------------------------------------------------#
v_count = 0
for v in variables:
    for bl, bt, name in zip(box_lons, box_lats, box_name):
        for s, mm in zip(seasons, mmonth):
            neutro = xr.open_dataset(
                '/pikachu/datos/luciano.andrian/cases_fields/'
                + v + '_neutros_' + s + '.nc').rename(
                {v: 'var'})
            neutro *= fix_factor
            mask = MakeMask(neutro, 'var')

            aux = np.zeros((9, 9))
            aux_num = np.zeros((9, 9))
            # aux = np.zeros((16, 16))
            # aux_num = np.zeros((16, 16))

            for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
                cases_bin, num_bin = BinsByCases(
                    v=v, v_name=v, fix_factor=fix_factor, s=s, mm=mm,
                    c=cases[c_count], c_count=c_count,
                    box_lat=bt, box_lon=bl,
                    bin_limits=bin_limits_1,
                    bins_by_cases_dmi=bins_by_cases_dmi_1,
                    bins_by_cases_n34=bins_by_cases_n34_1,
                    mask=True)

                bins_aux_dmi = bins_by_cases_dmi_1[c_count]
                bins_aux_n34 = bins_by_cases_n34_1[c_count]
                anom_bin_main = list()
                # loops en las bins para el dmi segun case
                for ba_dmi in range(0, len(bins_aux_dmi)):
                    # loop en las correspondientes al n34 segun case
                    for ba_n34 in range(0, len(bins_aux_n34)):
                        aux[bins_aux_n34[ba_n34],
                            bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                        aux_num[bins_aux_n34[ba_n34],
                                bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

            aux_num[np.where(aux_num == 0)] = np.nan
            Plot2D_1_1(aux, aux_num, dpi=dpi, color_thr=color_thr,
                       cmap=colorbars, vmin=vmin, vmax=vmax, levels=levels,
                       title='Precipitación - ' + name + ' - ' + s,
                       name_fig=v + '_box_' + name + '_cuad_' + s,
                       save=save)

################################################################################
