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
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/' # índices por estaciones
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/' # campos de las variables PP  ( y T cuando ande)
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/bins/2D/'
save = True
dpi = 400
# Funciones ############################################################################################################
def SelectBins(data, min, max):
    return (data >= min) & (data < max)

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
                mask=False, data_mask=None):

    # 1. se abren los archivos de los índices (completos y se pesan por su SD)
    # tambien los archivos de la variable en cuestion pero para cada "case" = c

    data_dates_dmi_or = xr.open_dataset(dates_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc')
    data_dates_dmi_or /= data_dates_dmi_or.mean('r').std() #corregido!

    data_dates_n34_or = xr.open_dataset(dates_dir + 'N34_' + s + '_Leads_r_CFSv2.nc')
    data_dates_n34_or /= data_dates_n34_or.mean('r').std()

    # 1.1 Climatología y case
    clim = xr.open_dataset(cases_dir + v + '_' + s.lower() + '.nc').rename({'prec': 'var'}) * fix_factor
    case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s.upper() + '.nc').rename({v_name: 'var'}) * fix_factor
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
        # PROBANDO!!!
        BOX = anom*data_mask
        BOX = BOX.sel(lat=slice(min(box_lat), max(box_lat)), lon=slice(min(box_lon), max(box_lon)))
        BOX = BOX.mean(['lon', 'lat'])
        BOX = BOX.sortby(BOX.time.dt.month)
    else:
        BOX = anom.sel(lat=slice(min(box_lat), max(box_lat)), lon=slice(min(box_lon), max(box_lon)))
        BOX = BOX.mean(['lon', 'lat'])
        BOX = BOX.sortby(BOX.time.dt.month)

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
            num_bin.append(len(np.where(~np.isnan(bin_f['var']))[0]))

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

# mascara para el oceano ----------------------------------------------------------------------------------------------#
aux = xr.open_dataset('/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/pp_gpcc_d_w_c_1920-2020_1.nc')
mask = aux.sel(lon=slice(270, 330), lat=slice(15, -60))
mask = aux.interp(lon=np.linspace(275, 330, 56), lat=np.linspace(-60, 15, 76))
mask = mask.mean('time')
mask = xr.where(np.isnan(mask), mask, 1)

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


color = ['#E4FDA9', '#F9FD7C', plt.cm.tab20c(7), plt.cm.tab20c(6), plt.cm.tab20c(5), plt.cm.tab20c(4),
             plt.cm.tab20b(15), plt.cm.tab20b(14), plt.cm.tab20b(13), plt.cm.tab20b(12)]
cmap_pos = colors.ListedColormap(color)
########################################################################################################################
mm = 7
for s in ['JJA', 'JAS', 'ASO', 'SON']:
    aux = np.zeros((16,16))
    aux_num = np.zeros((16, 16))
    for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        cases_bin, num_bin = BinsByCases(v='prec', v_name='prec', fix_factor=30,
                                         s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                         box_lat=[-35, -29], box_lon=[300, 310],
                                         bin_limits=bin_limits, bins_by_cases_dmi=bins_by_cases_dmi,
                                         bins_by_cases_n34=bins_by_cases_n34)

        bins_aux_dmi = bins_by_cases_dmi[c_count]
        bins_aux_n34 = bins_by_cases_n34[c_count]
        anom_bin_main = list()
        for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
            for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

    aux_num[np.where(aux_num == 0)] = np.nan
    Plot2D(aux, aux_num, dpi=dpi,
           cmap=cbar_Prueba, vmin=-44, vmax=44, levels=np.arange(-44, 44, 8),
           title='S-SESA - ' + s + ' PP',
           name_fig='S_SESA_cuad_' + s,
           save=save)

    Plot2D(aux_num, aux_num, dpi=300,
           cmap=cmap_pos, vmin=0, vmax=150, levels=np.linspace(0, 150, 11),
           title='Samples - S-SESA - ' + s + ' PP',
           name_fig='num_cuad_' + s + '.jpg',
           save=save, color_thr=120)
    mm += 1

# 1 _ 1 ################################################################################################################
# mm = 7
# for s in ['JJA', 'JAS', 'ASO', 'SON']:
#     aux = np.zeros((9,9))
#     aux_num = np.zeros((9, 9))
#     for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
#         cases_bin, num_bin = BinsByCases(v='prec', v_name='prec', fix_factor=30,
#                                          s=s, mm=mm, c=cases[c_count], c_count=c_count,
#                                          box_lat=[-35, -29], box_lon=[300, 310],
#                                          bin_limits=bin_limits_1, bins_by_cases_dmi=bins_by_cases_dmi_1,
#                                          bins_by_cases_n34=bins_by_cases_n34_1)
#
#         bins_aux_dmi = bins_by_cases_dmi_1[c_count]
#         bins_aux_n34 = bins_by_cases_n34_1[c_count]
#         anom_bin_main = list()
#         for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
#             for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
#                 aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
#                 aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]
#
#     aux_num[np.where(aux_num == 0)] = np.nan
#     Plot2D_1_1(aux, aux_num, dpi=dpi,
#            cmap=cbar_Prueba, vmin=-44, vmax=44, levels=np.arange(-44, 44, 8),
#            title='S-SESA - ' + s + ' PP',
#            name_fig='S_SESA_cuad_' + s,
#            save=save)
#     mm += 1
########################################################################################################################
mm = 7
for s in ['JJA', 'JAS', 'ASO', 'SON']:
    aux = np.zeros((16,16))
    aux_num = np.zeros((16, 16))
    for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        cases_bin, num_bin = BinsByCases(v='prec', v_name='prec', fix_factor=30,
                                         s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                         box_lat=[-29, -20], box_lon=[300, 320])

        bins_aux_dmi = bins_by_cases_dmi[c_count]
        bins_aux_n34 = bins_by_cases_n34[c_count]
        anom_bin_main = list()
        for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
            for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

    aux_num[np.where(aux_num == 0)] = np.nan
    Plot2D(aux, aux_num, dpi=dpi,
           cmap=cbar_Prueba, vmin=-44, vmax=44, levels=np.arange(-44, 44, 8),
           title='N-SESA - ' + s + ' PP',
           name_fig='N_SESA_cuad_' + s,
           save=save)
    mm +=1




# Esto sale de el código de kmeans
# Hay mas de un lugar con el mismo cluster
plt.imshow(mask_cluster.cluster);plt.show()
# voy a usar la mascara y seleccionar un cuadrado cerca de SESA
plt.imshow(mask_cluster.sel(lat=slice(-40,-20), lon=slice(290,320)).cluster);plt.show()
mask_cluster = xr.where(mask_cluster==7, 1, mask_cluster)
mm = 7
for s in ['JJA', 'JAS', 'ASO', 'SON']:
    aux = np.zeros((16,16))
    aux_num = np.zeros((16, 16))
    for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        cases_bin, num_bin = BinsByCases(v='prec', v_name='prec', fix_factor=30,
                                         s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                         box_lat=[-40, -20], box_lon=[290, 320], # cuadrado al rededor de la mascara
                                         bin_limits=bin_limits, bins_by_cases_dmi=bins_by_cases_dmi,
                                         bins_by_cases_n34=bins_by_cases_n34,
                                         mask=True, data_mask=mask_cluster.rename({'cluster':'var'}))

        bins_aux_dmi = bins_by_cases_dmi[c_count]
        bins_aux_n34 = bins_by_cases_n34[c_count]
        anom_bin_main = list()
        for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
            for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

    aux_num[np.where(aux_num == 0)] = np.nan
    Plot2D(aux, aux_num, dpi=dpi,
           cmap=cbar_Prueba, vmin=-44, vmax=44, levels=np.arange(-44, 44, 8),
           title='S-SESA - ' + s + ' PP',
           name_fig='kmean_PRUEBA_S_SESA_cuad_' + s,
           save=save)
    mm += 1


# 1 _ 1 ################################################################################################################
mm = 7
for s in ['JJA', 'JAS', 'ASO', 'SON']:
    aux = np.zeros((9,9))
    aux_num = np.zeros((9, 9))
    for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
        cases_bin, num_bin = BinsByCases(v='prec', v_name='prec', fix_factor=30,
                                         s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                         box_lat=[-40, -20], box_lon=[290, 320], # cuadrado al rededor de la mascara
                                         bin_limits=bin_limits, bins_by_cases_dmi=bins_by_cases_dmi,
                                         bins_by_cases_n34=bins_by_cases_n34,
                                         mask=True, data_mask=mask_cluster.rename({'cluster':'var'}))
        bins_aux_dmi = bins_by_cases_dmi_1[c_count]
        bins_aux_n34 = bins_by_cases_n34_1[c_count]
        anom_bin_main = list()
        for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
            for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

    aux_num[np.where(aux_num == 0)] = np.nan
    Plot2D_1_1(aux, aux_num, dpi=dpi,
           cmap=cbar_Prueba, vmin=-44, vmax=44, levels=np.arange(-44, 44, 8),
           title='S-SESA - ' + s + ' PP',
           name_fig='kmean_PRUEBA__S_SESA_cuad_1-1' + s,
           save=save)
    mm += 1
