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
save = False
dpi = 50
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
                mask=False, data_mask=None):

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
    anom = anom.sel(lat=slice(min(box_lat), max(box_lat)), lon=slice(min(box_lon), max(box_lon)))
    if mask:
        # PROBANDO!!!
        BOX = anom*data_mask
        #BOX = BOX.sel(lat=slice(min(box_lat), max(box_lat)), lon=slice(min(box_lon), max(box_lon)))
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

fix_factor = [30, 1]
variables = ['prec', 'tref']
color_thr = [20, 0.8]
vmin = [-44, -1]
vmax = [44, 1]
colorbars = [cbar_Prueba, cbar_t]
levels = [np.arange(-44, 44, 8), np.linspace(-1,1,12)]

v_count = 0
for v in variables:

    # SESA*
    data_cfsv2 = CFSv2Data(v)
    xr_cluster_label = KmMaks(data_cfsv2, v=v, l=0, pb=0)
    mask_cluster = xr.where(xr_cluster_label == 3, 1, np.nan).sel(lat=slice(-40, -20), lon=slice(290, 320))
    # plt.imshow(mask_cluster.cluster);
    # plt.show()

    mm = 7
    for s in ['JJA', 'JAS', 'ASO', 'SON']:
        aux = np.zeros((16, 16))
        aux_num = np.zeros((16, 16))
        for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
            cases_bin, num_bin = BinsByCases(v=v, v_name=v, fix_factor=fix_factor[v_count],
                                             s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                             box_lat=[-40, -20], box_lon=[290, 320],
                                             # cuadrado al rededor de la mascara
                                             bin_limits=bin_limits, bins_by_cases_dmi=bins_by_cases_dmi,
                                             bins_by_cases_n34=bins_by_cases_n34,
                                             mask=True, data_mask=mask_cluster.rename({'cluster': 'var'}))

            bins_aux_dmi = bins_by_cases_dmi[c_count]
            bins_aux_n34 = bins_by_cases_n34[c_count]
            anom_bin_main = list()
            for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
                for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                    aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                    aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

        aux_num[np.where(aux_num == 0)] = np.nan
        Plot2D(aux, aux_num, dpi=dpi, color_thr=color_thr[v_count],
               cmap=colorbars[v_count], vmin=vmin[v_count], vmax=vmax[v_count], levels=levels[v_count],
               title='SESA* - ' + s + ' ' + v,
               name_fig=v + '_kmean_SESA_cuad_' + s,
               save=save)
        mm += 1

    # Out SESA?
    mask_cluster = xr.where((xr_cluster_label == 2) | (xr_cluster_label == 0), 1, np.nan) \
        .sel(lat=slice(-40, -10), lon=slice(290, 320))
    # plt.imshow(mask_cluster.cluster);
    # plt.show()

    mm = 7
    for s in ['JJA', 'JAS', 'ASO', 'SON']:
        aux = np.zeros((16, 16))
        aux_num = np.zeros((16, 16))
        for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
            cases_bin, num_bin = BinsByCases(v=v, v_name=v, fix_factor=fix_factor[v_count],
                                             s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                             box_lat=[-40, -10], box_lon=[290, 320],
                                             # cuadrado al rededor de la mascara
                                             bin_limits=bin_limits, bins_by_cases_dmi=bins_by_cases_dmi,
                                             bins_by_cases_n34=bins_by_cases_n34,
                                             mask=True, data_mask=mask_cluster.rename({'cluster': 'var'}))

            bins_aux_dmi = bins_by_cases_dmi[c_count]
            bins_aux_n34 = bins_by_cases_n34[c_count]
            anom_bin_main = list()
            for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
                for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                    aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                    aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

        aux_num[np.where(aux_num == 0)] = np.nan
        Plot2D(aux, aux_num, dpi=dpi, color_thr=color_thr[v_count],
               cmap=colorbars[v_count], vmin=vmin[v_count], vmax=vmax[v_count], levels=levels[v_count],
               title='Out SESA? - ' + s + ' PP',
               name_fig=v + '_kmean_out_SESA_cuad_' + s,
               save=save)
        mm += 1

    # 1 _ 1 ############################################################################################################

    # SESA*
    mask_cluster = xr.where(xr_cluster_label == 3, 1, np.nan).sel(lat=slice(-40, -20), lon=slice(290, 320))
    # plt.imshow(mask_cluster.cluster);
    # plt.show()
    mm = 7
    for s in ['JJA', 'JAS', 'ASO', 'SON']:
        aux = np.zeros((9, 9))
        aux_num = np.zeros((9, 9))
        for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
            cases_bin, num_bin = BinsByCases(v=v, v_name=v, fix_factor=fix_factor[v_count],
                                             s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                             box_lat=[-40, -20], box_lon=[290, 320],
                                             # cuadrado al rededor de la mascara
                                             bin_limits=bin_limits_1, bins_by_cases_dmi=bins_by_cases_dmi_1,
                                             bins_by_cases_n34=bins_by_cases_n34_1,
                                             mask=True, data_mask=mask_cluster.rename({'cluster': 'var'}))
            bins_aux_dmi = bins_by_cases_dmi_1[c_count]
            bins_aux_n34 = bins_by_cases_n34_1[c_count]
            anom_bin_main = list()
            for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
                for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                    aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                    aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

        aux_num[np.where(aux_num == 0)] = np.nan
        Plot2D_1_1(aux, aux_num, dpi=dpi,
                   cmap=colorbars[v_count], vmin=vmin[v_count], vmax=vmax[v_count], levels=levels[v_count],
                   title='SESA* - ' + s + ' ' + v,
                   name_fig=v + '_kmean_SESA_cuad_1-1' + s,
                   save=save)
        mm += 1

    # Out SESA?
    mask_cluster = xr.where((xr_cluster_label == 2) | (xr_cluster_label == 0), 1, np.nan) \
        .sel(lat=slice(-40, -10), lon=slice(290, 320))
    # plt.imshow(mask_cluster.cluster);
    # plt.show()

    mm = 7
    for s in ['JJA', 'JAS', 'ASO', 'SON']:
        aux = np.zeros((9, 9))
        aux_num = np.zeros((9, 9))
        for c_count in [0, 1, 2, 3, 4, 5, 6, 7, 8]:
            cases_bin, num_bin = BinsByCases(v=v, v_name=v, fix_factor=fix_factor[v_count],
                                             s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                             box_lat=[-40, -10], box_lon=[290, 320],
                                             # cuadrado al rededor de la mascara
                                             bin_limits=bin_limits_1, bins_by_cases_dmi=bins_by_cases_dmi_1,
                                             bins_by_cases_n34=bins_by_cases_n34_1,
                                             mask=True, data_mask=mask_cluster.rename({'cluster': 'var'}))
            bins_aux_dmi = bins_by_cases_dmi_1[c_count]
            bins_aux_n34 = bins_by_cases_n34_1[c_count]
            anom_bin_main = list()
            for ba_dmi in range(0, len(bins_aux_dmi)):  # loops en las bins para el dmi segun case
                for ba_n34 in range(0, len(bins_aux_n34)):  # loop en las correspondientes al n34 segun case
                    aux[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = cases_bin[ba_dmi][ba_n34]
                    aux_num[bins_aux_n34[ba_n34], bins_aux_dmi[ba_dmi]] = num_bin[ba_dmi][ba_n34]

        aux_num[np.where(aux_num == 0)] = np.nan
        Plot2D_1_1(aux, aux_num, dpi=dpi,
                   cmap=colorbars[v_count], vmin=vmin[v_count], vmax=vmax[v_count], levels=levels[v_count],
                   title='Out SESA? - ' + s + ' ' + v,
                   name_fig=v + '_kmean_out_SESA_cuad_1-1' + s,
                   save=save)
        mm += 1
    v_count += 1
########################################################################################################################
