"""
Prueba Validación CFSv2
"""
################################################################################
save = False
################################################################################
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import os
#from ENSO_IOD_Funciones import OpenDatasets
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import matplotlib.path as mpath
import warnings
warnings.filterwarnings("ignore")
from scipy.stats import ttest_ind
from ENSO_IOD_Funciones import MakeMask
################################################################################
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
              'nc_composites_dates_no_ind_sst_anom/' #fechas
# obs
data_dir_t_pp = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'data_obs_d_w_c/' #T y PP ya procesados

# CFSv2 # fechas ya seleccionadas
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/'

out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Validation/'

#Plot
if save:
    dpi = 300
else:
    dpi = 100

# Functions ####################################################################
def CompositeSimple(original_data, index, mmin, mmax, mean=True):
    def is_months(month, mmin, mmax):
        return (month >= mmin) & (month <= mmax)

    if len(index) != 0:
        comp_field = original_data.sel(
            time=original_data.time.dt.year.isin([index]))
        comp_field = comp_field.sel(
            time=is_months(month=comp_field['time.month'], mmin=mmin, mmax=mmax))

        if len(comp_field.time) != 0:
            if mean:
                comp_field = comp_field.mean(['time'], skipna=True)
            else:
                pass
        else:  # si sólo hay un año
            comp_field = comp_field.drop_dims(['time'])

        return comp_field
    else:
        print(' len index = 0')

def CaseComp(data, s, mmonth, c, two_variables=False, data2=None):
    """
    Las fechas se toman del periodo 1920-2020 basados en el DMI y N34 con ERSSTv5
    Cuando se toman los periodos 1920-1949 y 1950_2020 las fechas que no pertencen
    se excluyen de los composites en CompositeSimple()
    """
    mmin = mmonth[0]
    mmax = mmonth[-1]

    aux = xr.open_dataset(nc_date_dir + '1920_2020' + '_' + s + '.nc')
    neutro = aux.Neutral

    try:
        case = aux[c]
        case = case.where(case >= 1982)
        neutro = neutro.where(neutro>=1982)
        aux.close()

        case_num = len(case.values[np.where(~np.isnan(case.values))])

        neutro_comp = CompositeSimple(original_data=data, index=neutro,
                                      mmin=mmin, mmax=mmax, mean=True)
        data_comp = CompositeSimple(original_data=data, index=case,
                                    mmin=mmin, mmax=mmax, mean=False)

        comp = data_comp - neutro_comp

        if two_variables:
            neutro_comp2 = CompositeSimple(original_data=data2, index=neutro,
                                           mmin=mmin, mmax=mmax)
            data_comp2 = CompositeSimple(original_data=data2, index=case,
                                         mmin=mmin, mmax=mmax)

            comp2 = data_comp2 - neutro_comp2
        else:
            comp2 = None
    except:
        print('Error en ' + s + c)

    if two_variables:
        return comp, case_num, comp2
    else:
        return comp, case_num, data_comp

def Plot(comp, levels, cmap, step1, contour1=True,
         two_variables=False, comp2=None, levels2=np.linspace(-1,1,13), step2=4,
         mapa='sa', title='title', name_fig='name_fig', dpi=100, save=save,
         comp_sig=None, color_sig='k', significance=True, linewidht2=.5,
         color_map='#4B4B4B', out_dir=out_dir, proj='eq'):

    import matplotlib.pyplot as plt

    if mapa.lower()=='sa':
        fig_size = (5, 6)
        extent= [270, 330, -60, 20]
        xticks = np.arange(270, 330, 10)
        yticks = np.arange(-60, 40, 20)

    elif mapa.lower()=='tropical':
        fig_size = (7, 2)
        extent = [50, 270, -20, 20]
        xticks = np.arange(50, 270, 60)
        yticks = np.arange(-20, 20, 20)

    elif mapa.lower()=='hs':
        fig_size = (8, 3)
        extent = [0, 359, -90, 20]
        xticks = np.arange(0, 330, 30)
        yticks = np.arange(-90, 20, 10)
        if proj != 'eq':
            fig_size = (5, 5)
    else:
        fig_size = (8, 3)
        extent = [30, 330, -80, 20]
        xticks = np.arange(30, 330, 30)
        yticks = np.arange(-80, 20, 10)
        if proj != 'eq':
            fig_size = (5, 5)

    levels_contour = levels.copy()
    comp_var = comp['var']
    if isinstance(levels, np.ndarray):
        levels_contour = levels[levels != 0]
    else:
        levels_contour.remove(0)

    if two_variables:
        levels_contour2 = levels2.copy()
        comp_var2=comp2['var']
        if isinstance(levels2, np.ndarray):
            levels_contour2 = levels2[levels2 != 0]
        else:
            levels_contour2.remove(0)

    crs_latlon = ccrs.PlateCarree()
    fig = plt.figure(figsize=fig_size, dpi=dpi)
    if proj=='eq':
        ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
        ax.set_extent(extent, crs=crs_latlon)
    else:
        ax = plt.axes(projection=ccrs.SouthPolarStereo(central_longitude=200))
        ax.set_extent([30, 340, -90, 0],
                          ccrs.PlateCarree(central_longitude=200))



    if two_variables:
        ax.contour(comp2.lon[::step2], comp2.lat[::step2],
                   comp_var2[::step2, ::step2],
                   linewidths=linewidht2, levels=levels_contour2,
                   transform=crs_latlon, colors='k')
    else:
        if contour1:
            ax.contour(comp.lon[::step1], comp.lat[::step1],
                       comp_var[::step1, ::step1],
                       linewidths=.8, levels=levels_contour,
                       transform=crs_latlon, colors='black')

    im = ax.contourf(comp.lon[::step1], comp.lat[::step1],
                     comp_var[::step1, ::step1],
                     levels=levels, transform=crs_latlon, cmap=cmap,
                     extend='both')

    if significance:
        colors_l = [color_sig, color_sig]
        comp_sig_var = comp_sig['var']
        cs = ax.contourf(comp_sig.lon, comp_sig.lat, comp_sig_var,
                         transform=crs_latlon, colors='none',
                         hatches=["..", ".."], extend='lower')
        for i, collection in enumerate(cs.collections):
            collection.set_edgecolor(colors_l[i % len(colors_l)])

        for collection in cs.collections:
            collection.set_linewidth(0.)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                   edgecolor=color_map)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)

    if proj=='eq':
        ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
        ax.set_xticks(xticks, crs=crs_latlon)
        ax.set_yticks(yticks, crs=crs_latlon)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        ax.xaxis.set_major_formatter(lon_formatter)
        ax.yaxis.set_major_formatter(lat_formatter)
    else:
        gls = ax.gridlines(draw_labels=True, crs=crs_latlon, lw=0.3,
                           color="gray", y_inline=True,
                           xlocs=range(-180, 180, 30),
                           ylocs=np.arange(-80, 0, 20))

        r_extent = 1.2e7
        ax.set_xlim(-r_extent, r_extent)
        ax.set_ylim(-r_extent, r_extent)
        circle_path = mpath.Path.unit_circle()
        circle_path = mpath.Path(circle_path.vertices.copy() * r_extent,
                                 circle_path.codes.copy())
        ax.set_boundary(circle_path)
        ax.set_frame_on(False)
        plt.draw()
        for ea in gls._labels:
            pos = ea[2].get_position()
            if (pos[0] == 150):
                ea[2].set_position([0, pos[1]])

    ax.tick_params(labelsize=7)
    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

def Weights(data):
    weights = np.transpose(np.tile(np.cos(data.lat * np.pi / 180),
                                   (len(data.lon), 1)))
    data_w = data * weights
    return data_w
################################################################################
seasons = ('JJA', 'SON')
min_max_months = [[6,8],[9,11]]

# cases = ['DMI_sim_pos', 'DMI_sim_neg',
#          'DMI_un_pos', 'DMI_un_neg',
#          'N34_un_pos', 'N34_un_neg']
#
# cases_cfs = ['sim_pos', 'sim_neg',
#              'dmi_puros_pos', 'dmi_puros_neg',
#              'n34_puros_pos', 'n34_puros_neg']
#
# title_case = ['DMI-ENSO simultaneous positive phase ',
#               'DMI-ENSO simultaneous negative phase ',
#               'DMI pure positive phase ',
#               'DMI pure negative phase ',
#               'ENSO pure positive phase ',
#               'ENSO pure negative phase ']

cases =  ['DMI_pos', 'DMI_neg', 'N34_pos', 'N34_neg']
cases_cfs = ['dmi_pos', 'dmi_neg', 'n34_pos', 'n34_neg']

# cases = ['DMI_un_pos', 'DMI_un_neg', 'N34_un_pos', 'N34_un_neg']
# cases_cfs = ['dmi_puros_pos', 'dmi_puros_neg','n34_puros_pos', 'n34_puros_neg']

v_name_fig=['pp_gpcc', 'pp_prec']
#------------------------------------------------------------------------------#
# colorbars
cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC',
                                 '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07',
                                 '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89',
                                '#FFCECC', 'white', '#B3DBFF', '#83B9EB',
                                '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')

cbar_t_r = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55',
                                  '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                                  '#2064AF', '#014A9B'])
cbar_t_r.set_under('#691800')
cbar_t_r.set_over('#013774')
cbar_t_r.set_bad(color='white')

cbar_sst = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89',
                                  '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3',
                                  '#2064AF'][::-1])
cbar_sst.set_over('#9B1C00')
cbar_sst.set_under('#014A9B')
cbar_sst.set_bad(color='white')
################################################################################
# HGT--------------------------------------------------------------------------#
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/' \
                'ERA5/1940_2020/'
scale = [-400,-250,-200,-150,-100,-50,0,50,100,150,200,250,400]
v_name = ['HGT200hPa - ERA5']

for v, v_count in zip(['HGT200_SON_mer_d_w'], [0]): # agregar viento
    data = xr.open_dataset(data_dir_era5 + v + '.nc')

    for s, s_count in zip(['SON'], [1]):
        # # neutro CFSv2 --------------------------------------------------------
        # neutro = xr.open_dataset(
        #     cases_dir + 'hgt_neutros_' + s + '.nc').rename({'hgt': 'var'})
        # neutro = neutro.__mul__(9.80665)

        # ENSO_CFSv2 --------------------------------------------------------#
        en_cfs = Weights(xr.open_dataset(
            cases_dir + 'hgt_' + cases_cfs[2] + '_' + s + '.nc')\
            .rename({'hgt': 'var'}).__mul__(9.8))
        ln_cfs = Weights(xr.open_dataset(
            cases_dir + 'hgt_' + cases_cfs[3] + '_' + s + '.nc')\
            .rename({'hgt': 'var'}).__mul__(9.8))

        enso_cfs_test = ttest_ind(
            en_cfs['var'].values, ln_cfs['var'].values, equal_var=False)[1]

        enso_dif_cfs = en_cfs.mean('time') - ln_cfs.mean('time')

        # IOD_CFSv2 -----------------------------------------------------------
        iodp_cfs = Weights(xr.open_dataset(
            cases_dir + 'hgt_' + cases_cfs[0] + '_' + s + '.nc')\
            .rename({'hgt': 'var'}).__mul__(9.8))
        iodn_cfs = Weights(xr.open_dataset(
            cases_dir + 'hgt_' + cases_cfs[1] + '_' + s + '.nc')\
            .rename({'hgt': 'var'}).__mul__(9.8))

        iod_cfs_test = ttest_ind(
            iodp_cfs['var'].values, iodn_cfs['var'].values, equal_var=False)[1]
        iod_dif_cfs = iodp_cfs.mean('time') - iodn_cfs.mean('time')

        # ENSO_obs ------------------------------------------------------------#
        x, y, en_obs = CaseComp(
            data, s, mmonth=min_max_months[s_count], c=cases[2],
            two_variables=False)
        en_obs = en_obs.interp(lon=enso_dif_cfs.lon.values,
                               lat=enso_dif_cfs.lat.values)

        x, y, ln_obs = CaseComp(
            data, s, mmonth=min_max_months[s_count], c=cases[3],
            two_variables=False)
        ln_obs = ln_obs.interp(lon=enso_dif_cfs.lon.values,
                               lat=enso_dif_cfs.lat.values)
        enso_obs_test = ttest_ind(
            en_obs['var'].values, ln_obs['var'].values, equal_var=False)[1]

        enso_dif_obs = en_obs.mean('time') - ln_obs.mean('time')

        # IOD_obs -------------------------------------------------------------#
        x, y, iodp_obs = CaseComp(data, s, mmonth=min_max_months[s_count],
                                  c=cases[0], two_variables=False)
        iodp_obs = iodp_obs.interp(lon=enso_dif_cfs.lon.values,
                                   lat=enso_dif_cfs.lat.values)

        x, y, iodn_obs = CaseComp(data, s, mmonth=min_max_months[s_count],
                                  c=cases[1], two_variables=False)
        iodn_obs = iodn_obs.interp(lon=enso_dif_cfs.lon.values,
                                   lat=enso_dif_cfs.lat.values)

        iod_obs_test = ttest_ind(
            iodp_obs['var'].values, iodn_obs['var'].values, equal_var=False)[1]

        iod_dif_obs = iodp_obs.mean('time') - iodn_obs.mean('time')

        fig = plt.figure(figsize=(8, 3), dpi=dpi)
        gs = fig.add_gridspec(2, 2)

        title = ['EN-LN_ERA5', 'EN-LN_CFSv2',
                 'IOD_pos-neg_ERA5', 'IOD_pos-neg_CFSv2']
        comps = [enso_dif_obs, enso_dif_cfs, iod_dif_obs, iod_dif_cfs]
        tests = [enso_obs_test, enso_cfs_test, iod_obs_test, iod_cfs_test ]
        for field, f_count, f_count_c, t, test in \
                zip(comps,[0,1,0,1],[0,0,1,1], title, tests):

            ax = fig.add_subplot(gs[f_count, f_count_c],
                             projection=ccrs.PlateCarree(central_longitude=180))
            crs_latlon = ccrs.PlateCarree()
            ax.set_extent([30, 340, -80, 10], crs=crs_latlon)

            ax.contour(field.lon, field.lat, field['var'],
                       transform=crs_latlon, colors='k',
                       levels=scale, linewidths=1)

            im = ax.contourf(field.lon, field.lat,
                             field.where(test<0.05)['var'],
                             levels=scale,
                             transform=crs_latlon, cmap=cbar_t, extend='both')

            cb = plt.colorbar(im, fraction=0.017, pad=0.025, shrink=0.8)
            cb.ax.tick_params(labelsize=5)
            ax.add_feature(cartopy.feature.LAND, facecolor='white',
                           edgecolor='#4B4B4B')
            ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
            ax.coastlines(color='#4B4B4B', linestyle='-', alpha=1)
            ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
            ax.set_xticks(np.arange(30, 340, 25), crs=crs_latlon)
            ax.set_yticks(np.arange(-80, 20, 10), crs=crs_latlon)
            lon_formatter = LongitudeFormatter(zero_direction_label=True)
            lat_formatter = LatitudeFormatter()
            ax.xaxis.set_major_formatter(lon_formatter)
            ax.yaxis.set_major_formatter(lat_formatter)
            ax.tick_params(labelsize=5)
            ax.set_title(t, fontsize=10)
            plt.tight_layout()

        fig.suptitle('HGT200 - ' + s)
        if save:
            plt.savefig(out_dir + 'HGT_dif_' + s + '.jpg', dpi=dpi)
        else:
            plt.show()

################################################################################
################################################################################
# PP --------------------------------------------------------------------------#
end_nc_file = '.nc'
s = 'SON'
v_name = ['Precipitation - GPCC']
scale = np.linspace(-30, 30, 13)

data = xr.open_dataset(data_dir_t_pp + 'ppgpcc_w_c_d_1_SON.nc')

# neutro CFSv2 --------------------------------------------------------------- #
neutro = (xr.open_dataset(cases_dir + 'prec_neutros_' + s + '.nc')
          .rename({'prec': 'var'}).__mul__(30))

# ENSO_CFSv2 ----------------------------------------------------------------- #
en_cfs = (xr.open_dataset(cases_dir + 'prec_' + cases_cfs[2] + '_' + s + '.nc')
          .rename({'prec': 'var'}).__mul__(30))
ln_cfs = (xr.open_dataset(cases_dir + 'prec_' + cases_cfs[3] + '_' + s + '.nc')
          .rename({'prec': 'var'}).__mul__(30))

enso_cfs_test = ttest_ind(
    en_cfs['var'].values, ln_cfs['var'].values, equal_var=False)[1]

enso_dif_cfs = en_cfs.mean('time') - ln_cfs.mean('time')
enso_dif_cfs = enso_dif_cfs.sel(lat=slice(None, None, -1))

# IOD_CFSv2 ---------------------------------------------------------------------------------------------------#
iodp_cfs = (xr.open_dataset(
    cases_dir + 'prec_' + cases_cfs[0] + '_' + s + '.nc')
            .rename({'prec': 'var'}).__mul__(30))
iodn_cfs = (xr.open_dataset
            (cases_dir + 'prec_' + cases_cfs[1] + '_' + s + '.nc')
            .rename({'prec': 'var'}).__mul__(30))

iod_cfs_test = ttest_ind(
    iodp_cfs['var'].values, iodn_cfs['var'].values, equal_var=False)[1]

iod_dif_cfs = iodp_cfs.mean('time') - iodn_cfs.mean('time')
iod_dif_cfs = iod_dif_cfs.sel(lat=slice(None, None, -1))

# ENSO_obs --------------------------------------------------------------------#
x, y, en_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[2],
                        two_variables=False)
#en_obs = en_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)
x, y, ln_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[3],
                        two_variables=False)
#ln_obs = ln_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)
enso_obs_test = ttest_ind(
    en_obs['var'].values, ln_obs['var'].values, equal_var=False)[1]

enso_dif_obs = en_obs.mean('time') - ln_obs.mean('time')
#
# IOD_obs ---------------------------------------------------------------------#
x, y, iodp_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[0],
                          two_variables=False)
#iodp_obs = iodp_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)
x, y, iodn_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[1],
                          two_variables=False)
#iodn_obs = iodn_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)

iod_obs_test = ttest_ind(
    iodp_obs['var'].values, iodn_obs['var'].values, equal_var=False)[1]

iod_dif_obs = iodp_obs.mean('time') - iodn_obs.mean('time')

fig = plt.figure(figsize=(5, 6), dpi=dpi)
gs = fig.add_gridspec(2, 2)

title = ['EN-LN_' + v_name[0].split()[2], 'EN-LN_CFSv2',
         'IOD_pos-neg_'  + v_name[0].split()[2], 'IOD_pos-neg_CFSv2']
comps = [enso_dif_obs, enso_dif_cfs, iod_dif_obs, iod_dif_cfs]
tests = [enso_obs_test, enso_cfs_test, iod_obs_test, iod_cfs_test]

for field, f_count, f_count_c, t, test in zip(comps,[0,1,0,1],[0,0,1,1],
                                              title, tests):
    field = field*MakeMask(field, 'var')
    ax = fig.add_subplot(gs[f_count, f_count_c],
                         projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270, 330, -60, 20], crs=crs_latlon)

    im = ax.contourf(field.lon, field.lat, field['var'],#.where(test<0.5)['var'],
                     levels=scale,
                     transform=crs_latlon, cmap=cbar_pp, extend='both')

    # aux = field.where(test < 0.05, np.nan)
    # ax.contourf(field.lon, field.lat, aux['var'],
    #             levels=scale, hatches=['***'],
    #             transform=crs_latlon, color='k')

    cb = plt.colorbar(im, fraction=0.040, pad=0.025, shrink=0.8)
    cb.ax.tick_params(labelsize=5)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                   edgecolor='#4B4B4B')
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color='#4B4B4B', linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=5)
    ax.set_title(t, fontsize=10)
    plt.tight_layout()


fig.suptitle('PP - ' + s)
if save:
    plt.savefig(out_dir + 'prec_dif_' + s + '_' + v_name[0].split()[2] +
                '.jpg', dpi=dpi)
else:
    plt.show()



fig = plt.figure(figsize=(5, 6), dpi=dpi)
gs = fig.add_gridspec(2, 2)

title = ['EN-LN_' + v_name[0].split()[2], 'EN-LN_CFSv2',
         'IOD_pos-neg_'  + v_name[0].split()[2], 'IOD_pos-neg_CFSv2']
comps = [enso_dif_obs, enso_dif_cfs, iod_dif_obs, iod_dif_cfs]
tests = [enso_obs_test, enso_cfs_test, iod_obs_test, iod_cfs_test]

for field, f_count, f_count_c, t, test in zip(comps,[0,1,0,1],[0,0,1,1],
                                              title, tests):
    field = field*MakeMask(field, 'var')
    ax = fig.add_subplot(gs[f_count, f_count_c],
                         projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270, 330, -60, 20], crs=crs_latlon)

    im = ax.contourf(field.lon, field.lat, field.where(test<0.1)['var'],
                     levels=scale,
                     transform=crs_latlon, cmap=cbar_pp, extend='both')

    # aux = field.where(test < 0.05, np.nan)
    # ax.contourf(field.lon, field.lat, aux['var'],
    #             levels=scale, hatches=['***'],
    #             transform=crs_latlon, color='k')

    cb = plt.colorbar(im, fraction=0.040, pad=0.025, shrink=0.8)
    cb.ax.tick_params(labelsize=5)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                   edgecolor='#4B4B4B')
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color='#4B4B4B', linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=5)
    ax.set_title(t, fontsize=10)
    plt.tight_layout()


fig.suptitle('PP - ' + s + ' - sig. 0.1')
if save:
    plt.savefig(out_dir + 'prec_dif_test_' + s + '_' + v_name[0].split()[2] +
                '.jpg', dpi=dpi)
else:
    plt.show()
#
# ########################################################################################################################
# T -------------------------------------------------------------------------------------------------------------------#
end_nc_file = '.nc'
s = 'SON'
v_name = ['trefipitation - GPCC']
scale = np.linspace(-1.2, 1.2, 13)

data = xr.open_dataset(data_dir_t_pp + 'tcru_w_c_d_0.25_SON.nc')

# neutro CFSv2 --------------------------------------------------------------- #
neutro = (xr.open_dataset(cases_dir + 'tref_neutros_' + s + '.nc')
          .rename({'tref': 'var'})-273)

# ENSO_CFSv2 ----------------------------------------------------------------- #
en_cfs = (xr.open_dataset(cases_dir + 'tref_' + cases_cfs[2] + '_' + s + '.nc')
          .rename({'tref': 'var'})-273)
ln_cfs = (xr.open_dataset(cases_dir + 'tref_' + cases_cfs[3] + '_' + s + '.nc')
          .rename({'tref': 'var'})-273)

enso_cfs_test = ttest_ind(
    en_cfs['var'].values, ln_cfs['var'].values, equal_var=False)[1]

enso_dif_cfs = en_cfs.mean('time') - ln_cfs.mean('time')
enso_dif_cfs = enso_dif_cfs.sel(lat=slice(None, None, -1))

# IOD_CFSv2 ---------------------------------------------------------------------------------------------------#
iodp_cfs = (xr.open_dataset(
    cases_dir + 'tref_' + cases_cfs[0] + '_' + s + '.nc')
            .rename({'tref': 'var'})-273)
iodn_cfs = (xr.open_dataset
            (cases_dir + 'tref_' + cases_cfs[1] + '_' + s + '.nc')
            .rename({'tref': 'var'})-273)

iod_cfs_test = ttest_ind(
    iodp_cfs['var'].values, iodn_cfs['var'].values, equal_var=False)[1]

iod_dif_cfs = iodp_cfs.mean('time') - iodn_cfs.mean('time')
iod_dif_cfs = iod_dif_cfs.sel(lat=slice(None, None, -1))

# ENSO_obs --------------------------------------------------------------------#
x, y, en_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[2],
                        two_variables=False)
#en_obs = en_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)
x, y, ln_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[3],
                        two_variables=False)
#ln_obs = ln_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)
enso_obs_test = ttest_ind(
    en_obs['var'].values, ln_obs['var'].values, equal_var=False)[1]

enso_dif_obs = en_obs.mean('time') - ln_obs.mean('time')
#
# IOD_obs ---------------------------------------------------------------------#
x, y, iodp_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[0],
                          two_variables=False)
#iodp_obs = iodp_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)
x, y, iodn_obs = CaseComp(data, s, mmonth=[9, 11], c=cases[1],
                          two_variables=False)
#iodn_obs = iodn_obs.interp(lon=enso_dif_cfs.lon.values, lat=enso_dif_cfs.lat.values)

iod_obs_test = ttest_ind(
    iodp_obs['var'].values, iodn_obs['var'].values, equal_var=False)[1]

iod_dif_obs = iodp_obs.mean('time') - iodn_obs.mean('time')

fig = plt.figure(figsize=(5, 6), dpi=dpi)
gs = fig.add_gridspec(2, 2)

title = ['EN-LN_' + v_name[0].split()[2], 'EN-LN_CFSv2',
         'IOD_pos-neg_'  + v_name[0].split()[2], 'IOD_pos-neg_CFSv2']
comps = [enso_dif_obs, enso_dif_cfs, iod_dif_obs, iod_dif_cfs]
tests = [enso_obs_test, enso_cfs_test, iod_obs_test, iod_cfs_test]

for field, f_count, f_count_c, t, test in zip(comps,[0,1,0,1],[0,0,1,1],
                                              title, tests):
    field = field*MakeMask(field, 'var')
    ax = fig.add_subplot(gs[f_count, f_count_c],
                         projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270, 330, -60, 20], crs=crs_latlon)

    im = ax.contourf(field.lon, field.lat, field['var'],#.where(test<0.5)['var'],
                     levels=scale,
                     transform=crs_latlon, cmap=cbar_t, extend='both')

    # aux = field.where(test < 0.05, np.nan)
    # ax.contourf(field.lon, field.lat, aux['var'],
    #             levels=scale, hatches=['***'],
    #             transform=crs_latlon, color='k')

    cb = plt.colorbar(im, fraction=0.040, pad=0.025, shrink=0.8)
    cb.ax.tick_params(labelsize=5)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                   edgecolor='#4B4B4B')
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color='#4B4B4B', linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=5)
    ax.set_title(t, fontsize=10)
    plt.tight_layout()


fig.suptitle('Temp - ' + s)
if save:
    plt.savefig(out_dir + 'tref_dif_' + s + '_' + v_name[0].split()[2] +
                '.jpg', dpi=dpi)
else:
    plt.show()
#

fig = plt.figure(figsize=(5, 6), dpi=dpi)
gs = fig.add_gridspec(2, 2)

title = ['EN-LN_' + v_name[0].split()[2], 'EN-LN_CFSv2',
         'IOD_pos-neg_'  + v_name[0].split()[2], 'IOD_pos-neg_CFSv2']
comps = [enso_dif_obs, enso_dif_cfs, iod_dif_obs, iod_dif_cfs]
tests = [enso_obs_test, enso_cfs_test, iod_obs_test, iod_cfs_test]

for field, f_count, f_count_c, t, test in zip(comps,[0,1,0,1],[0,0,1,1],
                                              title, tests):
    field = field*MakeMask(field, 'var')
    ax = fig.add_subplot(gs[f_count, f_count_c],
                         projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent([270, 330, -60, 20], crs=crs_latlon)

    im = ax.contourf(field.lon, field.lat, field.where(test<0.1)['var'],
                     levels=scale,
                     transform=crs_latlon, cmap=cbar_t, extend='both')

    # aux = field.where(test < 0.05, np.nan)
    # ax.contourf(field.lon, field.lat, aux['var'],
    #             levels=scale, hatches=['***'],
    #             transform=crs_latlon, color='k')

    cb = plt.colorbar(im, fraction=0.040, pad=0.025, shrink=0.8)
    cb.ax.tick_params(labelsize=5)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey',
                   edgecolor='#4B4B4B')
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color='#4B4B4B', linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
    ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=5)
    ax.set_title(t, fontsize=10)
    plt.tight_layout()


fig.suptitle('Temp - ' + s + ' - sig. 0.1')
if save:
    plt.savefig(out_dir + 'tref_dif_test_' + s + '_' + v_name[0].split()[2] +
                '.jpg', dpi=dpi)
else:
    plt.show()