"""
Composites observado con filtro de tendencia corregida en las variables T, PP y las de ERA5
T, PP y HGT con sig. MC
SST, VP y dig sin sig.
"""
########################################################################################################################
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
import warnings
warnings.filterwarnings("ignore")
########################################################################################################################
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/' #fechas
data_dir_t_pp = '/pikachu/datos/luciano.andrian/observado/ncfiles/data_obs_d_w_c/' #T y PP ya procesados
data_dir_era5 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/' # ERA5 ya procesados
out_dir_w_sig = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/w_sig/DMIbase/'
out_dir_no_sig = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/no_sig/DMIbase/'

#Plot
save = True
dpi = 300
sig = True
sig_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_quantiles/' # resultados de MC
# Functions ############################################################################################################
def CompositeSimple(original_data, index, mmin, mmax):
    def is_months(month, mmin, mmax):
        return (month >= mmin) & (month <= mmax)

    if len(index) != 0:
        comp_field = original_data.sel(time=original_data.time.dt.year.isin([index]))
        comp_field = comp_field.sel(
            time=is_months(month=comp_field['time.month'], mmin=mmin, mmax=mmax))
        if len(comp_field.time) != 0:
            comp_field = comp_field.mean(['time'], skipna=True)
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
        case = case.where(case >= 1950)
        aux.close()

        case_num = len(case.values[np.where(~np.isnan(case.values))])

        neutro_comp = CompositeSimple(original_data=data, index=neutro, mmin=mmin, mmax=mmax)
        data_comp = CompositeSimple(original_data=data, index=case, mmin=mmin, mmax=mmax)

        comp = data_comp - neutro_comp

        if two_variables:
            neutro_comp2 = CompositeSimple(original_data=data2, index=neutro, mmin=mmin, mmax=mmax)
            data_comp2 = CompositeSimple(original_data=data2, index=case, mmin=mmin, mmax=mmax)

            comp2 = data_comp2 - neutro_comp2
        else:
            comp2 = None
    except:
        print('Error en ' + s + c)

    if two_variables:
        return comp, case_num, comp2
    else:
        return comp, case_num


def Plot(comp, levels, cmap, step1, contour1=True,
         two_variables=False, comp2=None, levels2=np.linspace(-1,1,13), step2=4,
         mapa='sa', title='title', name_fig='name_fig', dpi=100, save=save,
         comp_sig=None, color_sig='k', significance=True, linewidht2=.5, color_map='#4B4B4B',
         out_dir=out_dir_w_sig):

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

    else:
        fig_size = (8, 3)
        extent = [30, 330, -80, 20]
        xticks = np.arange(30, 330, 30)
        yticks = np.arange(-80, 20, 10)

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

    fig = plt.figure(figsize=fig_size, dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    ax.set_extent(extent, crs=crs_latlon)
    if two_variables:
        ax.contour(comp2.lon[::step2], comp2.lat[::step2], comp_var2[::step2, ::step2],
                   linewidths=linewidht2, levels=levels_contour2, transform=crs_latlon, colors='k')
    else:
        if contour1:
            ax.contour(comp.lon[::step1], comp.lat[::step1], comp_var[::step1, ::step1],
                       linewidths=.8, levels=levels_contour, transform=crs_latlon, colors='black')

    im = ax.contourf(comp.lon[::step1], comp.lat[::step1], comp_var[::step1, ::step1],
                     levels=levels, transform=crs_latlon, cmap=cmap, extend='both')

    if significance:
        colors_l = [color_sig, color_sig]
        comp_sig_var = comp_sig['var']
        cs = ax.contourf(comp_sig.lon, comp_sig.lat, comp_sig_var,
                         transform=crs_latlon, colors='none', hatches=["..", ".."], extend='lower')
        for i, collection in enumerate(cs.collections):
            collection.set_edgecolor(colors_l[i % len(colors_l)])

        for collection in cs.collections:
            collection.set_linewidth(0.)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='lightgrey', edgecolor=color_map)
    ax.add_feature(cartopy.feature.COASTLINE, linewidth=0.5)
    ax.coastlines(color=color_map, linestyle='-', alpha=1)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    ax.set_xticks(xticks, crs=crs_latlon)
    ax.set_yticks(yticks, crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)
    plt.title(title, fontsize=10)
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()

########################################################################################################################
# seasons = ('JJA', 'JAS', 'ASO', 'SON')
# min_max_months = [[6,8],[7,9],[8,10],[9,11]]

# seasons = ('SON')
# min_max_months = [[9,11]]

seasons = ('JJA', 'SON')
min_max_months = [[6,8],[9,11]]

variables_t_p = ['t_cru_d_w_c_1950-2020_0.25.nc', 'pp_gpcc_d_w_c_1950-2020_0.25.nc', 'pp_prec_d_w_c_1950-2020_2.5.nc']
variables_ERA5 = ['hgt200_mer_d_w', 'div200_mer_d_w', 'vp200_mer_d_w']

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos',
         'DMI_un_neg','N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

scales = [np.linspace(-1, 1 ,17),  #t
          [-300, -250, -200, -150, -100, -50, -25, 0, 25, 50, 100, 150, 200, 250, 300],  # hgt
          np.linspace(-45, 45, 13),  #pp
          [-300,-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250,300],  #hgt
          np.linspace(-45, 45, 13),  # pp
          [-300, -250, -200, -150, -100, -50, -25, 0, 25, 50, 100, 150, 200, 250, 300],  # hgt
          np.linspace(-0.5e-5, 0.5e-5, 17)]  #div200

title_case = ['DMI-ENSO simultaneous positive phase ',
              'DMI-ENSO simultaneous negative phase ',
              'DMI negative phase ',
              'DMI positive phase ',
              'DMI pure positive phase ',
              'DMI pure negative phase ',
              'ENSO positive phase ',
              'ENSO negative phase ',
              'ENSO pure positive phase ',
              'ENSO pure negative phase ']

v_name = [ 'Temperature - Cru', 'Precipitation - GPCC', 'Precipitation - PREC',
           'HGT200hPa','Div200hPa', 'Potential Velocity']
v_name_fig=['temp', 'pp_gpcc', 'pp_prec','hgt200', 'div', 'pv']

# colorbars
cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')


cbar_t = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar_t.set_over('#691800')
cbar_t.set_under('#013774')
cbar_t.set_bad(color='white')

cbar_t_r = colors.ListedColormap(['#9B1C00', '#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'])
cbar_t_r.set_under('#691800')
cbar_t_r.set_over('#013774')
cbar_t_r.set_bad(color='white')

cbar_sst = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_sst.set_over('#9B1C00')
cbar_sst.set_under('#014A9B')
cbar_sst.set_bad(color='white')

cmap_t_pp = [cbar_t, cbar_pp, cbar_pp]
cmap_era5 = [cbar_t, cbar_t_r]
########################################################################################################################
#T y PP con contornos de HGT200
v_count = 0
plt.rcParams['hatch.linewidth'] = 2
for v in variables_t_p:
    data = xr.open_dataset(data_dir_t_pp + v)
    data2 = xr.open_dataset(data_dir_era5 + variables_ERA5[0] + '.nc')
    data2 = data2.sel(lon=slice(270, 330), lat=slice(15, -60))
    #data2 = data2.interp(lon=data.lon.values, lat=data.lat.values)

    c_count = 0
    for c in cases:
        s_count = 0
        for s in seasons:
            comp1, num_case, comp2 = CaseComp(data, s, mmonth=min_max_months[s_count], c=c,
                                              two_variables=True, data2=data2)

            data_sig = xr.open_dataset(sig_dir + v.split('_')[0] + '_' + v.split('_')[1] +
                                       '_' + c + '1950_2020_' + s + '_DMIbase.nc')

            comp1_i=comp1.interp(lon=data_sig.lon.values, lat=data_sig.lat.values)
            sig = comp1_i.where((comp1_i < data_sig['var'][0]) | (comp1_i > data_sig['var'][1]))
            sig = sig.where(np.isnan(sig['var']), 0)

            if v_count != 0:
                v_count_sc = 2
            else:
                v_count_sc = 0

            #MakeMask(pp, dataname='cluster')
            Plot(comp=comp1, levels=scales[v_count_sc], cmap = cmap_t_pp[v_count], step1=1, contour1=False,
                 two_variables=True, comp2=comp2, levels2=scales[v_count_sc + 1], step2=4,
                 mapa='sa', significance=True,
                 title=v_name[v_count] + '\n' + title_case[c_count] + '\n' + s + ' - Events: ' + str(num_case) ,
                 name_fig=v_name_fig[v_count] + s + '_' + cases[c_count] + '_mer_d_w_NSA',
                 dpi=dpi, save=save, comp_sig=sig, color_sig='k')

            s_count += 1
        c_count += 1
    v_count += 1

# HGT -----------------------------------------------------------------------------------------------------------------#
plt.rcParams['hatch.linewidth'] = 1.5
tw=[False, False] # por ahora sin vp
sig2 = [True, False]
steps = [1, 4]
contours1 = [True, False]
v_count = 2
for v in variables_ERA5:
    if v_count != 2:
        break #provisorio
    data = xr.open_dataset(data_dir_era5 + v + '.nc')

    c_count = 0
    for c in cases:
        s_count = 0
        for s in seasons:
            comp1, num_case = CaseComp(data, s, mmonth=min_max_months[s_count], c=c,
                                              two_variables=False, data2=None)

            data_sig = xr.open_dataset(sig_dir +  v.split('_')[0] +
                                       '_' + c + '1950_2020_' + s + '_DMIbase.nc')

            comp1_i=comp1.interp(lon=data_sig.lon.values, lat=data_sig.lat.values)
            sig = comp1_i.where((comp1_i < data_sig['var'][0]) | (comp1_i > data_sig['var'][1]))
            sig = sig.where(np.isnan(sig['var']), 0)

            Plot(comp=comp1, levels=scales[v_count + 1], cmap = cmap_era5[v_count-2], step1=steps[v_count-2],
                 contour1=contours1[v_count-2], two_variables=False,
                 mapa='hs', significance=True,
                 title=v_name[v_count] + '\n' + title_case[c_count] + '\n' + s + ' - Events: ' + str(num_case) ,
                 name_fig=v_name_fig[v_count]  + s + '_' + cases[c_count] + '_mer_d_w_NSA',
                 dpi=dpi, save=save, comp_sig=sig, color_sig='k')

            s_count += 1
        c_count += 1
    v_count += 1

# SST -----------------------------------------------------------------------------------------------------------------#
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

    # ax.set_xticks(np.arange(50, 270, 30), crs=crs_latlon)
    # ax.set_yticks(np.arange(-20, 20, 10), crs=crs_latlon)
v_count = 0
v = 'sst'
data = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
data = data.rename({'sst':'var'})
data = Detrend(data, 'time')
data = data.sel(lon=slice(50, 270), lat=slice(20, -20))

c_count = 0
for c in cases:
    s_count = 0
    for s in seasons:
        comp1, num_case = CaseComp(data, s, mmonth=min_max_months[s_count], c=c,
                                          two_variables=False)

        Plot(comp=comp1, levels=[-1.5, -1, -0.5, -0.25, 0, 0.25, 0.5, 1, 1.5], cmap=cbar_sst, step1=1, contour1=False,
             two_variables=False, mapa='tropical', significance=False,
             title='SST' + '\n' + title_case[c_count] + '\n' + s + ' - Events: ' + str(num_case),
             name_fig='SST_' + s + '_' + cases[c_count] + '_d_NSA',
             dpi=dpi, save=save, out_dir=out_dir_no_sig)

        s_count += 1
    c_count += 1

# Vp y Div ------------------------------------------------------------------------------------------------------------#
v_from_w = ['div_from_w', 'vp_from_w'] # creadas a partir de windphsere
def Detrend(xrda, dim):
    aux = xrda.polyfit(dim=dim, deg=1)
    try:
        trend = xr.polyval(xrda[dim], aux.var_polyfit_coefficients)
    except:
        trend = xr.polyval(xrda[dim], aux.polyfit_coefficients)
    dt = xrda - trend
    return dt

def OrdenarNC_wTime_fromW(data):

    newdata = xr.Dataset(
        data_vars=dict(
            var=(['time', 'lat', 'lon'], data['var'][0, :, :, :].values)

        ),
        coords=dict(
            lon=(['lon'], data.lon),
            lat=(['lat'], data.lat),
            time=(['time'], data.time)
        )
    )

    return newdata

data1 = xr.open_dataset(data_dir_era5 + v_from_w[0] + '.nc')
data1 = Detrend(OrdenarNC_wTime_fromW(data1.rename({'divergence':'var'})), 'time')
data2 = xr.open_dataset(data_dir_era5 + v_from_w[1] + '.nc')
data2 = Detrend(OrdenarNC_wTime_fromW(data2.rename({'velocity_potential':'var'})), 'time')

c_count = 0
for c in cases:
    s_count = 0
    for s in seasons:
        comp1, num_case = CaseComp(data1, s, mmonth=min_max_months[s_count], c=c,
                                   two_variables=False, data2=None)

        comp2, num_case = CaseComp(data2, s, mmonth=min_max_months[s_count], c=c,
                                   two_variables=False, data2=None)

        Plot(comp=comp1, levels=np.linspace(-0.5e-5, 0.5e-5, 13), cmap=cbar_sst, step1=1, contour1=True,
             two_variables=True, comp2=comp2, levels2=np.linspace(-4.5e6, 4.5e6, 13), significance=False,
             mapa='HS',
             title='Div200hpa [shade] - VP [cont.]' + '\n' + title_case[c_count] + '\n' + s + ' - Events: ' + str(num_case),
             name_fig='divp_' + s + '_' + cases[c_count] + '_d_NSA',
             dpi=dpi, save=save, linewidht2=.8, out_dir=out_dir_no_sig)



        s_count += 1
    c_count += 1

########################################################################################################################