"""
Anomalía de PP (y T cuando ande) en SA según la magnitud de los índices
(Similar a *2D_bins_DMI_N34_3.0.py, pero en lugar de regiones para todo SA)
"""
########################################################################################################################
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np
import os
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import regionmask
from ENSO_IOD_Funciones import SelectNMMEFiles
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
########################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/' # índices por estaciones
cases_dir = '/pikachu/datos/luciano.andrian/cases_fields/' # campos de las variables PP  ( y T cuando ande)
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Composites/CuadByCases/'
out_data_dir = '/pikachu/datos/luciano.andrian/cases_fields/'

save = False
dpi = 500
# Funciones ############################################################################################################
def SelectBins(data, min, max, sd=1):
    # sd opcional en caso de no estar escalado
    if np.abs(min) > np.abs(max):
        return (data >= min*sd) & (data < max*sd)
    elif np.abs(min) < np.abs(max):
        return (data > min*sd) & (data <= max*sd)
    elif np.abs(min) == np.abs(max):
        return (data >= min*sd) & (data <= max*sd)

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

def BinsByCases(v, v_name, fix_factor, s, mm, c, c_count,
                bin_limits, bins_by_cases_dmi, bins_by_cases_n34, snr=False,
                neutro_clim=False):

    # 1. se abren los archivos de los índices (completos y se pesan por su SD)
    # tambien los archivos de la variable en cuestion pero para cada "case" = c

    data_dates_dmi_or = xr.open_dataset(dates_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc')
    data_dates_dmi_or /=  data_dates_dmi_or.mean('r').std()

    data_dates_n34_or = xr.open_dataset(dates_dir + 'N34_' + s + '_Leads_r_CFSv2.nc')
    data_dates_n34_or /= data_dates_n34_or.mean('r').std()

    # 1.1 Climatología y case
    if neutro_clim:
        clim = xr.open_dataset(cases_dir + v + '_neutros' + '_' + s.upper() + '.nc').rename({v_name: 'var'}) * fix_factor
    else:
        clim = xr.open_dataset(cases_dir + v + '_' + s.lower() + '.nc').rename({v_name: 'var'}) * fix_factor
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

    # 1.2
    anom = anom.sortby(anom.time.dt.month)

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
    # spread = case - comp
    # spread = spread.std('time')

    case_sel_dmi = case_sel_dmi.sortby(case_sel_dmi.time.dt.month)

    case_sel_dmi_n34 = SelectVariables(aux_cases, data_dates_n34_or)
    case_sel_dmi_n34 = case_sel_dmi_n34.sortby(case_sel_dmi_n34.time.dt.month)
    # 2.1 uniendo var, dmi y n34
    data_merged = xr.Dataset(
        data_vars=dict(
            var=(['time', 'lat', 'lon'], anom['var'].values),
            dmi=(['time'], case_sel_dmi.sst.values),
            n34=(['time'], case_sel_dmi_n34.sst.values),
        ),
        coords=dict(
            time=anom.time.values
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

            if snr:
                spread = bin_f - bin_f.mean(['time'])
                spread = spread.std('time')
                SNR = bin_f.mean(['time']) / spread
                anom_bin.append(SNR)
            else:
                anom_bin.append(bin_f.mean('time')['var'])

            num_bin.append(len(np.where(~np.isnan(bin_f['dmi']))[0]))

        anom_bin_main.append(anom_bin)
        num_bin_main.append(num_bin)

    return anom_bin_main, num_bin_main, clim


def fix_calendar(ds, timevar='time'):
    """
    agrega los dias a los archivos nc de NMME
    """
    if ds[timevar].attrs['calendar'] == '360':
        ds[timevar].attrs['calendar'] = '360_day'
    return ds


def DetrendClim(data, mm, v_name='prec'):
    # la diferencia es mínima en fitlrar o no tendencia para hacer una climatología,
    # pero para no perder la costumbre de complicar las cosas...

    for l in [0, 1, 2, 3]:
        season_data = data.sel(time=data.time.dt.month.isin(mm - l), L=l)
        aux = season_data.polyfit(dim='time', deg=1)
        if v_name == 'prec':
            aux_trend = xr.polyval(season_data['time'], aux.prec_polyfit_coefficients[0])  # al rededor de la media
        elif v_name == 'tref':
            aux_trend = xr.polyval(season_data['time'], aux.tref_polyfit_coefficients[0])  # al rededor de la media
        elif v_name == 'hgt':
            aux_trend = xr.polyval(season_data['time'], aux.hgt_polyfit_coefficients[0])  # al rededor de la media



        if l == 0:
            season_anom_detrend = season_data - aux_trend
        else:
            aux_detrend = season_data - aux_trend
            season_anom_detrend = xr.concat([season_anom_detrend, aux_detrend], dim='time')

    return season_anom_detrend.mean(['r', 'time'])


def ComputeFieldsByCases(v, v_name, fix_factor, snr, data,
                         levels_main, cbar_main, levels_clim, cbar_clim,
                         title_var, name_fig, dpi,
                         x_lon=np.arange(280, 330, 10), x_lat=np.arange(-60, 20, 20),
                         figsize=[16,17], usemask=True, hcolorbar=False):
    # no, una genialidad...
    #sec_plot = [7, 2, 22, 17, 13, 14, 10, 11, 8, 9, 3, 4, 20, 21, 15, 16, 23, 24, 18, 19, 5, 6, 0, 1]  # , 12]
    sec_plot = [13, 14, 10, 11,
                7, 2, 22, 17,
                8, 3, 9, 4,
                20, 15, 21, 16,
                5, 0, 6, 1,
                23, 18, 24, 19]
    row_titles = ['Strong El Niño', None, None, None, None,
                  'Moderate El Niño', None, None, None, None,
                  'Neutro ENSO', None, None, None, None,
                  'Moderate La Niña', None, None, None, None,
                  'Strong La Niña', None, None, None, None]
    col_titles = ['Strong IOD - ', 'Moderate IOD - ', 'Neutro IOD', 'Moderate IOD + ', 'Strong IOD + ']
    # ------------------------------------------------------------------------------------------------------------------#

    mm = 7
    for s in ['JJA', 'JAS', 'ASO', 'SON']:
        sec_count = 0
        ct = 0
        cr = 0

        crs_latlon = ccrs.PlateCarree()
        fig, axs = plt.subplots(nrows=5, ncols=5,
                                subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)},
                                figsize=(figsize[0], figsize[1]))
        axs = axs.flatten()
        for c_count in [0, 1, 2, 3, 4, 5, 6, 7]:  # , 8]: # Loop en los cases -{neutro}
            cases_bin, num_bin, aux = BinsByCases(v=v, v_name=v_name, fix_factor=fix_factor,
                                                  s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                                  bin_limits=bin_limits, bins_by_cases_dmi=bins_by_cases_dmi,
                                                  bins_by_cases_n34=bins_by_cases_n34, snr=snr)

            mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(aux)
            mask = xr.where(np.isnan(mask), mask, 1)

            bins_aux_dmi = bins_by_cases_dmi[c_count]
            bins_aux_n34 = bins_by_cases_n34[c_count]
            for b_dmi in range(0, len(bins_aux_dmi)):
                for b_n34 in range(0, len(bins_aux_n34)):
                    n = sec_plot[sec_count]
                    comp_case = cases_bin[b_dmi][b_n34]
                    if usemask:
                        comp_case *= mask
                    if snr:
                        comp_case = comp_case['var']

                    if num_bin[b_dmi][b_n34] > 1:
                        im = axs[n].contourf(aux.lon, aux.lat, comp_case,
                                             levels=levels_main, transform=crs_latlon, cmap=cbar_main, extend='both')
                    else:
                        im = axs[n].contourf(aux.lon, aux.lat, comp_case*0,
                                             levels=levels_main, transform=crs_latlon, cmap=cbar_main, extend='both')

                    axs[n].add_feature(cartopy.feature.LAND, facecolor='lightgrey')
                    axs[n].add_feature(cartopy.feature.COASTLINE)
                    axs[n].gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', color='gray')
                    axs[n].set_xticks(x_lon, crs=crs_latlon)
                    axs[n].set_yticks(x_lat, crs=crs_latlon)
                    lon_formatter = LongitudeFormatter(zero_direction_label=True)
                    lat_formatter = LatitudeFormatter()
                    axs[n].xaxis.set_major_formatter(lon_formatter)
                    axs[n].yaxis.set_major_formatter(lat_formatter)
                    axs[n].tick_params(labelsize=5)
                    if n == 0 or n == 1 or n == 2 or n == 3 or n == 4:
                        axs[n].set_title(col_titles[n], fontsize=15)
                        # ct += 1
                    if n == 0 or n == 5 or n == 10 or n == 15 or n == 20:
                        axs[n].set_ylabel(row_titles[n], fontsize=15)
                        # cr += 1

                    axs[n].xaxis.set_label_position('top')
                    axs[n].set_xlabel('N=' + str(num_bin[b_dmi][b_n34]), fontsize=12, loc='left', fontweight="bold")

                    sec_count += 1
        # en el lugar del neutro -> climatología de la variable (data)

        comp_case_clim = DetrendClim(data, mm, v_name=v_name)

        if usemask:
            comp_case_clim = comp_case_clim[v_name] * fix_factor * mask
        else:
            comp_case_clim = comp_case_clim[v_name] * fix_factor

        if v_name=='hgt':
            comp_case_clim *= 0

        im2 = axs[12].contourf(aux.lon, aux.lat, comp_case_clim,
                               levels=levels_clim, transform=crs_latlon, cmap=cbar_clim, extend='max')
        if v_name != 'hgt':
            cb = plt.colorbar(im2, fraction=0.042, pad=0.032, shrink=1, ax=axs[12], aspect=20)
            cb.ax.tick_params(labelsize=10)

        axs[12].add_feature(cartopy.feature.LAND, facecolor='lightgrey')
        axs[12].add_feature(cartopy.feature.COASTLINE)
        axs[12].gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', color='gray')
        axs[12].set_xticks(x_lon, crs=crs_latlon)
        axs[12].set_yticks(x_lat, crs=crs_latlon)
        lon_formatter = LongitudeFormatter(zero_direction_label=True)
        lat_formatter = LatitudeFormatter()
        axs[12].xaxis.set_major_formatter(lon_formatter)
        axs[12].yaxis.set_major_formatter(lat_formatter)
        axs[12].tick_params(labelsize=5)

        # Colorbar, por ahora vertical. tratando de optimizar el espacio en la verticar de la figura...
        if hcolorbar:
            pos = fig.add_axes([0.2, 0.05, 0.6, 0.01])
            cbar = fig.colorbar(im, cax=pos, pad=0.1, orientation='horizontal')
        else:
            pos = fig.add_axes([0.90, 0.2, 0.012, 0.6])
            cbar = fig.colorbar(im, cax=pos, pad=0.1)

        cbar.ax.tick_params(labelsize=18)
        if snr:
            fig.suptitle('SNR:' + title_var + ' - ' + s, fontsize=20)
        else:
            fig.suptitle(title_var + ' - ' + s, fontsize=20)
        # fig.tight_layout() #BUG matplotlib 3.5.0 #Solucionado definitivamnete en 3.6 ?
        if hcolorbar:
            fig.subplots_adjust(top=0.93, bottom=0.07)
        else:
            fig.subplots_adjust(top=0.93)
        if snr:
            plt.savefig(out_dir + 'SNR_' + name_fig + '_' + s + '.jpg', bbox_inches='tight', dpi=dpi)
        else:
            plt.savefig(out_dir + name_fig + '_' + s + '.jpg', bbox_inches='tight', dpi=dpi)

        plt.close('all')
        mm += 1


########################################################################################################################
cases = ['dmi_puros_pos', 'dmi_puros_neg',
        'n34_puros_pos', 'n34_puros_neg',
        'sim_pos', 'sim_neg',
        'dmi_neg_n34_pos', 'dmi_pos_n34_neg',
        'neutros']

bin_limits = [[-4.5,-1], [-1, -0.5], #1
              [-0.5, 0.5], #2
              [0.5, 1], [1, 4.5]] #4

bins_by_cases_dmi = [[3, 4], [0, 1],
                     [2], [2],
                     [3, 4], [0, 1],
                     [0, 1],[3, 4],
                     [2]]

bins_by_cases_n34 = [[2], [2],
                     [3, 4], [0,1],
                     [3, 4], [0, 1],
                     [3, 4], [0, 1],
                     [2]]

cbar_pp = colors.ListedColormap(['#003C30', '#004C42', '#0C7169', '#79C8BC', '#B4E2DB',
                                'white',
                                '#F1DFB3', '#DCBC75', '#995D13', '#6A3D07', '#543005', ][::-1])
cbar_pp.set_under('#3F2404')
cbar_pp.set_over('#00221A')
cbar_pp.set_bad(color='white')

cbar_snr = colors.ListedColormap(['#002A3D','#074D4F', '#1E6D5A' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#DCBC75', '#995D13','#6A3D07','#543005','#3F2404'][::-1])
cbar_snr.set_under('#3F2404')
cbar_snr.set_over('#002A3D')
cbar_snr.set_bad(color='white')

cbar_t = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar_t.set_over('#9B1C00')
cbar_t.set_under('#014A9B')
cbar_t.set_bad(color='white')

cbar_snr_t = colors.ListedColormap(['#070B4F','#2E07AC', '#387AE4' ,'#52C39D','#6FFE9B',
                                  '#FFFFFF',
                                  '#FEB77E', '#FB8761','#CA3E72','#782281','#251255'])
cbar_snr_t.set_over('#251255')
cbar_snr_t.set_under('#070B4F')
cbar_snr_t.set_bad(color='white')

cbar_hgt = colors.ListedColormap(['#9B1C00','#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF', '#014A9B'][::-1])
cbar_hgt.set_over('#641B00')
cbar_hgt.set_under('#012A52')
cbar_hgt.set_bad(color='white')

########################################################################################################################
# Prec #################################################################################################################
# HINDCAST para climatología ------------------------------------------------------------------------------------------#
try:
    data = xr.open_dataset(out_data_dir + 'prec_data_full.nc')
except:
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable='prec',
                            dir='/pikachu/datos/osman/nmme/monthly/hindcast/', All=True)
    files = sorted(files, key=lambda x: x.split()[0])

    # abriendo todos los archivos
    data = xr.open_mfdataset(files, decode_times=False)  # xr no entiende la codificacion de Leads, r y las fechas
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
    data = data.sel(r=slice(1, 24))
    data = data.rolling(time=3, center=True).mean()
    #data.to_netcdf(out_data_dir + 'prec_data_full.nc')
    data = data.compute()
#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# Composite
levels = np.linspace(-30, 30, 13)
levels_clim = np.linspace(0,250,11)


ComputeFieldsByCases(v='prec', v_name='prec', fix_factor=30, snr=False,
                     data=data,
                     levels_main=levels, cbar_main=cbar_pp, 
                     levels_clim=levels_clim, cbar_clim='YlGnBu',
                     title_var='Prec', name_fig='prec', dpi=500)

# Signal-to-Noise ratio
levels = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]

ComputeFieldsByCases(v='prec', v_name='prec', fix_factor=30, snr=True,
                     data=data,
                     levels_main=levels, cbar_main=cbar_snr,
                     levels_clim=levels_clim, cbar_clim='YlGnBu',
                     title_var='Prec', name_fig='prec', dpi=500)

########################################################################################################################
# Temp #################################################################################################################
# HINDCAST para climatología ------------------------------------------------------------------------------------------#
try:
    data = xr.open_dataset(out_data_dir + 'temp_data_full.nc')
except:
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable='tref',
                            dir='/pikachu/datos/osman/nmme/monthly/hindcast/', All=True)
    files = sorted(files, key=lambda x: x.split()[0])

    # abriendo todos los archivos
    data = xr.open_mfdataset(files, decode_times=False)  # xr no entiende la codificacion de Leads, r y las fechas
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lon=slice(275, 330), lat=slice(-60, 15))
    data = data.sel(r=slice(1, 24))

    # media movil de 3 meses para separar en estaciones
    data = data.rolling(time=3, center=True).mean()
    #data = data.to_netcdf(out_data_dir + 'temp_data_full.nc')
    data = data.compute()

#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#
# Composite
levels =  np.linspace(-1.2,1.2,13)
levels_clim = np.linspace(0,25,11)

ComputeFieldsByCases(v='tref', v_name='tref', fix_factor=1, snr=False,
                     data=data-273,
                     levels_main=levels, cbar_main=cbar_t,
                     levels_clim=levels_clim, cbar_clim='YlOrRd',
                     title_var='Temp', name_fig='temp', dpi=500)

# Signal-to-Noise ratio
levels = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]

ComputeFieldsByCases(v='tref', v_name='tref', fix_factor=1, snr=True,
                     data=data-273,
                     levels_main=levels, cbar_main=cbar_snr_t,
                     levels_clim=levels_clim, cbar_clim='YlOrRd',
                     title_var='Tref', name_fig='tref', dpi=500)

########################################################################################################################
# HGT #################################################################################################################
# HINDCAST para climatología ------------------------------------------------------------------------------------------#

try:
    data = xr.open_dataset(out_data_dir + 'hgt_data_full.nc')
except:
    files = SelectNMMEFiles(model_name='NCEP-CFSv2', variable='hgt',
                            dir='/pikachu/datos/luciano.andrian/hindcast/', All=True)
    files = sorted(files, key=lambda x: x.split()[0])

    # abriendo todos los archivos
    data = xr.open_mfdataset(files, decode_times=False)  # xr no entiende la codificacion de Leads, r y las fechas
    data = data.rename({'X': 'lon', 'Y': 'lat', 'M': 'r', 'S': 'time'})
    data = data.sel(L=[0.5, 1.5, 2.5, 3.5])  # Solo leads 0 1 2 3
    data['L'] = [0, 1, 2, 3]
    data = xr.decode_cf(fix_calendar(data))  # corrigiendo fechas
    data = data.sel(lon=slice(30, 340), lat=slice(-80, 20), P=200)
    data = data.drop('P')
    data = data.sel(r=slice(1, 24))

    # media movil de 3 meses para separar en estaciones
    data = data.rolling(time=3, center=True).mean()
    #data = data.to_netcdf(out_data_dir + 'hgt_data_full.nc')
    data = data.compute()

#----------------------------------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------------------------------#

levels = [-300,-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250,300]
levels_clim = np.linspace(10000,15000,11)

ComputeFieldsByCases(v='hgt', v_name='hgt', fix_factor=9.8, snr=False,
                     data=data*9.8,
                     levels_main=levels, cbar_main=cbar_hgt,
                     levels_clim=levels_clim, cbar_clim='Spectral',
                     title_var='HGT200_NC', name_fig='hgt200_neutro_clim', dpi=500,
                     x_lon=np.arange(30, 340, 25), x_lat=np.arange(-80, 20, 10),
                     figsize=[20,12], usemask=False, hcolorbar=True)

# Signal-to-Noise ratio
levels = [-1,-.8,-.6,-.4,-.2,-.1,0,0.1,0.2,0.4,0.6,0.8,1]

ComputeFieldsByCases(v='hgt', v_name='hgt', fix_factor=9.8, snr=True,
                     data=data,
                     levels_main=levels, cbar_main=cbar_snr_t,
                     levels_clim=levels_clim, cbar_clim='YlGnBu',
                     title_var='HGT200_NC', name_fig='hgt200_neutro_clim', dpi=500,
                     x_lon=np.arange(30, 340, 25), x_lat=np.arange(-80, 20, 10),
                     figsize=[20, 13], usemask=False, hcolorbar=True)
########################################################################################################################