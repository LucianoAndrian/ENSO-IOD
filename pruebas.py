from itertools import groupby
import xarray as xr
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import statsmodels.formula.api as sm
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import regionmask
import matplotlib.pyplot as plt
import matplotlib.path as mpath

import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

out_dir = '~/'
levels = [-300,-250,-200,-150,-100,-50,-25,0,25,50,100,150,200,250,300]
levels_clim = np.linspace(10000,15000,11)
v='hgt'; v_name='hgt'; fix_factor=9.8; snr=False; data=None; cases=cases; bin_limits=bin_limits;
bins_by_cases_dmi=bins_by_cases_dmi;bins_by_cases_n34=bins_by_cases_n34;cases_dir=cases_dir; dates_dir=dates_dir;levels_main=levels; cbar_main=cbar_hgt;levels_clim=levels_clim; cbar_clim='Spectral';title_var='HGT200_NC'; name_fig='hgt200_neutro_clim'; dpi=100;x_lon=np.arange(30, 340, 25); x_lat=np.arange(-80, 20, 10);figsize=[20,10]; usemask=False; hcolorbar=True; save=save;out_dir=out_dir

from ENSO_IOD_Funciones import SelectVariables

def SelectBins(data, min, max, sd=1):
    # sd opcional en caso de no estar escalado
    if np.abs(min) > np.abs(max):
        return (data >= min*sd) & (data < max*sd)
    elif np.abs(min) < np.abs(max):
        return (data > min*sd) & (data <= max*sd)
    elif np.abs(min) == np.abs(max):
        return (data >= min*sd) & (data <= max*sd)

def BinsByCases(v, v_name, fix_factor, s, mm, c, c_count,
                bin_limits, bins_by_cases_dmi, bins_by_cases_n34, dates_dir, cases_dir,
                snr=False, neutro_clim=False, obsdates=False):

    # 1. se abren los archivos de los índices (completos y se pesan por su SD)
    # tambien los archivos de la variable en cuestion pero para cada "case" = c

    data_dates_dmi_or = xr.open_dataset(dates_dir + 'DMI_' + s + '_Leads_r_CFSv2.nc')
    data_dates_dmi_or /=  data_dates_dmi_or.mean('r').std()

    data_dates_n34_or = xr.open_dataset(dates_dir + 'N34_' + s + '_Leads_r_CFSv2.nc')
    data_dates_n34_or /= data_dates_n34_or.mean('r').std()

    # 1.1 Climatología y case
    end_nc_file = '.nc' if v != 'tref' else '_nodetrend.nc'

    if neutro_clim:
        clim = xr.open_dataset(cases_dir + v + '_neutros' + '_' + s.upper() + end_nc_file).rename({v_name: 'var'}) * fix_factor
    else:
        clim = xr.open_dataset(cases_dir + v + '_' + s.lower() + end_nc_file).rename({v_name: 'var'}) * fix_factor

    case = xr.open_dataset(cases_dir + v + '_' + c + '_' + s.upper() + end_nc_file).rename({v_name: 'var'}) * fix_factor

    # Anomalía
    for l in [0, 1, 2, 3]:
        try:
            clim_aux = clim.sel(time=clim.time.dt.month.isin(mm - l)).mean(['r', 'time'])
        except:
            clim_aux = clim.sel(time=clim.time.dt.month.isin(mm - l)).mean(['time'])

        if l==0:
            anom = case.sel(time=case.time.dt.month.isin(mm - l)) - clim_aux
        else:
            anom2 = case.sel(time=case.time.dt.month.isin(mm - l)) - clim_aux
            anom = xr.concat([anom, anom2], dim='time')

    # 1.2
    anom = anom.sortby(anom.time.dt.month)

    # 2. Vinculo fechas case -> índices DMI y N34 para poder clasificarlos
    # las fechas entre el case variable y el case indices COINCIDEN,
    # DE ESA FORMA SE ELIGIERON LOS CASES VARIABLE
    # pero diferen en orden. Para evitar complicar la selección usando r y L
    # con .sortby(..time.dt.month) en cada caso se simplifica el problema
    # y coinciden todos los eventos en fecha, r y L

    if obsdates: # por ahora no funca
        print('Fechas Observadas deshabilitado')
        return
        # aux_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
        # aux_cases = xr.open_dataset(aux_dir + v + '_' + c + '_' + s + '_CFSv2_obsDates.nc')\
        #     .rename({v: 'index'})
        # aux_cases['index'] = aux_cases.time
        # aux_cases = aux_cases.drop(['lon', 'lat'])
    else:
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

def ComputeFieldsByCases(v, v_name, fix_factor, snr, data,
                         levels_main, cbar_main, levels_clim, cbar_clim,
                         title_var, name_fig, dpi,
                         cases, bin_limits, bins_by_cases_dmi, bins_by_cases_n34,
                         cases_dir, dates_dir,
                         figsize=[16,17], usemask=True, hcolorbar=False, save=True,
                         proj='eq', obsdates=False,  out_dir='~/'):

    # no, una genialidad... -------------------------------------------------------------------------------------------#
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
    num_neutros = [483, 585, 676, 673]
    porcentaje = 0.1
    # ------------------------------------------------------------------------------------------------------------------#
    print('Only SON')
    print('No climatology')
    mm = 10
    for s in ['SON']:
        n_check = []
        sec_count = 0
        # esto no tiene sentido
        # comp_case_clim = DetrendClim(data, mm, v_name=v_name)

        crs_latlon = ccrs.PlateCarree()
        if proj=='eq':
            fig, axs = plt.subplots(nrows=5, ncols=5,
                                    subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)},
                                    figsize=(figsize[0], figsize[1]))
        else:
            fig, axs = plt.subplots(nrows=5, ncols=5,
                                    subplot_kw={'projection': ccrs.SouthPolarStereo(central_longitude=200)},
                                    figsize=(figsize[0], figsize[1]))
        axs = axs.flatten()
        # Loop en los cases -{neutro} ---------------------------------------------------------------------------------#
        for c_count in [0, 1, 2, 3, 4, 5, 6, 7]:  # , 8]:
            cases_bin, num_bin, aux = BinsByCases(v=v, v_name=v_name, fix_factor=fix_factor,
                                                  s=s, mm=mm, c=cases[c_count], c_count=c_count,
                                                  bin_limits=bin_limits, bins_by_cases_dmi=bins_by_cases_dmi,
                                                  bins_by_cases_n34=bins_by_cases_n34, snr=snr,
                                                  cases_dir=cases_dir, dates_dir=dates_dir, obsdates=obsdates)

            bins_aux_dmi = bins_by_cases_dmi[c_count]
            bins_aux_n34 = bins_by_cases_n34[c_count]
            for b_dmi in range(0, len(bins_aux_dmi)):
                for b_n34 in range(0, len(bins_aux_n34)):
                    n = sec_plot[sec_count]
                    if proj != 'eq':
                        axs[n].set_extent([0, 360, -80, 20],
                                          ccrs.PlateCarree(central_longitude=200))
                    comp_case = cases_bin[b_dmi][b_n34]

                    # if v == 'prec' and s == 'JJA':
                    #
                    #     mask2 = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(comp_case_clim)
                    #     mask2 = xr.where(np.isnan(mask2), mask2, 1)
                    #     mask2 = mask2.to_dataset(name='prec')
                    #
                    #     dry_season_mask = comp_case_clim.where(comp_case_clim.prec>30)
                    #     dry_season_mask = xr.where(np.isnan(dry_season_mask), dry_season_mask, 1)
                    #     dry_season_mask *= mask2
                    #
                    #     if snr:
                    #         comp_case['var'] *= dry_season_mask.prec
                    #     else:
                    #         comp_case *= dry_season_mask.prec.values

                    if usemask:
                        mask = regionmask.defined_regions.natural_earth_v5_0_0.countries_110.mask(aux)
                        mask = xr.where(np.isnan(mask), mask, 1)
                        comp_case *= mask
                    if snr:
                        comp_case = comp_case['var']

                    if num_bin[b_dmi][b_n34] > num_neutros[mm-7]*porcentaje:
                        im = axs[n].contourf(aux.lon, aux.lat, comp_case,
                                             levels=levels_main, transform=crs_latlon,
                                             cmap=cbar_main, extend='both')

                        axs[n].add_feature(cartopy.feature.LAND, facecolor='lightgrey')
                        axs[n].add_feature(cartopy.feature.COASTLINE)
                        if proj == 'eq':
                            axs[n].gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', color='gray')
                            axs[n].set_xticks([])
                            axs[n].set_yticks([])
                            # axs[n].set_xticks(x_lon, crs=crs_latlon)
                            # axs[n].set_yticks(x_lat, crs=crs_latlon)
                            # lon_formatter = LongitudeFormatter(zero_direction_label=True)
                            # lat_formatter = LatitudeFormatter()
                            # axs[n].xaxis.set_major_formatter(lon_formatter)
                            # axs[n].yaxis.set_major_formatter(lat_formatter)
                        else:
                            #polar
                            gls = axs[n].gridlines(draw_labels=True, crs=crs_latlon, lw=0.3, color="gray",
                                                   y_inline=True, xlocs=range(-180, 180, 30),
                                                   ylocs=np.arange(-80, 20, 20))
                            r_extent = 1.2e7
                            axs[n].set_xlim(-r_extent, r_extent)
                            axs[n].set_ylim(-r_extent, r_extent)
                            circle_path = mpath.Path.unit_circle()
                            circle_path = mpath.Path(circle_path.vertices.copy() * r_extent,
                                                     circle_path.codes.copy())
                            axs[n].set_boundary(circle_path)
                            axs[n].set_frame_on(False)
                            plt.draw()
                            for ea in gls._labels:
                                pos = ea[2].get_position()
                                if (pos[0] == 150):
                                    ea[2].set_position([0, pos[1]])

                        axs[n].tick_params(labelsize=0)

                        if n == 0 or n == 1 or n == 2 or n == 3 or n == 4:
                            axs[n].set_title(col_titles[n], fontsize=15)

                        if n == 0 or n == 5 or n == 10 or n == 15 or n == 20:
                            axs[n].set_ylabel(row_titles[n], fontsize=15)

                        axs[n].xaxis.set_label_position('top')
                        axs[n].set_xlabel('N=' + str(num_bin[b_dmi][b_n34]), fontsize=12, loc='left', fontweight="bold")

                    else:
                        n_check.append(n)
                        axs[n].axis('off')

                    sec_count += 1

        # subtitulos columnas de no ploteados -------------------------------------------------------------------------#
        for n_aux in [0, 1, 2, 3, 4]:
            if n_aux in n_check:
                if n_aux+5 in n_check:
                    axs[n_aux+10].set_title(col_titles[n_aux], fontsize=15)
                else:
                    axs[n_aux+5].set_title(col_titles[n_aux], fontsize=15)

        for n_aux in [0, 5, 10, 15, 20]:
            if n_aux in n_check:
                if n_aux+1 in n_check:
                    axs[n_aux+2].set_ylabel(row_titles[n_aux], fontsize=15)
                else:
                    axs[n_aux+1].set_ylabel(row_titles[n_aux], fontsize=15)

        # Climatologia = NADA en HGT ----------------------------------------------------------------------------------#
        # en el lugar del neutro -> climatología de la variable (data)

        # if usemask:
        #     comp_case_clim = comp_case_clim[v_name] * fix_factor * mask
        # else:
        #     comp_case_clim = comp_case_clim[v_name] * fix_factor

        # if v_name=='hgt':
        #     comp_case_clim = 0

        aux0 = aux.sel(r=1, time='1982-10-01').drop(['r','L','time'])
        im2 = axs[12].contourf(aux.lon, aux.lat, aux0['var'][0,:,:],
                               levels=levels_clim, transform=crs_latlon, cmap=cbar_clim, extend='max')
        #--------------------------------------------------------------------------------------------------------------#
        axs[12].add_feature(cartopy.feature.LAND, facecolor='grey')
        axs[12].add_feature(cartopy.feature.COASTLINE)

        if v_name != 'hgt':
            cb = plt.colorbar(im2, fraction=0.042, pad=0.032, shrink=1, ax=axs[12], aspect=20)
            cb.ax.tick_params(labelsize=5)

        if proj == 'eq':
            axs[12].gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-', color='gray')
            axs[12].set_xticks([])
            axs[12].set_yticks([])
        else:
            #polar
            gls = axs[12].gridlines(draw_labels=True, crs=crs_latlon, lw=0.3, color="gray",
                                   y_inline=True, xlocs=range(-180, 180, 30), ylocs=np.arange(-80, 0, 20))
            r_extent = 1.2e7
            axs[12].set_xlim(-r_extent, r_extent)
            axs[12].set_ylim(-r_extent, r_extent)
            circle_path = mpath.Path.unit_circle()
            circle_path = mpath.Path(circle_path.vertices.copy() * r_extent,
                                     circle_path.codes.copy())
            axs[12].set_boundary(circle_path)
            axs[12].set_frame_on(False)
            axs[12].tick_params(labelsize=0)
            plt.draw()
            for ea in gls._labels:
                pos = ea[2].get_position()
                if (pos[0] == 150):
                    ea[2].set_position([0, pos[1]])

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
        if save:
            if snr:
                plt.savefig(out_dir + 'SNR_' + name_fig + '_' + s + '.jpg', bbox_inches='tight', dpi=dpi)
            else:
                plt.savefig(out_dir + name_fig + '_' + s + '.jpg', bbox_inches='tight', dpi=dpi)

            plt.close('all')
        else:
            plt.show()
        mm += 1
########################################################################################################################