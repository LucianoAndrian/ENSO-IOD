import xarray as xr
from ENSO_IOD_Funciones import WAF
from ENSO_IOD_Funciones import PlotWAFCountours
########################################################################################################################
data_dir2 = '/pikachu/datos/luciano.andrian/observado/ncfiles/ERA5/mer_d_w/'
nc_date_dir = '/pikachu/datos/luciano.andrian/observado/ncfiles/nc_composites_dates_no_ind_sst_anom/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/WAF/no_sstanoms/'

save=True
## Functions ###########################################################################################################
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

def CompositeOfWAF(original_data, index, mmin, mmax, neutro_comp):
    def is_months(month, mmin, mmax):
        return (month >= mmin) & (month <= mmax)

    if len(index) != 0:
        comp_field = original_data.sel(time=original_data.time.dt.year.isin([index]))
        comp_field = comp_field.sel(
            time=is_months(month=comp_field['time.month'], mmin=mmin, mmax=mmax))

        #WAF
        comp_field = comp_field.groupby('time.year').mean()
        comp_field = comp_field - neutro_comp

        y_count=0
        for y in comp_field.year.values:
            comp = comp_field.sel(year=y)

            px, py = WAF(neutro_comp, comp, neutro_comp.lon, neutro_comp.lat, reshape=True, variable='var')

            if y_count==0:
                px_xr = xr.DataArray(px)
                py_xr = xr.DataArray(py)
            else:
                px_xr_next = xr.DataArray(px)
                py_xr_next = xr.DataArray(py)
                px_xr = xr.concat([px_xr, px_xr_next], dim='dim_0')
                py_xr = xr.concat([py_xr, py_xr_next], dim='dim_0')

            del px, py
            y_count += 1


        #composite:
        px_final = px_xr.mean('dim_0')
        px_final = np.expand_dims(px_final, axis=0) #para graficar
        py_final = py_xr.mean('dim_0')
        py_final = np.expand_dims(py_final, axis=0)  # para graficar

        return px_final, py_final
########################################################################################################################

data = xr.open_dataset(data_dir2 + 'sf.nc')
data = data.rename({'streamfunction': 'var'})

data2 = xr.open_dataset(data_dir2 + 'RWS.nc')

data3 = xr.open_dataset(data_dir2 + 'hgt200_mer_d_w.nc')

from matplotlib import colors
cbar_rws = colors.ListedColormap(['#4FCDF4','white', '#F4C847'])
cbar_rws.set_over('#F49927')
cbar_rws.set_under('#3489F4')
cbar_rws.set_bad(color='white')
########################################################################################################################
#
# seasons = ("Full_Season", 'JJA', 'ASO', 'SON')
# min_max_months = [[7,11], [6,8],[8,10],[9,11]]

seasons = ('JJA', 'SON')
min_max_months = [[6,8], [9,11]]

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

text = False

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

########################################################################################################################
#for comp_of in [False, True]:
comp_of = False
c_cases = 0
for c in cases:
    print(c)
    count = 0
    for s in seasons:
        # Las fechas se toman del periodo 1920-2020 basados en el DMI y N34 con ERSSTv5
        # Cuando se toman los periodos 1920-1949 y 1950_2020 las fechas que no pertencen
        # se excluyen de los composites en CompositeSimple()
        aux = xr.open_dataset(nc_date_dir + '1920_2020' + '_' + s + '.nc')

        neutro = aux.Neutral
        case = aux[c]
        aux.close()
        mmonth = min_max_months[count]

        mmin = mmonth[0]
        mmax = mmonth[-1]

        neutro_comp = CompositeSimple(original_data=data, index=neutro,
                                      mmin=mmin, mmax=mmax)
        data_comp = CompositeSimple(original_data=data, index=case,
                                    mmin=mmin, mmax=mmax)

        comp = data_comp - neutro_comp

        comp = xr.Dataset(
            data_vars=dict(
                var=(['y', 'x'], comp['var'][0, :, :].values)

            ),
            coords=dict(
                lon=(['x'], comp.lon),
                lat=(['y'], comp.lat),
            )
        )

        if comp_of:
            px, py = CompositeOfWAF(original_data=data, index=case,
                                    mmin=mmin, mmax=mmax,
                                    neutro_comp=neutro_comp)
            waf_scale = 1 / 100
        else:
            px, py = WAF(neutro_comp, comp, data_comp.lon, data_comp.lat, reshape=True, variable='var')
            waf_scale = None

        neutro_comp = CompositeSimple(original_data=data2, index=neutro,
                                      mmin=mmin, mmax=mmax)
        data_comp = CompositeSimple(original_data=data2, index=case,
                                    mmin=mmin, mmax=mmax)
        comp = data_comp - neutro_comp
        comp = xr.Dataset(
            data_vars=dict(
                var=(['y', 'x'], comp['var'][0, :, :].values)

            ),
            coords=dict(
                lon=(['x'], comp.lon),
                lat=(['y'], comp.lat),
            )
        )

        neutro_comp2 = CompositeSimple(original_data=data3, index=neutro,
                                       mmin=mmin, mmax=mmax)
        data_comp2 = CompositeSimple(original_data=data3, index=case,
                                     mmin=mmin, mmax=mmax)
        comp2 = data_comp2 - neutro_comp2

        import numpy as np

        if comp_of:
            title = 'RWS [somb], WAF [->] y HGT200 [cont] ' + '\n' + \
                    title_case[c_cases] + ' ' + s + '\n' + ' Composite of WAF'

            name_fig = out_dir + 'Comp_of_WAF_' + c + '_' + s
        else:
            title = 'RWS [somb], WAF [->] y HGT200 [cont] ' + '\n' + \
                    title_case[c_cases] + ' ' + s + '\n' + ' WAF of Composite'

            name_fig = out_dir + 'WAF_' + c + '_' + s + '_NSA'

        case = len(case[c].where(case[c] > 1949, drop=True).values)
        weights = np.transpose(np.tile(-2 * np.cos(np.arange(-90, 89) * 1 * np.pi / 180) + 2.1, (359, 1)))
        weights_arr = np.zeros_like(px)
        weights_arr[0, :, :] = weights

        PlotWAFCountours(comp, comp['var'], title=title, name_fig=name_fig,
                         save=save, dpi=200, levels=np.linspace(-1.2e-10, 1.2e-10, 4),
                         contour=False, cmap=cbar_rws, number_events=case,
                         waf=True, px=px*weights_arr , py=py*weights_arr , text=False, waf_scale=waf_scale,
                         two_variables=True, comp2=comp2, step=1, step_waf=5,
                         levels2=np.linspace(-450, 450, 13), contour0=False)

        count += 1
    c_cases += 1
########################################################################################################################