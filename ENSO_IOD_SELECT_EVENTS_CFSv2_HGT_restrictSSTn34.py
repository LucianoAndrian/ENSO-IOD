"""
A partir de LoadLeads_CFSv2.py y ENSO_IOD_DMI_N34_CFSv2.py
Selecciona para cada miembro de ensambe y lead los campos en los que se dan los eventos
"""
########################################################################################################################
import xarray as xr
import numpy as np

########################################################################################################################
dates_dir = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/DMI_N34_Leads_r/'
dir_leads = '/datos/luciano.andrian/ncfiles/NMME_CFSv2/'
out_dir = '/pikachu/datos/luciano.andrian/cases/cases_restrict_SST/'
########################################################################################################################
def xrClassifierEvents(index, r=None, by_r=True):
    if by_r:
        index_r = index.sel(r=r)
        aux_index_r = index_r.time[np.where(~np.isnan(index_r.index))]
        index_r_f = index_r.sel(time=index_r.time.isin(aux_index_r))

        index_pos = index_r_f.index.time[index_r_f.index > 0]
        index_neg = index_r_f.index.time[index_r_f.index < 0]

        return index_pos, index_neg, index_r_f
    else:
        # print('by_r = False')
        index_pos = index.index.time[index.index > 0]
        index_neg = index.index.time[index.index < 0]
        return index_pos, index_neg


def ConcatEvent(xr_original, xr_to_concat, dim='time'):
    if (len(xr_to_concat.time) != 0) and (len(xr_original.time) != 0):
        xr_concat = xr.concat([xr_original, xr_to_concat], dim=dim)
    elif (len(xr_to_concat.time) == 0) and (len(xr_original.time) != 0):
        xr_concat = xr_original
    elif (len(xr_to_concat.time) != 0) and (len(xr_original.time) == 0):
        xr_concat = xr_to_concat
    elif (len(xr_to_concat.time) == 0) and (len(xr_original.time) == 0):
        return xr_original

    return xr_concat

# def xrFieldTimeDetrend_HGT(xrda, mes, dim, deg=1):
#     # detrend along a single dimension
#     # xrda = xrda.mean('r')
#     # aux = xrda.drop('L')
#     # xrdaL = aux.expand_dims({'L': xrda.L.values})
#     xrdaL = xrda
#
#     for l in [0, 1, 2, 3]:
#         xrda2 = xrdaL.sel(time=xrdaL.time.dt.month.isin(mes - l), L=xrdaL.L.isin(l)).mean('r')
#         aux = xrda2.polyfit(dim=dim, deg=deg)
#         if l == 0:
#             trend = aux
#         else:
#             trend = xr.concat([trend, aux], dim='L')
#
#     trend_f = xr.polyval(xrdaL[dim], trend.hgt_polyfit_coefficients[0]) #+ aux.hgt_polyfit_coefficients[1]
#
#     return trend_f
########################################################################################################################
mmonth_seasons_names = [0,1,2,3] #son las de abajo...
seasons = ['JJA', 'JAS', 'ASO', 'SON']
mmonth_seasons = [7, 8, 9, 10]
sets = [[0,1,2,3]]

for m_name in mmonth_seasons_names:
    print(seasons[m_name])
    data_dmi_s = xr.open_dataset(dates_dir + seasons[m_name] + '_DMI_Leads_r_CFSv2.nc')
    data_dmi_s = data_dmi_s.rename({'sst':'index'})
    data_n34_s = xr.open_dataset(dates_dir + seasons[m_name] + '_N34_Leads_r_CFSv2.nc')
    data_n34_s = data_n34_s.rename({'sst':'index'})
    data_hgt_s = xr.open_dataset(dir_leads + seasons[m_name] + '_hgt_Leads_r_CFSv2.nc')

    data_hgt_s['L'] = [0,1,2,3]
    # aux0 = data_hgt_s
    # aux0 = aux0.drop('L')
    # aux0 = aux0.expand_dims({'L':data_hgt_s.L.values})
    #plt.plot(data_hgt_s.mean(['r', 'lat', 'lon']).hgt);
    #plt.show()

    #Filtrando tendencia para cada Lead
    for l in [0,1,2,3]:
        for ms in [mmonth_seasons[m_name]-l]:
            aux = data_hgt_s.sel(time=data_hgt_s.time.dt.month.isin(ms), L=data_hgt_s.L.isin(l))
            aux = aux.drop('L')

            # #Calculando la tendencia de la media del ensamble (para cada lead)
            aux_mean = aux.mean('r')
            aux1 = aux_mean.polyfit(dim='time', deg=1, full=True)
            trend = xr.polyval(coord=aux.time, coeffs=aux1.hgt_polyfit_coefficients)
            aux3 = aux - trend
            if l == 0:
                aux2 = aux3
            else:
                aux2 = xr.concat([aux2, aux3], dim='L')

    data_hgt_s = aux2

    for s in sets:
        print('Set: ' + str(s))
        l = np.arange(len(s))
        ms = mmonth_seasons[m_name]
        data_dmi = data_dmi_s.sel(time=data_dmi_s.time.dt.month.isin(ms - l))
        data_n34 = data_n34_s.sel(time=data_n34_s.time.dt.month.isin(ms - l))

        data_hgt = data_hgt_s.sel(time=data_hgt_s.time.dt.month.isin(ms - l))

        data_dmi = data_dmi_s.where(np.abs(data_dmi) > 0.75 * data_dmi.std(['time', 'r']))

        data_n34_restricted = data_n34_s.where(np.abs(data_n34) < 0.5*data_n34.std(['time','r']))
        data_n34 = data_n34_s.where(np.abs(data_n34) > data_n34.std(['time', 'r']))

        r_count = 0
        sim_DMIpos_N34neg = -1
        sim_DMIneg_N34pos = -1
        for r in range(1, 25):
            DMI_sim_pos_N34_neg = []
            DMI_sim_neg_N34_pos = []
            DMI_pos, DMI_neg, DMI = xrClassifierEvents(data_dmi.drop('L'), r)
            N34_pos, N34_neg, N34 = xrClassifierEvents(data_n34.drop('L'), r)

            N34_rest_pos, N34_rest_neg, N34_rest = xrClassifierEvents(data_n34_restricted.drop('L'), r)

            # Simultaneous events
            sim_events = np.intersect1d(N34.time, DMI.time)
            DMI_sim = DMI.sel(time=DMI.time.isin(sim_events))
            N34_sim = N34.sel(time=N34.time.isin(sim_events))

            DMI_sim_pos, DMI_sim_neg = xrClassifierEvents(DMI_sim, by_r=False)
            N34_sim_pos, N34_sim_neg = xrClassifierEvents(N34_sim, by_r=False)

            if len(DMI_sim_neg_N34_pos) != 0:
                sim_DMIneg_N34pos += 1
                DMI_sim_neg = DMI_sim_neg[np.in1d(DMI_sim_neg, DMI_sim_neg_N34_pos, invert=True)]
            if len(DMI_sim_pos_N34_neg) != 0:
                sim_DMIpos_N34neg += 1
                DMI_sim_pos = DMI_sim_pos[np.in1d(DMI_sim_pos, DMI_sim_pos_N34_neg, invert=True)]

            # Unique events
            DMI_un = DMI.sel(time=~DMI.time.isin(sim_events))
            N34_un = N34.sel(time=~N34.time.isin(sim_events))

            # Restringed N34 to DMI
            DMI_un = DMI_un.sel(time=DMI_un.time.isin(N34_rest.time))

            DMI_un_pos, DMI_un_neg = xrClassifierEvents(DMI_un, by_r=False)
            N34_un_pos, N34_un_neg = xrClassifierEvents(N34_un, by_r=False)

            """
            los xr seleccionados van a quedar con campos de NaN que corresponden
            a las fechas de otros leads.
            Esto no es un problema para los códigos de graficado ya que antes tambien existian 
            algunos pocos campos de NaN por errores en los archivos de CFSv2
            """
            aux_hgt = data_hgt.sel(r=r)

            data_dmi_pos_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_pos))
            data_dmi_neg_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_neg))

            data_dmi_un_pos_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_un_pos))
            data_dmi_un_neg_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_un_neg))

            data_n34_pos_hgt = aux_hgt.sel(time=aux_hgt.time.isin(N34_pos))
            data_n34_neg_hgt = aux_hgt.sel(time=aux_hgt.time.isin(N34_neg))

            data_n34_un_pos_hgt = aux_hgt.sel(time=aux_hgt.time.isin(N34_un_pos))
            data_n34_un_neg_hgt = aux_hgt.sel(time=aux_hgt.time.isin(N34_un_neg))

            data_sim_pos_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_sim_pos))
            data_sim_neg_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_sim_neg))

            if len(DMI_sim_pos_N34_neg) != 0:
                data_sim_DMIpos_N34neg_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_sim_pos_N34_neg))

            if len(DMI_sim_neg_N34_pos) != 0:
                data_sim_DMIneg_N34pos_hgt = aux_hgt.sel(time=aux_hgt.time.isin(DMI_sim_neg_N34_pos))

            dates_ref = aux_hgt.time
            mask = np.in1d(dates_ref, DMI.time, invert=True)
            neutro = aux_hgt.sel(time=aux_hgt.time.isin(dates_ref[mask]))
            mask = np.in1d(dates_ref, N34_rest.time, invert=True)
            neutro_hgt = neutro.sel(time=neutro.time.isin(dates_ref[mask]))
            del neutro

            if r_count == 0:
                data_dmi_pos_f_hgt = data_dmi_pos_hgt
                data_dmi_neg_f_hgt = data_dmi_neg_hgt
                data_dmi_un_pos_f_hgt = data_dmi_un_pos_hgt
                data_dmi_un_neg_f_hgt = data_dmi_un_neg_hgt
                data_n34_pos_f_hgt = data_n34_pos_hgt
                data_n34_neg_f_hgt = data_n34_neg_hgt
                data_n34_un_pos_f_hgt = data_n34_un_pos_hgt
                data_n34_un_neg_f_hgt = data_n34_un_neg_hgt
                data_sim_pos_f_hgt = data_sim_pos_hgt
                data_sim_neg_f_hgt = data_sim_neg_hgt
                r_count = 1

                neutro_hgt_f = neutro_hgt


            else:
                data_dmi_pos_f_hgt = ConcatEvent(data_dmi_pos_f_hgt, data_dmi_pos_hgt)
                data_dmi_neg_f_hgt = ConcatEvent(data_dmi_neg_f_hgt, data_dmi_neg_hgt)
                data_dmi_un_pos_f_hgt = ConcatEvent(data_dmi_un_pos_f_hgt, data_dmi_un_pos_hgt)
                data_dmi_un_neg_f_hgt = ConcatEvent(data_dmi_un_neg_f_hgt, data_dmi_un_neg_hgt)
                data_n34_pos_f_hgt = ConcatEvent(data_n34_pos_f_hgt, data_n34_pos_hgt)
                data_n34_neg_f_hgt = ConcatEvent(data_n34_neg_f_hgt, data_n34_neg_hgt)
                data_n34_un_pos_f_hgt = ConcatEvent(data_n34_un_pos_f_hgt, data_n34_un_pos_hgt)
                data_n34_un_neg_f_hgt = ConcatEvent(data_n34_un_neg_f_hgt, data_n34_un_neg_hgt)
                data_sim_pos_f_hgt = ConcatEvent(data_sim_pos_f_hgt, data_sim_pos_hgt)
                data_sim_neg_f_hgt = ConcatEvent(data_sim_neg_f_hgt, data_sim_neg_hgt)

                neutro_hgt_f = ConcatEvent(neutro_hgt_f, neutro_hgt)

            if (len(DMI_sim_pos_N34_neg) != 0) and (sim_DMIpos_N34neg == 0):
                data_sim_DMIpos_N34neg_f_hgt = data_sim_DMIpos_N34neg_hgt
                sim_DMIpos_N34neg_anterior = 0
            elif (len(DMI_sim_pos_N34_neg) != 0):
                data_sim_DMIpos_N34neg_f_hgt = ConcatEvent(data_sim_DMIpos_N34neg_f_hgt, data_sim_DMIpos_N34neg_hgt)

            if (len(DMI_sim_neg_N34_pos) != 0) and (sim_DMIneg_N34pos == 0):
                data_sim_DMIneg_N34pos_f_hgt = data_sim_DMIneg_N34pos_hgt
                sim_DMIneg_N34pos_anterior = 0
            elif (len(DMI_sim_neg_N34_pos) != 0):
                data_sim_DMIneg_N34pos_f_hgt = ConcatEvent(data_sim_DMIneg_N34pos_f_hgt, data_sim_DMIneg_N34pos_hgt)

        print('saving...')
        data_dmi_pos_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_DMI_pos.nc')
        data_dmi_neg_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_DMI_neg.nc')
        data_dmi_un_pos_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_DMI_un_pos.nc')
        data_dmi_un_neg_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_DMI_un_neg.nc')
        data_n34_pos_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_N34_pos.nc')
        data_n34_neg_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_N34_neg.nc')
        data_n34_un_pos_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_N34_un_pos.nc')
        data_n34_un_neg_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_N34_un_neg.nc')
        data_sim_pos_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_sim_pos.nc')
        data_sim_neg_f_hgt.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_sim_neg.nc')

        if len(DMI_sim_pos_N34_neg) != 0:
            data_sim_DMIpos_N34neg_f_hgt.to_netcdf(
                out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_sim_DMIpos_N34neg.nc')

        if len(DMI_sim_neg_N34_pos) != 0:
            data_sim_DMIneg_N34pos_f_hgt.to_netcdf(
                out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_sim_DMIneg_N34pos.nc')

        neutro_hgt_f.to_netcdf(out_dir + 'hgt_' + seasons[m_name] + '_Set' + str(s[-1]) + '_NEUTRO.nc')
