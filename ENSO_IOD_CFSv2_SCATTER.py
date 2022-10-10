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
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/Modelos/Scatter/'

neutro_dir = '/pikachu/datos/luciano.andrian/cases/cases_restrict_SST/'
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
        #print('by_r = False')
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
########################################################################################################################
mmonth_seasons_names = [0,1,2,3] #son las de abajo...
seasons = ['JJA', 'JAS', 'ASO', 'SON']
mmonth_seasons = [7, 8, 9, 10]
sets = [[0,1,2,3]]

variables = ['n34', 'dmi']

for m_name in mmonth_seasons_names:
    print(seasons[m_name])
    data_dmi_s = xr.open_dataset(dates_dir + seasons[m_name] + '_DMI_Leads_r_CFSv2.nc')
    data_n34_s = xr.open_dataset(dates_dir + seasons[m_name] + '_N34_Leads_r_CFSv2.nc')

    for s in sets:
        print('Set: ' + str(s))
        l = np.arange(len(s))
        ms = mmonth_seasons[m_name]
        data_dmi = data_dmi_s.sel(time=data_dmi_s.time.dt.month.isin(ms-l))
        data_n34 = data_n34_s.sel(time=data_n34_s.time.dt.month.isin(ms-l))
        aux_n342 = data_n34.copy()
        aux_dmi2 = data_dmi.copy()

        dmi_sd = data_dmi.std(['time','r'])
        n34_sd = data_n34.std(['time','r'])

        data_dmi_no_iod = data_dmi_s.where(np.abs(data_dmi) < 0.75*dmi_sd)
        data_n34_no_enso = data_n34_s.where(np.abs(data_n34) < n34_sd)

        data_n34_no_enso2 = data_n34_s.where(np.abs(data_n34) < 0.5*n34_sd)

        data_n34_restricted = data_n34_s.where(np.abs(data_n34) < 0.5*n34_sd)

        data_dmi = data_dmi_s.where(np.abs(data_dmi) > 0.75*dmi_sd)
        data_n34 = data_n34_s.where(np.abs(data_n34) > n34_sd)

        data_dmi = data_dmi.__mul__(1 / dmi_sd)
        aux_dmi2 = aux_dmi2.__mul__(1 / dmi_sd)
        data_n34 = data_n34.__mul__(1 / n34_sd)
        aux_n342 = aux_n342.__mul__(1 / n34_sd)

        data_dmi_no_iod = data_dmi_no_iod.__mul__(1 / dmi_sd)
        data_n34_no_enso = data_n34_no_enso.__mul__(1 / n34_sd)
        data_n34_no_enso2 = data_n34_no_enso2.__mul__(1 / n34_sd)

        r_count = 0
        sim_DMIpos_N34neg=-1
        sim_DMIneg_N34pos=-1
        data_sim_DMIneg_N34pos_f_dmi = []
        data_sim_DMIpos_N34neg_f_dmi = []
        for r in range(1, 25):
            DMI_sim_pos_N34_neg = []
            DMI_sim_neg_N34_pos = []
            DMI_pos, DMI_neg, DMI = xrClassifierEvents(data_dmi.drop('L'), r)
            N34_pos, N34_neg, N34 = xrClassifierEvents(data_n34.drop('L'), r)
            x, y, aux_n3422 = xrClassifierEvents(aux_n342.drop('L'), r)
            x, y, aux_dmi22 = xrClassifierEvents(aux_dmi2.drop('L'), r)

            x, y, aux_n34_no_enso = xrClassifierEvents(data_n34_no_enso.drop('L'), r)
            x, y, aux_n34_no_enso2 = xrClassifierEvents(data_n34_no_enso2.drop('L'), r)
            x, y, aux_dmi_no_iod = xrClassifierEvents(data_dmi_no_iod.drop('L'), r)

            x,y, N34_rest = xrClassifierEvents(data_n34_restricted.drop('L'), r)



            # Simultaneous events
            sim_events = np.intersect1d(N34.time, DMI.time)
            DMI_sim = DMI.sel(time=DMI.time.isin(sim_events))
            N34_sim = N34.sel(time=N34.time.isin(sim_events))

            DMI_sim_pos, DMI_sim_neg = xrClassifierEvents(DMI_sim, by_r=False)
            N34_sim_pos, N34_sim_neg = xrClassifierEvents(N34_sim, by_r=False)

            DMI_sim_pos_N34_neg = np.intersect1d(DMI_sim_pos, N34_sim_neg)
            DMI_sim_neg_N34_pos = np.intersect1d(DMI_sim_neg, N34_sim_pos)

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

            aux_dmi = data_dmi.sel(r=r)

            data_dmi_un_pos_dmi = aux_dmi.sel(time=aux_dmi.time.isin(DMI_un_pos))
            data_dmi_un_pos_dmi_n34_values = aux_n3422.sel(time=aux_n3422.time.isin(DMI_un_pos))

            data_dmi_un_neg_dmi = aux_dmi.sel(time=aux_dmi.time.isin(DMI_un_neg))
            data_dmi_un_neg_dmi_n34_values = aux_n3422.sel(time=aux_n3422.time.isin(DMI_un_neg))

            data_sim_pos_dmi = aux_dmi.sel(time=aux_dmi.time.isin(DMI_sim_pos))
            data_sim_neg_dmi = aux_dmi.sel(time=aux_dmi.time.isin(DMI_sim_neg))

            aux_n34 = data_n34.sel(r=r)

            data_dmi_un_pos_n34 = aux_n34.sel(time=aux_n34.time.isin(DMI_un_pos))
            data_dmi_un_neg_n34 = aux_n34.sel(time=aux_n34.time.isin(DMI_un_neg))

            data_n34_pos_n34 = aux_n34.sel(time=aux_n34.time.isin(N34_pos))
            data_n34_neg_n34 = aux_n34.sel(time=aux_n34.time.isin(N34_neg))

            data_n34_un_pos_n34 = aux_n34.sel(time=aux_n34.time.isin(N34_un_pos))
            data_n34_un_pos_n34_dmi_values = aux_dmi22.sel(time=aux_dmi22.time.isin(N34_un_pos))
            data_n34_un_neg_n34 = aux_n34.sel(time=aux_n34.time.isin(N34_un_neg))
            data_n34_un_neg_n34_dmi_values = aux_dmi22.sel(time=aux_dmi22.time.isin(N34_un_neg))

            data_sim_pos_n34 = aux_n34.sel(time=aux_n34.time.isin(DMI_sim_pos))
            data_sim_neg_n34 = aux_n34.sel(time=aux_n34.time.isin(DMI_sim_neg))



            print(r)
            print(np.where(data_sim_pos_n34.index.values<0))
            if len(DMI_sim_pos_N34_neg) != 0:
                data_sim_DMIpos_N34neg_n34 = aux_n34.sel(time=aux_n34.time.isin(DMI_sim_pos_N34_neg))
                data_sim_DMIpos_N34neg_dmi = aux_dmi.sel(time=aux_dmi.time.isin(DMI_sim_pos_N34_neg))

            if len(DMI_sim_neg_N34_pos) != 0:
                data_sim_DMIneg_N34pos_n34 = aux_n34.sel(time=aux_n34.time.isin(DMI_sim_neg_N34_pos))
                data_sim_DMIneg_N34pos_dmi = aux_dmi.sel(time=aux_dmi.time.isin(DMI_sim_neg_N34_pos))

            # dates_ref = aux_dmi.time
            # mask = np.in1d(dates_ref, DMI.time, invert=True)
            # neutro = aux_dmi.sel(time=aux_dmi.time.isin(dates_ref[mask]))
            # mask = np.in1d(dates_ref, N34_rest.time, invert=True)
            # neutro_dmi = neutro.sel(time=neutro.time.isin(dates_ref[mask]))
            # del neutro
            #
            # dates_ref = aux_n34.time
            # mask = np.in1d(dates_ref, DMI.time, invert=True)
            # neutro = aux_n34.sel(time=aux_n34.time.isin(dates_ref[mask]))
            # mask = np.in1d(dates_ref, N34_rest.time, invert=True)
            # neutro_n34 = neutro.sel(time=neutro.time.isin(dates_ref[mask]))
            #
            # neutro_dmi_no_n34 = neutro_dmi.sel(
            #     time=neutro_dmi.time.isin(
            #         neutro_dmi.time[np.in1d(neutro_dmi.time, neutro_n34.time, invert=True)]))
            #
            # neutro_n34_no_dmi = neutro_n34.sel(
            #     time=neutro_n34.time.isin(
            #         neutro_n34.time[np.in1d(neutro_n34.time, neutro_dmi.time, invert=True)]))


            aux_n34_no_enso = aux_n34_no_enso.sel(
                time=aux_n34_no_enso.time.isin(
                    aux_n34_no_enso.time[np.in1d(aux_n34_no_enso.time, data_dmi_un_pos_n34.time, invert=True)]))

            aux_n34_no_enso = aux_n34_no_enso.sel(
                time=aux_n34_no_enso.time.isin(
                    aux_n34_no_enso.time[np.in1d(aux_n34_no_enso.time, data_dmi_un_neg_n34.time, invert=True)]))

            aux_dmi_no_iod = aux_dmi_no_iod.sel(
                time=aux_dmi_no_iod.time.isin(
                    aux_dmi_no_iod.time[np.in1d(aux_dmi_no_iod.time, data_n34_un_pos_n34.time, invert=True)]))

            aux_dmi_no_iod = aux_dmi_no_iod.sel(
                time=aux_dmi_no_iod.time.isin(
                    aux_dmi_no_iod.time[np.in1d(aux_dmi_no_iod.time, data_n34_un_neg_n34.time, invert=True)]))

            aux_n34_no_enso2 = aux_n34_no_enso2.sel(
                time=aux_n34_no_enso2.time.isin(
                    aux_n34_no_enso2.time[np.in1d(aux_n34_no_enso2.time, data_dmi_un_pos_n34.time, invert=True)]))

            aux_n34_no_enso2 = aux_n34_no_enso2.sel(
                time=aux_n34_no_enso2.time.isin(
                    aux_n34_no_enso2.time[np.in1d(aux_n34_no_enso2.time, data_dmi_un_neg_n34.time, invert=True)]))


            aux_n34_no_enso_dmi = aux_dmi22.sel(time=aux_dmi22.time.isin(aux_n34_no_enso.time.values))
            aux_dmi_no_iod_n34 = aux_n3422.sel(time=aux_n3422.time.isin(aux_dmi_no_iod.time.values))

            aux_n34_no_enso_dmi2 = aux_dmi22.sel(time=aux_dmi22.time.isin(aux_n34_no_enso2.time.values))


            if r_count == 0:

                data_dmi_un_pos_f_dmi = data_dmi_un_pos_dmi
                data_dmi_un_neg_f_dmi = data_dmi_un_neg_dmi
                data_sim_pos_f_dmi = data_sim_pos_dmi
                data_sim_neg_f_dmi = data_sim_neg_dmi
                # que elegancia la de francia...

                data_dmi_un_pos_f_n34 = data_dmi_un_pos_n34
                data_dmi_un_neg_f_n34 = data_dmi_un_neg_n34
                data_n34_pos_f_n34 = data_n34_pos_n34
                data_n34_neg_f_n34 = data_n34_neg_n34
                data_n34_un_pos_f_n34 = data_n34_un_pos_n34
                data_n34_un_neg_f_n34 = data_n34_un_neg_n34
                data_sim_pos_f_n34 = data_sim_pos_n34
                data_sim_neg_f_n34 = data_sim_neg_n34
                r_count = 1

                # neutro_dmi_f = neutro_dmi
                # neutro_n34_f = neutro_n34

                data_dmi_un_pos_dmi_n34_values_f = data_dmi_un_pos_dmi_n34_values
                data_dmi_un_neg_dmi_n34_values_f = data_dmi_un_neg_dmi_n34_values

                data_n34_un_pos_n34_dmi_values_f = data_n34_un_pos_n34_dmi_values
                data_n34_un_neg_n34_dmi_values_f = data_n34_un_neg_n34_dmi_values

                data_no_enso = aux_n34_no_enso
                data_no_enso2 = aux_n34_no_enso2
                data_no_enso_dmi = aux_n34_no_enso_dmi
                data_no_enso_dmi2 = aux_n34_no_enso_dmi2
                data_no_iod = aux_dmi_no_iod
                data_no_iod_n34 = aux_dmi_no_iod_n34

            else:

                data_dmi_un_pos_f_dmi = ConcatEvent(data_dmi_un_pos_f_dmi, data_dmi_un_pos_dmi)
                data_dmi_un_neg_f_dmi = ConcatEvent(data_dmi_un_neg_f_dmi, data_dmi_un_neg_dmi)
                data_sim_pos_f_dmi = ConcatEvent(data_sim_pos_f_dmi, data_sim_pos_dmi)
                data_sim_neg_f_dmi = ConcatEvent(data_sim_neg_f_dmi, data_sim_neg_dmi)
                data_dmi_un_pos_f_n34 = ConcatEvent(data_dmi_un_pos_f_n34, data_dmi_un_pos_n34)
                data_dmi_un_neg_f_n34 = ConcatEvent(data_dmi_un_neg_f_n34, data_dmi_un_neg_n34)
                data_n34_pos_f_n34 = ConcatEvent(data_n34_pos_f_n34, data_n34_pos_n34)
                data_n34_neg_f_n34 = ConcatEvent(data_n34_neg_f_n34, data_n34_neg_n34)
                data_n34_un_pos_f_n34 = ConcatEvent(data_n34_un_pos_f_n34, data_n34_un_pos_n34)
                data_n34_un_neg_f_n34 = ConcatEvent(data_n34_un_neg_f_n34, data_n34_un_neg_n34)
                data_sim_pos_f_n34 = ConcatEvent(data_sim_pos_f_n34, data_sim_pos_n34)
                data_sim_neg_f_n34 = ConcatEvent(data_sim_neg_f_n34, data_sim_neg_n34)

                # neutro_dmi_f = ConcatEvent(neutro_dmi_f, neutro_dmi)
                # neutro_n34_f = ConcatEvent(neutro_n34_f, neutro_n34)

                data_dmi_un_pos_dmi_n34_values_f = ConcatEvent(data_dmi_un_pos_dmi_n34_values_f, data_dmi_un_pos_dmi_n34_values)
                data_dmi_un_neg_dmi_n34_values_f = ConcatEvent(data_dmi_un_neg_dmi_n34_values_f, data_dmi_un_neg_dmi_n34_values)

                data_n34_un_pos_n34_dmi_values_f = ConcatEvent(data_n34_un_pos_n34_dmi_values_f, data_n34_un_pos_n34_dmi_values)
                data_n34_un_neg_n34_dmi_values_f = ConcatEvent(data_n34_un_neg_n34_dmi_values_f, data_n34_un_neg_n34_dmi_values)

                data_no_enso = ConcatEvent(data_no_enso, aux_n34_no_enso)
                data_no_enso_dmi = ConcatEvent(data_no_enso_dmi, aux_n34_no_enso_dmi)

                data_no_enso2 = ConcatEvent(data_no_enso2, aux_n34_no_enso2)
                data_no_enso_dmi2 = ConcatEvent(data_no_enso_dmi2, aux_n34_no_enso_dmi2)

                data_no_iod = ConcatEvent(data_no_iod, aux_dmi_no_iod)
                data_no_iod_n34 = ConcatEvent(data_no_iod_n34, aux_dmi_no_iod_n34)

            if (len(DMI_sim_pos_N34_neg) != 0) and (sim_DMIpos_N34neg == 0):
                data_sim_DMIpos_N34neg_f_n34 = data_sim_DMIpos_N34neg_n34
                data_sim_DMIpos_N34neg_f_dmi = data_sim_DMIpos_N34neg_dmi
                sim_DMIpos_N34neg_anterior = 0
            elif (len(DMI_sim_pos_N34_neg) !=0):
                data_sim_DMIpos_N34neg_f_n34 = ConcatEvent(data_sim_DMIpos_N34neg_f_n34, data_sim_DMIpos_N34neg_n34)
                data_sim_DMIpos_N34neg_f_dmi = ConcatEvent(data_sim_DMIpos_N34neg_f_dmi, data_sim_DMIpos_N34neg_dmi)

            if (len(DMI_sim_neg_N34_pos) != 0) and (sim_DMIneg_N34pos == 0):
                data_sim_DMIneg_N34pos_f_n34 = data_sim_DMIneg_N34pos_n34
                data_sim_DMIneg_N34pos_f_dmi = data_sim_DMIneg_N34pos_dmi
                sim_DMIneg_N34pos_anterior = 0
            elif (len(DMI_sim_neg_N34_pos) != 0):
                data_sim_DMIneg_N34pos_f_n34 = ConcatEvent(data_sim_DMIneg_N34pos_f_n34, data_sim_DMIneg_N34pos_n34)
                data_sim_DMIneg_N34pos_f_dmi = ConcatEvent(data_sim_DMIneg_N34pos_f_dmi, data_sim_DMIneg_N34pos_dmi)


        if s == sets[0]:
            s_name = '0'
        elif s == sets[1]:
            s_name = '1'
        elif s == sets[2]:
            s_name = '3'
        else:
            s_name = '7'


        save = True
        dpi = 400
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(dpi=dpi)
        #
        # plt.scatter(x=data_no_iod.index.values, y=data_no_iod_n34.index.values, marker='o',
        #             s=20, edgecolor='k', color='gray', alpha=.2)
        plt.scatter(y=data_no_enso.index.values, x=data_no_enso_dmi.index.values, marker='o',
                    s=20, edgecolor='k', color='gray', alpha=.2)

        # plt.scatter(x=data_no_iod.index.values, y=data_no_iod_n34.index.values, marker='o',
        #             s=20, edgecolor='magenta', color='gray', alpha=.5)
        plt.scatter(y=data_no_enso2.index.values, x=data_no_enso_dmi2.index.values, marker='o',
                    s=20, edgecolor='k', color='gray', alpha=.6)
        # plt.scatter(y=data_no_enso2.index.values, x=data_no_enso_dmi2.index.values, marker='x',
        #             s=5, edgecolor='k', color='k', alpha=.2)



        plt.scatter(x=data_dmi_un_pos_f_dmi.index.values, y=data_dmi_un_pos_dmi_n34_values_f.index.values, marker='>',
                    s=20, edgecolor='k', color='firebrick', alpha=.5)
        plt.scatter(x=data_dmi_un_neg_f_dmi.index.values, y=data_dmi_un_neg_dmi_n34_values_f.index.values, marker='<',
                    s=20, edgecolor='k', color='lime', alpha=0.5)

        plt.scatter(y=data_n34_un_pos_f_n34.index.values, x=data_n34_un_pos_n34_dmi_values_f.index.values, marker='^',
                    s=20, edgecolor='k', color='darkorange', alpha=0.5)
        plt.scatter(y=data_n34_un_neg_f_n34.index.values, x=data_n34_un_neg_n34_dmi_values_f.index.values, marker='v',
                    s=20, edgecolor='k', color='blue', alpha=0.5)

        if len(data_sim_DMIpos_N34neg_f_dmi) != 0:
            plt.scatter(x=data_sim_DMIpos_N34neg_f_dmi.index.values, y=data_sim_DMIpos_N34neg_f_n34.index.values,
                        marker='o', s=20, edgecolor='k', color='purple', alpha=0.5)

        if len(data_sim_DMIneg_N34pos_f_dmi) != 0:
            plt.scatter(x=data_sim_DMIneg_N34pos_f_dmi.index.values, y=data_sim_DMIneg_N34pos_f_n34.index.values,
                        marker='8', s=20, edgecolor='k', color='orange', alpha=0.5)

        plt.scatter(x=data_sim_pos_f_dmi.index.values, y=data_sim_pos_f_n34.index.values, marker='h', s=20,
                    edgecolor='k',
                    color='red', alpha=0.5)
        plt.scatter(x=data_sim_neg_f_dmi.index.values, y=data_sim_neg_f_n34.index.values, marker='s', s=20,
                    edgecolor='k',
                    color='lightseagreen', alpha=0.5)

        plt.ylim((-4, 4));
        plt.xlim((-4, 4))
        plt.axhspan(-1, 1, alpha=0.2, color='black', zorder=0)
        plt.axvspan(-.75, .75, alpha=0.2, color='black', zorder=0)
        # ax.grid(True)
        fig.set_size_inches(6, 6)
        plt.xlabel('IOD', size=15)
        plt.ylabel('Niño 3.4', size=15)

        plt.text(-3.8, 3.6, 'EN/IOD-', dict(size=15))
        plt.text(-.3, 3.6, 'EN', dict(size=15))
        plt.text(+2.3, 3.6, 'EN/IOD+', dict(size=15))
        plt.text(+3, -.1, 'IOD+', dict(size=15))
        plt.text(+2.6, -3.9, 'LN/IOD+', dict(size=15))
        plt.text(-.3, -3.9, 'LN', dict(size=15))
        plt.text(-3.8, -3.9, ' LN/IOD-', dict(size=15))
        plt.text(-3.8, -.1, 'IOD-', dict(size=15))
        plt.title(seasons[m_name] + ' - ' + 'Leads: ' + str(s))
        plt.tight_layout()
        if save:
            plt.savefig(out_dir + 'ENSO_IOD_CFSv2_Scatter_restSST' + seasons[m_name] + '_Set' + s_name + '.jpg')
        else:
            plt.show()




