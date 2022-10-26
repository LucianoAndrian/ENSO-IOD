
cbar_r = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                                'white',
                                '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'])
cbar_r.set_under('#9B1C00')
cbar_r.set_over('#014A9B')
cbar_r.set_bad(color='white')

cbar = colors.ListedColormap(['#B9391B', '#CD4838', '#E25E55', '#F28C89', '#FFCECC',
                              'white',
                              '#B3DBFF', '#83B9EB', '#5E9AD7', '#3C7DC3', '#2064AF'][::-1])
cbar.set_over('#9B1C00')
cbar.set_under('#014A9B')
cbar.set_bad(color='white')


variables = ['hgt200', 'pp_20CR-V3', 't_20CR-V3']#, 'psl', 'sf', 'div', 'vp']
var_name = ['z','pp', 't', 'psl', 'streamfunction', 'divergence','velocity_potential']
title_var = ['HGT', 'PP', 'Temp', 'PSL', 'Psi', 'Divergence', 'Potential Velocity']

# cmap_revert = [False, True, False, False, True, False, False]
# rescale = [False, False, False, True, False, False, False]
# rescale_value = [1,1,1,1/100,1,1,1]
# wafs = [False, False, False, False, False, False, False]
# two_variables = [False, True, True, False, False, True, False]
# contours = [True, True, True, True, True, False, True]
# SA = [False, True, True, False, False, False, False]
# step = [1, 1, 1, 1, 1, 6, 1]
# contour0 = [True,True,True,False,False,False,False]

two_variables = True
contour = True
SA = True
step = 1
contour0 = False

cmap = ['BrBG', cbar]
# scales = [np.linspace(-450, 450, 21),#hgt
#           np.linspace(-30, 30, 13), #pp
#           np.linspace(-1, 1, 11),#t
#           np.linspace(-3, 3, 13), #psl
#           np.linspace(-4.5e6, 4.5e6, 13),#sf
#           np.linspace(-0.45e-5, 0.45e-5,13), #div
#           np.linspace(-4.5e6, 4.5e6, 13)]#vp

scales_pp_t = [np.linspace(30, 30, 13), # pp
               np.linspace(-1, 1 ,17)] # t
scales = np.linspace(-450, 450, 15)



if bwa:
    bwa_title = 'wo. Pacific Effect - ' + str(start) + ' - 2015'
    bwa_fig = 'woPc_' + str(start)
else:
    bwa_title = 'w. Pacific Effect - ' + str(start) + ' - 2015'
    bwa_fig = 'wPc_' + str(start)


save = True
full_season=True


for v in range(len(variables)): # Terminar!

    if v == 0:
        for v_0 in [0,1,2]:

            if v_0 == 0:
                var = xr.open_dataset(pwd + variables[v_0] + '.nc')
                var = var.sel(time=slice(str(start) + '-01-01', str(end) + '-12-01'))
                DMI_sim, DMI_un, N34_un, dmi_comp, N34_comp, neutral, DMI_sim_pos, \
                DMI_sim_neg, DMI_un_pos, DMI_un_neg, N34_un_pos, N34_un_neg, DMI_pos, \
                DMI_neg, N34_pos, N34_neg = MultipleComposite(var, n34, dmi, season, start=1920, full_season=full_season)

                DMI_sim2 = None
                DMI_un2 = None
                N34_un2 = None
                dmi_comp2 = None
                N34_comp2 = None
                neutral2 = None

            if v_0 == 1:
                var2 = xr.open_dataset(pwd + variables[v_0] + '.nc')
                var2 = var2.sel(time=slice(str(start) + '-01-01', str(end) + '-12-01'))
                var2 = var2.rename({'prate':'var'})
                var2 = var2.__mul__(86400*(365/12))

                DMI_sim2, DMI_un2, N34_un2, dmi_comp2, N34_comp2, neutral2, DMI_sim_pos2, \
                DMI_sim_neg2, DMI_un_pos2, DMI_un_neg2, N34_un_pos2, N34_un_neg2, DMI_pos2, \
                DMI_neg2, N34_pos2, N34_neg2 = MultipleComposite(var2, n34, dmi, season, start=1920, full_season=full_season)
            elif v_0 == 2:
                var2 = xr.open_dataset(pwd + variables[v_0] + '.nc')
                var2 = var2.sel(time=slice(str(start) + '-01-01', str(end) + '-12-01'))
                var2 = var2.rename({'air':'var'})

                DMI_sim2, DMI_un2, N34_un2, dmi_comp2, N34_comp2, neutral2, DMI_sim_pos2, \
                DMI_sim_neg2, DMI_un_pos2, DMI_un_neg2, N34_un_pos2, N34_un_neg2, DMI_pos2, \
                DMI_neg2, N34_pos2, N34_neg2 = MultipleComposite(var2, n34, dmi, season, start=1920, full_season=full_season)

            print("open: " + variables[v_0] + '.nc')


            # Sim
            Plots(data=DMI_sim2, neutral=neutral2, variable='var',
                  data2=DMI_sim, neutral2=neutral,
                  DMI_pos=DMI_sim_pos, DMI_neg=DMI_sim_neg,
                  N34_pos=DMI_sim_pos, N34_neg=DMI_sim_neg,
                  mode='Simultaneus IODs-ENSOs',
                  title=title_var[v_0] + ' - ' + bwa_title + ' - ',
                  neutral_name='All_Neutral',
                  levels=scales[v_0],
                  name_fig=out_dir + 'SIM_' + bwa_fig + '_' + variables[v_0],
                  save=save, contour=False, waf=wafs[v_0],
                  two_variables=two_variables[v_0], levels2=scales[v_0-v_0],
                  season='JJASON', SA=SA[v_0],
                  step=step[v_0], contour0=contour0[v_0],
                  cmap=colors[v_0])

            #un
            Plots(data=DMI_un2, neutral=neutral2, variable='var',
                  data2=DMI_un, neutral2=neutral,
                  DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
                  N34_pos=N34_un_pos, N34_neg=N34_un_neg,
                  mode='Isolated IODs',
                  title=title_var[v_0] + ' - ' + bwa_title + ' - ',
                  neutral_name='All_Neutral',
                  levels=scales[v_0],
                  name_fig=out_dir + 'UN_DMI' + bwa_fig + '_' + variables[v_0],
                  save=save, contour=False, waf=wafs[v_0],
                  two_variables=two_variables[v_0], levels2=scales[v_0 - v_0],
                  season='JJASON', SA=SA[v_0],
                  step=step[v_0], contour0=contour0[v_0],
                  cmap=colors[v_0])

            Plots(data=N34_un2, neutral=neutral2, variable='var',
                  data2=N34_un, neutral2=neutral,
                  DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
                  N34_pos=N34_un_pos, N34_neg=N34_un_neg,
                  mode='Isolated ENSOs',
                  title=title_var[v_0] + ' - ' + bwa_title + ' - ',
                  neutral_name='All_Neutral',
                  levels=scales[v_0],
                  name_fig=out_dir + 'UN_N34' + bwa_fig + '_' + variables[v_0],
                  save=save, contour=False, waf=wafs[v_0],
                  two_variables=two_variables[v_0], levels2=scales[v_0 - v_0],
                  season='JJASON', SA=SA[v_0],
                  step=step[v_0], contour0=contour0[v_0],
                  cmap=colors[v_0],text=False)

            Plots(data=dmi_comp2, neutral=neutral2, variable='var',
                  data2=dmi_comp, neutral2=neutral,
                  DMI_pos=DMI_pos, DMI_neg=DMI_neg,
                  N34_pos=N34_pos, N34_neg=N34_neg,
                  mode='All IODs',
                  title=title_var[v_0] + ' - ' + bwa_title + ' - ',
                  neutral_name='All_Neutral',
                  levels=scales[v_0],
                  name_fig=out_dir + 'All_IOD' + bwa_fig + '_' + variables[v_0],
                  save=save, contour=False, waf=wafs[v_0],
                  two_variables=two_variables[v_0], levels2=scales[v_0 - v_0],
                  season='JJASON', SA=SA[v_0],
                  step=step[v_0], contour0=contour0[v_0],
                  cmap=colors[v_0],text=False)

            Plots(data=N34_comp2, neutral=neutral2, variable='var',
                  data2=N34_comp, neutral2=neutral,
                  DMI_pos=DMI_pos, DMI_neg=DMI_neg,
                  N34_pos=N34_pos, N34_neg=N34_neg,
                  mode='All ENSOs',
                  title=title_var[v_0] + ' - ' + bwa_title + ' - ',
                  neutral_name='All_Neutral',
                  levels=scales[v_0],
                  name_fig=out_dir + 'All_ENSO' + bwa_fig + '_' + variables[v_0],
                  save=save, contour=False, waf=wafs[v_0],
                  two_variables=two_variables[v_0], levels2=scales[v_0 - v_0],
                  season='JJASON', SA=SA[v_0],
                  step=step[v_0], contour0=contour0[v_0],
                  cmap=colors[v_0],text=False)

        v = v_0 + 1
    elif v == 5:
        v_0 = v
        var = xr.open_dataset(pwd + variables[v_0] + '.nc')
        var = var.sel(time=slice(str(start) + '-01-01', str(end) + '-12-01'))
        DMI_sim, DMI_un, N34_un, dmi_comp, N34_comp, neutral, DMI_sim_pos, \
        DMI_sim_neg, DMI_un_pos, DMI_un_neg, N34_un_pos, N34_un_neg, DMI_pos, \
        DMI_neg, N34_pos, N34_neg = MultipleComposite(var, n34, dmi, season, start=1920, full_season=full_season)

        var2 = xr.open_dataset(pwd + variables[v_0+1] + '.nc')
        var2 = var2.sel(time=slice(str(start) + '-01-01', str(end) + '-12-01'))

        DMI_sim2, DMI_un2, N34_un2, dmi_comp2, N34_comp2, neutral2, DMI_sim_pos2, \
        DMI_sim_neg2, DMI_un_pos2, DMI_un_neg2, N34_un_pos2, N34_un_neg2, DMI_pos2, \
        DMI_neg2, N34_pos2, N34_neg2 = MultipleComposite(var2, n34, dmi, season, start=1920, full_season=full_season)

        # Sim
        Plots(data=DMI_sim, neutral=neutral, variable='var',
              data2=DMI_sim2, neutral2=neutral2,
              DMI_pos=DMI_sim_pos, DMI_neg=DMI_sim_neg,
              N34_pos=DMI_sim_pos, N34_neg=DMI_sim_neg,
              mode='Simultaneus IODs-ENSOs',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'SIM_' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season='ASO', SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

        # un
        Plots(data=DMI_un, neutral=neutral, variable='var',
              data2=DMI_un2, neutral2=neutral2,
              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
              mode='Isolated IODs',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'UN_DMI' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

        Plots(data=N34_un, neutral=neutral, variable='var',
              data2=N34_un2, neutral2=neutral2,
              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
              mode='Isolated ENSO',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'UN_N34' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])


        Plots(data=dmi_comp, neutral=neutral, variable='var',
              data2=dmi_comp2, neutral2=neutral2,
              DMI_pos=DMI_pos, DMI_neg=DMI_neg,
              N34_pos=N34_pos, N34_neg=N34_neg,
              mode='All IODs',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'All_IOD' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 +1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

        Plots(data=N34_comp, neutral=neutral, variable='var',
              data2=N34_comp2, neutral2=neutral2,
              DMI_pos=DMI_pos, DMI_neg=DMI_neg,
              N34_pos=N34_pos, N34_neg=N34_neg,
              mode='All ENSOs',
              title=title_var[v_0 + 1] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'All_ENSO' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

    else:
        var = xr.open_dataset(pwd + variables[v] + '.nc')
        var = var.sel(time=slice(str(start) + '-01-01', str(end) + '-12-01'))
        print('open ' + variables[v])
        if rescale[v]:
            print('Rescale')
            var = var.__mul__(rescale_value[v])
        DMI_sim, DMI_un, N34_un, dmi_comp, N34_comp, neutral, DMI_sim_pos, \
        DMI_sim_neg, DMI_un_pos, DMI_un_neg, N34_un_pos, N34_un_neg, DMI_pos, \
        DMI_neg, N34_pos, N34_neg = MultipleComposite(var, n34, dmi, season, start=1920, full_season=full_season)

        DMI_sim2 = None
        DMI_un2 = None
        N34_un2 = None
        dmi_comp2 = None
        N34_comp2 = None
        neutral2 = None

        v_0 = v



        # Sim
        Plots(data=DMI_sim, neutral=neutral, variable='var',
              data2=DMI_sim2, neutral2=neutral2,
              DMI_pos=DMI_sim_pos, DMI_neg=DMI_sim_neg,
              N34_pos=DMI_sim_pos, N34_neg=DMI_sim_neg,
              mode='Simultaneus IODs-ENSOs',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'SIM_' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season='ASO', SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

        # un
        Plots(data=DMI_un, neutral=neutral, variable='var',
              data2=DMI_un2, neutral2=neutral2,
              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
              mode='Isolated IODs',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'UN_DMI' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

        Plots(data=N34_un, neutral=neutral, variable='var',
              data2=N34_un2, neutral2=neutral2,
              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
              mode='Isolated ENSOs',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'UN_DMI' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

        Plots(data=dmi_comp, neutral=neutral, variable='var',
              data2=dmi_comp2, neutral2=neutral2,
              DMI_pos=DMI_pos, DMI_neg=DMI_neg,
              N34_pos=N34_pos, N34_neg=N34_neg,
              mode='All IODs',
              title=title_var[v_0] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'All_IOD' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 +1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])

        Plots(data=N34_comp, neutral=neutral, variable='var',
              data2=N34_comp2, neutral2=neutral2,
              DMI_pos=DMI_pos, DMI_neg=DMI_neg,
              N34_pos=N34_pos, N34_neg=N34_neg,
              mode='All ENSOs',
              title=title_var[v_0 + 1] + ' - ' + bwa_title + ' - ',
              neutral_name='All_Neutral',
              levels=scales[v_0],
              name_fig=out_dir + 'All_ENSO' + bwa_fig + '_' + variables[v_0],
              save=save, contour=False, waf=wafs[v_0],
              two_variables=two_variables[v_0], levels2=scales[v_0 + 1],
              season=season, SA=SA[v_0],
              step=step[v_0], contour0=contour0[v_0],
              cmap=colors[v_0])



