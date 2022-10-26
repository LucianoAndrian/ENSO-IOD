main_month_name = ''
def PlotComp(comp, comp_var, title='Fig', fase=None, name_fig='Fig',
             save=False, dpi=200, levels=np.linspace(-1.5, 1.5, 13),
             contour=False,cmap='RdBu_r', number_events='', main_month_name=main_month_name,
             waf=False, px=None, py=None, text=True, SA=False):

    from numpy import ma
    import matplotlib.pyplot as plt

    if SA:
        fig = plt.figure(figsize=(5, 5), dpi=dpi)
    else:
        fig = plt.figure(figsize=(7, 3.5), dpi=dpi)

    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    if SA:
        ax.set_extent([270,330, -60,20],crs_latlon)
    else:
        ax.set_extent([0, 359, -80, 40], crs=crs_latlon)

    im = ax.contourf(comp.lon, comp.lat, comp_var, levels=levels,
                     transform=crs_latlon, cmap=cmap, extend='both')

    if contour:
        values2 = ax.contour(comp.lon, comp.lat, comp_var, levels=levels,
                            transform=crs_latlon, colors='k', linewidths=0.7)


        values = ax.contour(comp.lon, comp.lat, comp_var, levels=0,
                            transform=crs_latlon, colors='darkgray', linewidths=1)
        ax.clabel(values, inline=1, fontsize=5, fmt='%1.1f')


    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    if SA:
        ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
        ax.set_yticks(np.arange(-60, 40, 20), crs=crs_latlon)
    else:
        ax.set_xticks(np.arange(30, 330, 60), crs=crs_latlon)
        ax.set_yticks(np.arange(-80, 20, 20), crs=crs_latlon)
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)

    if waf:
        Q60 = np.percentile(np.sqrt(np.add(np.power(px, 2), np.power(py, 2))), 0)
        M = np.sqrt(np.add(np.power(px, 2), np.power(py, 2))) < Q60
        # mask array
        px_mask = ma.array(px, mask=M)
        py_mask = ma.array(py, mask=M)
        # plot vectors
        lons, lats = np.meshgrid(comp.lon.values, comp.lat.values)
        ax.quiver(lons[::20, ::20], lats[::20, ::20], px_mask[0, ::20, ::20],
                  py_mask[0, ::20, ::20], transform=crs_latlon,pivot='mid'
                  , scale=1/50)#, width=1.5e-3, headwidth=3.1,  # headwidht (default3)
                  #headlength=2.2)  # (default5))

    plt.title(str(title) + ' - ' + str(main_month_name) + '  ' + str(fase.split(' ', 1)[1]), fontsize=10)
    if text:
        plt.figtext(0.5, 0.01, number_events, ha="center", fontsize=10,
                bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5})
    plt.tight_layout()

    if save:
        plt.savefig(name_fig + str(main_month_name) + '_' + str(fase.split(' ', 1)[1]) + '.jpg')
        plt.close()
    else:
        plt.show()

def PlotComp2var(comp1, comp_var1, comp2, comp_var2, title='Fig', fase=None, name_fig='Fig',
                 save=False, dpi=200, levels=np.linspace(-1.5, 1.5, 13),
                 levels2=np.linspace(-1.5, 1.5, 13), main_month_name=main_month_name,
                 cmap='RdBu_r', number_events='',text=True, SA = False, step=6, contour0=True):

    from numpy import ma
    import matplotlib.pyplot as plt

    if SA:
        fig = plt.figure(figsize=(5, 5), dpi=dpi)
    else:
        fig = plt.figure(figsize=(7, 3.5), dpi=dpi)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    if SA:
        ax.set_extent([270,330, -60,20], crs=crs_latlon)
    else:
        ax.set_extent([0, 359, -80, 40], crs=crs_latlon)


    im = ax.contourf(comp1.lon[::step], comp1.lat[::step], comp_var1[::step,::step], levels=levels,
                     transform=crs_latlon, cmap=cmap, extend='both')
    if contour0:
        ax.contour(comp1.lon, comp1.lat, comp_var1, levels=0,
                   transform=crs_latlon, colors='grey', linewidths=1)


    values = ax.contour(comp2.lon, comp2.lat, comp_var2, levels=levels2,
                        transform=crs_latlon, colors='k', linewidths=1)



    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.COASTLINE)
    # ax.add_feature(cartopy.feature.BORDERS, linestyle='-', alpha=.5)
    ax.gridlines(crs=crs_latlon, linewidth=0.3, linestyle='-')
    if SA:
        ax.set_xticks(np.arange(270, 330, 10), crs=crs_latlon)
        ax.set_yticks(np.arange(-60, 20, 20), crs=crs_latlon)
    else:
        ax.set_xticks(np.arange(30, 330, 60), crs=crs_latlon)
        ax.set_yticks(np.arange(-80, 40, 20), crs=crs_latlon)


    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=7)
    plt.title(str(title) + ' - ' + str(main_month_name) + '  ' + str(fase.split(' ', 1)[1]), fontsize=10)
    if text:
        plt.figtext(0.5, 0.01, number_events, ha="center", fontsize=10,
                    bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5})
    plt.tight_layout()

    if save:
        plt.savefig(name_fig + str(main_month_name) + '_' + str(fase.split(' ', 1)[1]) + '.jpg')
        plt.close()
    else:
        plt.show()

def Plots(data, variable, neutral, DMI_pos, DMI_neg,
          N34_pos, N34_neg, neutral_name="",cmap='RdBu_r',
          dpi = 200, mode = "", levels=np.linspace(-1.5, 1.5, 13),
          name_fig='', save=False, contour=False, title = "",waf=False,
          two_variables=False, data2=None, neutral2=None, levels2=None,
          main_month_name=main_month_name,text=True, SA=False, contour0=False, step=6):

    def Title(DMI_phase, N34_phase, title=title):
        DMI_phase = set(DMI_phase)
        N34_phase = set(N34_phase)
        if mode.split(' ', 1)[0] != 'Simultaneus':
            if mode.split(' ', 1)[1] == 'IODs':
                title = title + mode + ': ' + str(len(DMI_phase)) + '\n' + 'against ' + clim
                number_events = str(DMI_phase)
            else:
                title = title + mode + ': ' + str(len(N34_phase)) + '\n' + 'against ' + clim
                number_events = str(N34_phase)

        elif mode.split(' ', 1)[0] == 'Simultaneus':
            title = title +mode + '\n' + 'IODs: ' + str(len(DMI_phase)) + \
                    ' - ENSOs: ' + str(len(N34_phase)) + '\n' + 'against ' + clim
            number_events = str(N34_phase)
        return title, number_events

    if two_variables:
        if data[0] != 0:
            comp = data[0] - neutral
            clim = neutral_name

            comp2= data2[0] - neutral2

            PlotComp2var(comp1=comp, comp_var1=comp[variable],
                         comp2=comp2, comp_var2=comp2[variable],
                         dpi=dpi, fase='- Positive',
                         title=Title(DMI_phase=DMI_pos, N34_phase=N34_pos)[0],
                         levels=levels, cmap=cmap, levels2 = levels2,
                         save=save, name_fig=name_fig,
                         number_events=Title(DMI_phase=DMI_pos, N34_phase=N34_pos)[1],
                         main_month_name=main_month_name,text=text, SA=SA,
                         contour0=contour0,step=step)

        if data[1] != 0:
            comp = data[1] - neutral
            clim = neutral_name
            comp2 = data2[1] - neutral2

            PlotComp2var(comp1=comp, comp_var1=comp[variable],
                         comp2=comp2, comp_var2=comp2[variable],
                         dpi=dpi, fase='- Negative',
                         title=Title(DMI_phase=DMI_neg, N34_phase=N34_neg)[0],
                         levels=levels, cmap=cmap,  levels2 = levels2,
                         save=save, name_fig=name_fig,
                         number_events=Title(DMI_phase=DMI_neg, N34_phase=N34_neg)[1],
                         main_month_name=main_month_name,text=text,SA=SA,
                         contour0=contour0, step=step)

    else:
        if data[0] != 0:
            comp = data[0] - neutral
            clim = neutral_name

            if waf:
                px, py = WAF(psiclm=neutral, psiaa=comp, lon=comp.lon, lat=comp.lat, reshape=True, variable='var')
            else:
                px = py = None

            PlotComp(comp=comp, comp_var=comp[variable], dpi=dpi,
                     fase='- Positive',
                     title=Title(DMI_phase=DMI_pos, N34_phase=N34_pos)[0],
                     levels=levels, cmap=cmap, save=save, name_fig=name_fig,
                     contour=contour, number_events=Title(DMI_phase=DMI_pos, N34_phase=N34_pos)[1],
                     waf=waf, px=px, py=py,main_month_name=main_month_name,text=text,SA=SA)

        if data[1] != 0:
            comp = data[1] - neutral
            clim = neutral_name

            if waf:
                px, py = WAF(psiclm=neutral, psiaa=comp, lon=comp.lon, lat=comp.lat, reshape=True, variable='var')
            else:
                px = py = None

            PlotComp(comp=comp, comp_var=comp[variable], dpi=dpi,
                     fase='- Negative',
                     title=Title(DMI_phase=DMI_neg, N34_phase=N34_neg)[0],
                     levels=levels, cmap=cmap, save=save, name_fig=name_fig,
                     contour=contour, number_events=Title(DMI_phase=DMI_neg, N34_phase=N34_neg)[1],
                     waf=waf, px=px, py=py,main_month_name=main_month_name,text=text, SA=SA)

def CompositeComputeAndPlot(variables, var_name, scales,
                            cmap_revert, contours, wafs,
                            rescale, rescale_value, start='1920',
                            season_start=8, season_end=10,
                            save=False, two_variables=None,
                            full_season=False,text=True, cmap_data=False,
                            n34=None, dmi=None,bwa=False
                            , SA=None,step=None,contour0=None):


    from matplotlib import colors
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

    seasons = ['DJF', 'JFM', 'FMA', 'MAM', 'AMJ', 'MJJ', 'JJA',
               'JAS', 'ASO', 'SON', 'OND', 'NDJ']
    title_var = ['PSL', 'PP', 'HGT', 'Temp', 'HGT', 'HGT', 'Psi', 'Divergence', 'Potential Velocity']

    pwd = '/datos/luciano.andrian/ncfiles/'
    out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/'

    if bwa:
        bwa_title = 'wo. Pacific Effect - ' + str(start) + ' - 2020'
        bwa_fig = 'woPc_' + str(start)
    else:
        bwa_title = 'w. Pacific Effect - ' + str(start) + ' - 2020'
        bwa_fig = 'wPc_' + str(start)

    if full_season:
        print('#----Full Season JJASON----#')
        print('main_month ignore')
        print('seasons ignore')
        text = False
        season_end = season_start +1

    for v in range(len(variables)):
        if v ==8:
            break
        if v == 2:
            v = 3
        elif v == 4:
            v = 5



        if two_variables[v]:
            print('Two variables')

            if (v == 1) | (v == 3):
                if cmap_data &(v == 1) & (start =='1979'):
                    variables[v] = 'pp_CMAP.nc'

                var = xr.open_dataset('/home/luciano.andrian/doc/scrips/ncfiles/' + variables[v] + '.nc')

                print("open: " + variables[v] + '.nc')
                if v == 1:
                    var = var.sel(time=slice(start + '-01-01', '2020-12-01'))
                    var = var.rename({'precip': 'var'})
                    if cmap_data:
                        var = var.__mul__(365/12)
                else:
                    start = '1950'
                    print('Start 1950 (GHCN-CAMS')
                    var = var.sel(time=slice(start + '-01-01', '2020-12-01'))
                    var = var.rename({'air': 'var'})

                var1 = xr.open_dataset(pwd + variables[v+1] + '.nc')
                var1 = var1.sel(time=slice(start + '-01-01', '2020-12-01'))
                print("open: " + variables[v+1] + '.nc')

            else:
                var = xr.open_dataset(pwd + variables[v] + '.nc')
                var = var.sel(time=slice(start + '-01-01', '2020-12-01'))
                print("open: " + variables[v] + '.nc')
                var1 = xr.open_dataset(pwd + variables[v + 1] + '.nc')
                var1 = var1.sel(time=slice(start + '-01-01', '2020-12-01'))
                print("open: " + variables[v + 1] + '.nc')

            start_in = int(start)
            aux = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
            n34 = Nino34CPC(aux, start=start_in)[2]
            del aux
            dmi = DMI(filter_bwa=bwa, per=0, start_per=str(start_in))[0]

            print('NO rescale')

            if cmap_revert[v]:
                cmap = cbar_r
                if v == 1:
                    cmap='BrBG'
            else:
                if v == 3:
                    cmap = 'jet'
                else:
                    cmap = cbar

            for season in range(season_start, season_end):

                main_month, main_month_name = len(seasons[:season]) + 1, seasons[season]
                if full_season:
                    main_month_name = 'JJASON'

                print(main_month_name)

                N34, N34_mmin, N34_mmax = SelectYears(df=n34, name_var='N34',
                                                      main_month=main_month, full_season=full_season)
                DMI2, DMI_mmin, DMI_mmax = SelectYears(df=dmi, name_var='DMI',
                                                       main_month=main_month, full_season=full_season)


                if (len(DMI2) != 0) & (len(N34) != 0):
                    DMI_pos, DMI_neg = ClassifierEvents(DMI2, full_season=False) # Full Season no necesario
                    N34_pos, N34_neg = ClassifierEvents(N34, full_season=False)
                    # Neutral events
                    # DMI_neutral = NeutralEvents(df=DMI2, mmin=DMI_mmin, mmax=DMI_mmax,
                    #                             var_original=var, start=start)
                    # N34_neutral = NeutralEvents(df=N34, mmin=N34_mmin, mmax=N34_mmax,
                    #                             var_original=var, start=start)
                    # both neutral, DMI and N34
                    All_neutral = NeutralEvents(df=DMI2, mmin=DMI_mmin, mmax=DMI_mmax, start=start,
                                                df2=N34, double=True, var_original=var)
                    All_neutral1 = NeutralEvents(df=DMI2, mmin=DMI_mmin, mmax=DMI_mmax, start=start,
                                                 df2=N34, double=True, var_original=var1)

                    # Simultaneous events
                    sim_events = np.intersect1d(N34.Años.values, DMI2.Años.values)

                    if (len(sim_events) != 0):
                        DMI_sim = DMI2.where(DMI2.Años.isin(sim_events)).dropna()
                        N34_sim = N34.where(N34.Años.isin(sim_events)).dropna()
                        DMI_sim_pos, DMI_sim_neg = ClassifierEvents(DMI_sim)
                        N34_sim_pos, N34_sim_neg = ClassifierEvents(N34_sim)

                        # Unique events
                        DMI_un = DMI2.where(-DMI2.Años.isin(sim_events)).dropna()
                        N34_un = N34.where(-N34.Años.isin(sim_events)).dropna()

                        DMI_un_pos, DMI_un_neg = ClassifierEvents(DMI_un)
                        N34_un_pos, N34_un_neg = ClassifierEvents(N34_un)

                        # ------------------------------------ SIMULTANEUS ---------------------------------------------#
                        DMI_sim = Composite(original_data=var, index_pos=DMI_sim_pos, index_neg=DMI_sim_neg,
                                            mmin=DMI_mmin, mmax=DMI_mmax)

                        DMI_sim1 = Composite(original_data=var1, index_pos=DMI_sim_pos, index_neg=DMI_sim_neg,
                                             mmin=DMI_mmin, mmax=DMI_mmax)

                        from matplotlib import colors
                        contour = False
                        Plots(data=DMI_sim, variable='var',
                              neutral=All_neutral,
                              DMI_pos=DMI_sim_pos, DMI_neg=DMI_sim_neg,
                              N34_pos=N34_sim_pos, N34_neg=N34_sim_neg,
                              mode='Simultaneus IODs-ENSOs',
                              title=title_var[v] + ' - ' + bwa_title + ' - ',
                              neutral_name='All_Neutral',
                              levels=scales[v], cmap=cmap,
                              name_fig=out_dir + 'SIM_' + bwa_fig + '_' + variables[v],
                              save=save, contour=contour, waf=wafs[v],
                              two_variables=True,  data2=DMI_sim1,
                              neutral2=All_neutral1, levels2=scales[v+1],
                              main_month_name=main_month_name,
                              SA=SA[v], step=step[v], contour0=contour0[v])

                        del DMI_sim
                        # ------------------------------------ UNIQUES ---------------------------------------------#
                        DMI_un = Composite(original_data=var, index_pos=DMI_un_pos, index_neg=DMI_un_neg,
                                           mmin=DMI_mmin, mmax=DMI_mmax)

                        N34_un = Composite(original_data=var, index_pos=N34_un_pos, index_neg=N34_un_neg,
                                           mmin=N34_mmin, mmax=N34_mmax)

                        DMI_un1 = Composite(original_data=var1, index_pos=DMI_un_pos, index_neg=DMI_un_neg,
                                            mmin=DMI_mmin, mmax=DMI_mmax)

                        N34_un1 = Composite(original_data=var1, index_pos=N34_un_pos, index_neg=N34_un_neg,
                                            mmin=N34_mmin, mmax=N34_mmax)

                        Plots(data=DMI_un, variable='var',
                              neutral=All_neutral,
                              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
                              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
                              mode='Isolated IODs',
                              title=title_var[v] + ' - ' + bwa_title + ' - ',
                              neutral_name='All_Neutral',
                              levels=scales[v], cmap=cmap,
                              name_fig=out_dir + 'UN_DMI_' + bwa_fig + '_' + variables[v],
                              save=save, contour=contour, waf=wafs[v],
                              two_variables=True,  data2=DMI_un1,
                              neutral2=All_neutral1, levels2=scales[v+1],
                              main_month_name=main_month_name,
                              SA=SA[v], step=step[v], contour0=contour0[v])

                        Plots(data=N34_un, variable='var',
                              neutral=All_neutral,
                              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
                              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
                              mode='Isolated ENSOs',
                              title=title_var[v] + ' - ' + bwa_title + ' - ',
                              neutral_name='All_Neutral',
                              levels=scales[v], cmap=cmap,
                              name_fig=out_dir + 'UN_ENSO_' + bwa_fig + '_' + variables[v],
                              save=save, contour=contour, waf=wafs[v],
                              two_variables=True,  data2=N34_un1,
                              neutral2=All_neutral1, levels2=scales[v+1],
                              main_month_name=main_month_name,text=text,
                              SA=SA[v],step=step[v], contour0=contour0[v])

                        del DMI_un, N34_un
                        # ------------------------------------ ALL ---------------------------------------------#
                        dmi_comp = Composite(original_data=var, index_pos=DMI_pos, index_neg=DMI_neg,
                                             mmin=DMI_mmin, mmax=DMI_mmax)
                        N34_comp = Composite(original_data=var, index_pos=N34_pos, index_neg=N34_neg,
                                             mmin=N34_mmin, mmax=N34_mmax)

                        dmi_comp1 = Composite(original_data=var1, index_pos=DMI_pos, index_neg=DMI_neg,
                                              mmin=DMI_mmin, mmax=DMI_mmax)
                        N34_comp1 = Composite(original_data=var1, index_pos=N34_pos, index_neg=N34_neg,
                                              mmin=N34_mmin, mmax=N34_mmax)

                    else:
                        print('NO simultanueos events, only uniques')

                        dmi_comp = Composite(original_data=var, index_pos=DMI_pos, index_neg=DMI_neg,
                                             mmin=DMI_mmin, mmax=DMI_mmax)
                        N34_comp = Composite(original_data=var, index_pos=N34_pos, index_neg=N34_neg,
                                             mmin=N34_mmin, mmax=N34_mmax)

                        dmi_comp1 = Composite(original_data=var1, index_pos=DMI_pos, index_neg=DMI_neg,
                                              mmin=DMI_mmin, mmax=DMI_mmax)
                        N34_comp1 = Composite(original_data=var1, index_pos=N34_pos, index_neg=N34_neg,
                                              mmin=N34_mmin, mmax=N34_mmax)

                    Plots(data=N34_comp, variable='var',
                          neutral=All_neutral,
                          DMI_pos=DMI_pos, DMI_neg=DMI_neg,
                          N34_pos=N34_pos, N34_neg=N34_neg,
                          mode='. ENSOs',
                          title=title_var[v] + ' - ' + bwa_title + ' - ',
                          neutral_name='All_Neutral',
                          levels=scales[v], cmap=cmap,
                          name_fig=out_dir + 'All_ENSO_' + bwa_fig + '_' + variables[v],
                          save=save, contour=contour, waf=wafs[v],
                          two_variables=True, data2=N34_comp1,
                          neutral2=All_neutral1, levels2=scales[v + 1],
                          main_month_name=main_month_name,text=text,
                          SA=SA[v],step=step[v], contour0=contour0[v])

                    Plots(data=dmi_comp, variable='var',
                          neutral=All_neutral,
                          DMI_pos=DMI_pos, DMI_neg=DMI_neg,
                          N34_pos=N34_pos, N34_neg=N34_neg,
                          mode='. IODs',
                          title=title_var[v] + ' - ' + bwa_title + ' - ',
                          neutral_name='All_Neutral',
                          levels=scales[v], cmap=cmap,
                          name_fig=out_dir + 'All_IOD_' + bwa_fig + '_' + variables[v],
                          save=save, contour=contour, waf=wafs[v],
                          two_variables=True, data2=dmi_comp1,
                          neutral2=All_neutral1, levels2=scales[v + 1],
                          main_month_name=main_month_name,
                          SA=SA[v],step=step[v], contour0=contour0[v])

                    del dmi_comp, N34_comp

                else:
                    print('skip' + main_month_name)

        else:

            var = xr.open_dataset(pwd + variables[v] + '.nc')
            print(variables[v] + '.nc')

            contour = contours[v]

            if rescale[v]:
                var = var.__mul__(rescale_value[v])

            if cmap_revert[v]:
                print('cmap revert')
                cmap = cbar_r
            else:
                print('cmap no revert')
                cmap = cbar


            start_in = int(start)
            aux = xr.open_dataset("/pikachu/datos4/Obs/sst/sst.mnmean_2020.nc")
            n34 = Nino34CPC(aux, start=start_in)[2]
            del aux
            dmi = DMI(filter_bwa=bwa, per=0, start_per=str(start_in))[0]

            for season in range(season_start, season_end):

                main_month, main_month_name = len(seasons[:season]) + 1, seasons[season]
                if full_season:
                    main_month_name = 'JJASON'

                print(main_month_name)

                N34, N34_mmin, N34_mmax = SelectYears(df=n34, name_var='N34',
                                                      main_month=main_month, full_season=full_season)
                DMI2, DMI_mmin, DMI_mmax = SelectYears(df=dmi, name_var='DMI',
                                                       main_month=main_month, full_season=full_season)

                if (len(DMI2) != 0) & (len(N34) != 0):
                    DMI_pos, DMI_neg = ClassifierEvents(DMI2, full_season=full_season)
                    N34_pos, N34_neg = ClassifierEvents(N34, full_season=full_season)

                    # Neutral events
                    # DMI_neutral = NeutralEvents(df=DMI2, mmin=DMI_mmin, mmax=DMI_mmax,
                    #                             var_original=var, start=start)
                    # N34_neutral = NeutralEvents(df=N34, mmin=N34_mmin, mmax=N34_mmax,
                    #                             var_original=var, start=start)
                    # both neutral, DMI and N34
                    All_neutral = NeutralEvents(df=DMI2, mmin=DMI_mmin, mmax=DMI_mmax, start=start,
                                                df2=N34, double=True, var_original=var)

                    # Simultaneous events
                    sim_events = np.intersect1d(N34.Años.values, DMI2.Años.values)

                    if (len(sim_events) != 0):
                        DMI_sim = DMI2.where(DMI2.Años.isin(sim_events)).dropna()
                        N34_sim = N34.where(N34.Años.isin(sim_events)).dropna()
                        DMI_sim_pos, DMI_sim_neg = ClassifierEvents(DMI_sim)
                        N34_sim_pos, N34_sim_neg = ClassifierEvents(N34_sim)

                        # Unique events
                        DMI_un = DMI2.where(-DMI2.Años.isin(sim_events)).dropna()
                        N34_un = N34.where(-N34.Años.isin(sim_events)).dropna()

                        DMI_un_pos, DMI_un_neg = ClassifierEvents(DMI_un)
                        N34_un_pos, N34_un_neg = ClassifierEvents(N34_un)

                        # ------------------------------------ SIMULTANEUS ---------------------------------------------#
                        DMI_sim = Composite(original_data=var, index_pos=DMI_sim_pos, index_neg=DMI_sim_neg,
                                            mmin=DMI_mmin, mmax=DMI_mmax)

                        from matplotlib import colors

                        Plots(data=DMI_sim, variable='var',
                              neutral=All_neutral,
                              DMI_pos=DMI_sim_pos, DMI_neg=DMI_sim_neg,
                              N34_pos=N34_sim_pos, N34_neg=N34_sim_neg,
                              mode='Simultaneus IODs-ENSOs',
                              title=title_var[v] + ' - ' + bwa_title + ' - ',
                              neutral_name='All_Neutral',
                              levels=scales[v], cmap=cmap,
                              name_fig=out_dir + 'SIM_' + bwa_fig + '_' + variables[v],
                              save=save, contour=contour, waf=wafs[v],
                              main_month_name=main_month_name)

                        del DMI_sim
                        # ------------------------------------ UNIQUES ---------------------------------------------#
                        DMI_un = Composite(original_data=var, index_pos=DMI_un_pos, index_neg=DMI_un_neg,
                                           mmin=DMI_mmin, mmax=DMI_mmax)

                        N34_un = Composite(original_data=var, index_pos=N34_un_pos, index_neg=N34_un_neg,
                                           mmin=N34_mmin, mmax=N34_mmax)

                        Plots(data=DMI_un, variable='var',
                              neutral=All_neutral,
                              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
                              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
                              mode='Isolated IODs',
                              title=title_var[v] + ' - ' + bwa_title + ' - ',
                              neutral_name='All_Neutral',
                              levels=scales[v], cmap=cmap,
                              name_fig=out_dir + 'UN_DMI_' + bwa_fig + '_' + variables[v],
                              save=save, contour=contour, waf=wafs[v],
                              main_month_name=main_month_name)

                        Plots(data=N34_un, variable='var',
                              neutral=All_neutral,
                              DMI_pos=DMI_un_pos, DMI_neg=DMI_un_neg,
                              N34_pos=N34_un_pos, N34_neg=N34_un_neg,
                              mode='Isolated ENSOs',
                              title=title_var[v] + ' - ' + bwa_title + ' - ',
                              neutral_name='All_Neutral',
                              levels=scales[v], cmap=cmap,
                              name_fig=out_dir + 'UN_ENSO_' + bwa_fig + '_' + variables[v],
                              save=save, contour=contour, waf=wafs[v],
                              main_month_name=main_month_name,text=text)

                        del DMI_un, N34_un
                        # ------------------------------------ ALL ---------------------------------------------#
                        dmi_comp = Composite(original_data=var, index_pos=DMI_pos, index_neg=DMI_neg,
                                             mmin=DMI_mmin, mmax=DMI_mmax)
                        N34_comp = Composite(original_data=var, index_pos=N34_pos, index_neg=N34_neg,
                                             mmin=N34_mmin, mmax=N34_mmax)

                    else:
                        print('NO simultanueos events, only uniques')

                        dmi_comp = Composite(original_data=var, index_pos=DMI_pos, index_neg=DMI_neg,
                                             mmin=DMI_mmin, mmax=DMI_mmax)
                        N34_comp = Composite(original_data=var, index_pos=N34_pos, index_neg=N34_neg,
                                             mmin=N34_mmin, mmax=N34_mmax)

                    Plots(data=N34_comp, variable='var',
                          neutral=All_neutral,
                          DMI_pos=DMI_pos, DMI_neg=DMI_neg,
                          N34_pos=N34_pos, N34_neg=N34_neg,
                          mode='. ENSOs',
                          title=title_var[v] + ' - ' + bwa_title + ' - ',
                          neutral_name='All_Neutral',
                          levels=scales[v], cmap=cmap,
                          name_fig=out_dir + 'All_ENSO_' + bwa_fig + '_' + variables[v],
                          save=save, contour=contour, waf=wafs[v],
                          main_month_name=main_month_name,text=text)

                    Plots(data=dmi_comp, variable='var',
                          neutral=All_neutral,
                          DMI_pos=DMI_pos, DMI_neg=DMI_neg,
                          N34_pos=N34_pos, N34_neg=N34_neg,
                          mode='. IODs',
                          title=title_var[v] + ' - ' + bwa_title + ' - ',
                          neutral_name='All_Neutral',
                          levels=scales[v], cmap=cmap,
                          name_fig=out_dir + 'All_IOD_' + bwa_fig + '_' + variables[v],
                          save=save, contour=contour, waf=wafs[v],
                          main_month_name=main_month_name)

                    del dmi_comp, N34_comp

                else:
                    print('skip' + main_month_name)

