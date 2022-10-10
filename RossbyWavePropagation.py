import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.feature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cartopy.crs as ccrs
import numpy as np
from matplotlib import colors
import xarray as xr

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
            print('1 year')
            comp_field = comp_field.drop_dims(['time'])

        return comp_field
    else:
        print(' len index = 0')

def Plot(comp, levels = np.linspace(-1,1,11), cmap='RdBu',
         SA=False, dpi=100, save=True, step=1,
         name_fig='fig', title='title', contour0=True,):

    from numpy import ma
    import matplotlib.pyplot as plt

    comp_var = comp['var']

    if SA:
        fig = plt.figure(figsize=(5, 6), dpi=dpi)
    else:
        fig = plt.figure(figsize=(7, 3.5), dpi=dpi)

    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    crs_latlon = ccrs.PlateCarree()
    if SA:
        ax.set_extent([270,330, -60,20], crs_latlon)
    else:
        ax.set_extent([0, 359, -80, 40], crs=crs_latlon)


    im = ax.contourf(comp.lon[::step], comp.lat[::step], comp_var[::step,::step],
                     levels=levels,transform=crs_latlon, cmap=cmap, extend='both')
    if contour0:
        ax.contour(comp.lon, comp.lat, comp_var, levels=0,
                   transform=crs_latlon, colors='green', linewidths=1)

    cb = plt.colorbar(im, fraction=0.042, pad=0.035,shrink=0.8)
    cb.ax.tick_params(labelsize=8)
    #ax.add_feature(cartopy.feature.LAND, facecolor='#d9d9d9')
    ax.add_feature(cartopy.feature.LAND, facecolor='k')
    ax.add_feature(cartopy.feature.COASTLINE)
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

    plt.title(title, fontsize=10)
    # if text:
    #     plt.figtext(0.5, 0.01, 'Number of events: ' + str(number_events), ha="center", fontsize=10,
    #             bbox={'facecolor': 'w', 'alpha': 0.5, 'pad': 5})
    plt.tight_layout()

    if save:
        plt.savefig(out_dir + name_fig + '.jpg')
        plt.close()
    else:
        plt.show()
########################################################################################################################

data_dir2 = '/datos/luciano.andrian/ncfiles/ERA5/'
nc_date_dir = '/datos/luciano.andrian/ncfiles/nc_composites_dates/'
out_dir = '/home/luciano.andrian/doc/salidas/ENSO_IOD/composite/DMI_ind/RW/'


seasons = ("Full_Season",'SON', 'ASO', 'JAS', 'JJA')
min_max_months = [[7,11],[9,11], [8,10], [7,9], [6,8]]

cases = ['DMI_sim_pos', 'DMI_sim_neg', 'DMI_neg', 'DMI_pos', 'DMI_un_pos', 'DMI_un_neg',
         'N34_pos', 'N34_neg', 'N34_un_pos', 'N34_un_neg']

variables = ['RWS', 'u', 'etay', 'div_from_W']

scales = [np.linspace(-2e-10, 2e-10, 13),
          np.linspace(-12, 12, 13),
          np.linspace(-3e-11, 3e-11, 13),
          np.linspace(-0.5e-5, 0.5e-5, 13)]

change_names = ['__xarray_dataarray_variable__', None,
                'meridional_gradient_of_absolute_vorticity',
                'divergence']

v_name=['RWS', 'U-200hPa', 'ETAy', 'Divergence-200hPa']

########################################################################################################################
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
########################################################################################################################


c_var=0
save=True
step=4
i = '1950'
for v in variables:
    print(v + '.nc')
    data = xr.open_dataset(data_dir2 + v + '.nc')
    if v != 'u':
        data = data.rename({change_names[c_var]: 'var'})
    count=0
    for s in seasons:
        print(s)

        aux = xr.open_dataset(nc_date_dir + 'Composite_' + i + '_2020_' + s + '.nc')
        neutro = aux.Neutral

        dmi_un_neg = aux.DMI_un_neg
        dmi_un_pos = aux.DMI_un_pos
        dmi_n34_sim_neg = aux.DMI_sim_neg
        dmi_n34_sim_pos = aux.DMI_sim_pos
        aux.close()

        mmonth = min_max_months[count]
        mmin = mmonth[0]
        mmax = mmonth[-1]

        print(mmin)
        print(mmax)

        neutro_comp = CompositeSimple(original_data=data, index=neutro,
                                      mmin=mmin, mmax=mmax)

        for y in dmi_un_neg.values:
            data_y = CompositeSimple(original_data=data, index=[y],
                                     mmin=mmin, mmax=mmax)
            comp = data_y - neutro_comp
            if len(comp) != 0:
                Plot(comp=comp, levels=scales[c_var], cmap=cbar,
                     SA=False, dpi=200, save=save, step=step,
                     name_fig=v + '_DMI_un_neg_' + str(int(y)) + '_' + s + '_' + i + '_2020',
                     title=str(int(y)) + ' - Negative Isolated DMI - ' + v_name[c_var] + '\n' + s + ' ' + i + '-2020',
                     contour0=False)

        for y in dmi_un_pos.values:
            data_y = CompositeSimple(original_data=data, index=[y],
                                     mmin=mmin, mmax=mmax)
            comp = data_y - neutro_comp
            if len(comp) != 0:
                Plot(comp=comp, levels=scales[c_var], cmap=cbar,
                     SA=False, dpi=200, save=save, step=step,
                     name_fig=v + '_DMI_un_pos_' + str(int(y)) + '_' + s + '_' + i + '_2020',
                     title=str(int(y)) + ' - Positive Isolated DMI - ' + v_name[c_var] + '\n' + s + ' ' + i + '-2020',
                     contour0=False)


        for y in dmi_n34_sim_pos.values:
            data_y = CompositeSimple(original_data=data, index=[y],
                                     mmin=mmin, mmax=mmax)
            comp = data_y - neutro_comp
            if len(comp) != 0:
                Plot(comp=comp, levels=scales[c_var], cmap=cbar,
                     SA=False, dpi=200, save=save, step=step,
                     name_fig=v + '_DMI_N34_sim_pos_' + str(int(y)) + '_' + s + '_' + i + '_2020',
                     title=str(int(y)) + ' - Positive DMI and N34 - ' + v_name[c_var] + '\n' + s + ' ' + i + '-2020',
                     contour0=False)


        for y in dmi_n34_sim_neg.values:
            data_y = CompositeSimple(original_data=data, index=[y],
                                     mmin=mmin, mmax=mmax)
            comp = data_y - neutro_comp
            if len(comp) != 0:
                Plot(comp=comp, levels=scales[c_var], cmap=cbar,
                     SA=False, dpi=200, save=save, step=step,
                     name_fig=v + '_DMI_N34_sim_neg_' + str(int(y)) + '_' + s + '_' + i + '_2020',
                     title=str(int(y)) + ' - Negative DMI and N34 - ' + v_name[c_var] + '\n' + s + ' ' + i + '-2020',
                     contour0=False)

        count += 1
    c_var += 1


