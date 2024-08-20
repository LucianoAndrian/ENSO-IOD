"""
Que tan bien el modelo simula el patrón del ENSO en comparación con su patrón
historico?
"""
################################################################################
save = True
################################################################################
out_dir = '/home/luciano.andrian/doc/RegreVerif/'
data_dir = '/pikachu/datos/luciano.andrian/cases_fields/'
index_dir = '/pikachu/datos/luciano.andrian/DMI_N34_Leads_r/'
################################################################################
# import
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from ENSO_IOD_Funciones import LinearReg

if save:
    dpi=300
else:
    dpi=100

################################################################################
def SpatialCorr(campo1, campo2):
    stacked_campo1 = campo1.stack(dim_0=('lon', 'lat'))
    stacked_campo2 = campo2.stack(dim_0=('lon', 'lat'))

    # mascara para identificar los NaN
    mask = np.isfinite(stacked_campo1) & np.isfinite(stacked_campo2)
    campo1_valid = stacked_campo1[mask]
    campo2_valid = stacked_campo2[mask]

    correlacion = np.corrcoef(campo1_valid, campo2_valid)[0, 1]

    return correlacion

def Compute_em(data_em, nombre_var, n34_em, name_fig, save=save, dpi=dpi):

    data_em = data_em.rename({nombre_var:'var'})

    rf = []
    fig = plt.figure(figsize=(7, 6), dpi=dpi)
    ax = fig.add_subplot(111)
    for l, c in zip([0, 1, 2, 3],
                    ['#D81A06', '#2C95FF', 'forestgreen', '#CB5FB1']):

        data_em_L = data_em.sel(time=data_em.time.dt.month.isin(10 - l))
        n34_em_L = n34_em.sel(time=n34_em.time.dt.month.isin(10 - l))

        rl = []
        for t in data_em_L.time.values:
            aux_data_em_L = data_em_L.drop_sel(time=[t])
            aux_n34_em_L = n34_em_L.drop_sel(time=[t]).dropna(dim='time')

            aux_data_em_L = aux_data_em_L / aux_data_em_L.std('time')
            aux_n34_em_L = aux_n34_em_L / aux_n34_em_L.std('time')

            data_to_verif = data_em_L.sel(time=t) / aux_data_em_L.std('time')
            n34_to_verif = n34_em_L.sel(time=t) / aux_n34_em_L.std('time')

            # aux_data_em_L['time'] = aux_n34_em_L.sst.values
            aux_data_em_L = aux_data_em_L.sel(time=aux_n34_em_L.time.values)
            aux_data_em_L['time'] = aux_n34_em_L.sst.values

            aux_regre = LinearReg(aux_data_em_L, 'time')

            aux_regre_pat = (aux_regre.var_polyfit_coefficients[0] +
                             aux_regre.var_polyfit_coefficients[1])

            rl.append(SpatialCorr(
                aux_regre_pat * np.abs(n34_to_verif.sst), data_to_verif['var']))

        rf.append(rl)

        ax.scatter(n34_em_L.sst.values, rl, s=10,
                   color=c, marker='o', alpha=0.8, label=f"Lead {l}")

    ax.set_ylim(-1, 1)
    ax.set_xlim(-3, 3)
    ax.grid()
    ax.hlines(y=0, xmin=-5, xmax=5, color='k', linestyle='--')
    ax.vlines(x=0, ymin=-1, ymax=1, color='k', linestyle='--')
    plt.legend()
    plt.xlabel('ONI')
    plt.ylabel(f" {nombre_var} pattern correlation")
    plt.tight_layout()
    if save:
        plt.savefig(f"{out_dir}{name_fig}.jpg")
        plt.close()
    else:
        plt.show()

def Compute_r(data, nombre_var, n34, name_fig, save=save, dpi=dpi):
    data = data.rename({nombre_var: 'var'})

    fig = plt.figure(figsize=(7, 6), dpi=dpi)
    ax = fig.add_subplot(111)
    for r in range(1, 25):
        print(r)
        data_r = data.sel(r=r)
        n34_r = n34.sel(r=r)

        for l, c in zip([0, 1, 2, 3],
                        ['#D81A06', '#2C95FF',
                         'forestgreen', '#CB5FB1']):

            data_r_L = data_r.sel(time=data_r.time.dt.month.isin(10 - l))
            n34_r_L = n34_r.sel(time=n34_r.time.dt.month.isin(10 - l))

            rl = []
            for t in data_r_L.time.values:
                aux_data_r_L = data_r_L.drop_sel(time=[t])
                aux_n34_r_L = n34_r_L.drop_sel(time=[t]).dropna(dim='time')

                aux_data_r_L = aux_data_r_L / aux_data_r_L.std('time')
                aux_n34_r_L = aux_n34_r_L / aux_n34_r_L.std('time')

                data_to_verif = data_r_L.sel(time=t) / aux_data_r_L.std('time')
                n34_to_verif = n34_r_L.sel(time=t) / aux_n34_r_L.std('time')

                aux_data_r_L = aux_data_r_L.sel(time=aux_n34_r_L.time.values)
                aux_data_r_L['time'] = aux_n34_r_L.sst.values
                # try:

                aux_regre = LinearReg(aux_data_r_L, 'time')

                aux_regre_pat = (aux_regre.var_polyfit_coefficients[0] +
                                 aux_regre.var_polyfit_coefficients[1])

                rl.append(SpatialCorr(
                    aux_regre_pat * np.abs(n34_to_verif.sst),
                    data_to_verif['var']))

            ax.scatter(n34_r_L.sst.values, rl, s=10,
                       color=c, marker='o', alpha=0.8)

    ax.set_ylim(-1, 1)
    ax.set_xlim(-3, 3)
    ax.grid()
    ax.hlines(y=0, xmin=-5, xmax=5, color='k', linestyle='--')
    ax.vlines(x=0, ymin=-1, ymax=1, color='k', linestyle='--')
    plt.xlabel('ONI')
    plt.ylabel(f"{nombre_var} pattern correlation")
    plt.tight_layout()
    if save:
        plt.savefig(f"{out_dir}{name_fig}.jpg")
        plt.close()
    else:
        plt.show()
################################################################################

n34 = xr.open_dataset(index_dir + 'N34_SON_Leads_r_CFSv2.nc')

# PREC ----------------------------------------------------------------------- #
data = xr.open_dataset(data_dir + 'prec_son.nc')
data_em = data.mean('r')
n34_em = n34.mean('r')

Compute_em(data_em, 'prec', n34_em, name_fig='prec_sacatter_em')
Compute_r(data, 'prec', n34, name_fig='prec_sacatter_r')

# Tref ----------------------------------------------------------------------- #
data = xr.open_dataset(data_dir + 'tref_son.nc').sel(lon=slice(275,330),
                                                     lat=slice(-60,15))

data_em = data.mean('r')
n34_em = n34.mean('r')

Compute_em(data_em, 'tref', n34_em, name_fig='tref_sacatter_em')
Compute_r(data, 'tref', n34, name_fig='tref_sacatter_r')

# HGT ------------------------------------------------------------------------ #
data = (xr.open_dataset(data_dir + 'hgt_son.nc').__mul__(9.80665)
        .sel(lon=slice(200,340), lat=slice(-80,-20)))

data_em = data.mean('r')
n34_em = n34.mean('r')

Compute_em(data_em, 'hgt', n34_em, name_fig='hgt_sacatter_em')
Compute_r(data, 'hgt', n34, name_fig='hgt_sacatter_r')

################################################################################