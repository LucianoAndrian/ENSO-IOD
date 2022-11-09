from ENSO_IOD_Funciones import DMI
import matplotlib.pyplot as plt
import numpy as np

i = 1920
end = 2020
fs = False
dmi, aux, dmi_aux = DMI(filter_bwa=False, start_per=i, end_per=end,
                        sst_anom_sd=False, opposite_signs_criteria=False)
dmi = dmi.iloc[np.where(dmi.Años>1949)]

annual_dmi =  np.zeros(shape=(70,12))
i_count=0
for i in np.arange(1950,2020):
    annual_dmi[i_count,:] = dmi_aux.sel(time=dmi_aux.time.dt.year.isin(i))
    i_count += 1

sd = np.std(annual_dmi, axis=0)
mean = np.mean(annual_dmi, axis=0)

sd_eventos = []
for i in np.arange(1,13):
    sd_eventos.append(np.std(dmi.iloc[np.where(dmi.Mes==i)[0]]['DMI']))


cantidad = []
for i in np.arange(1,13):
    cantidad.append(len(dmi.iloc[np.where(dmi.Mes==i)]))

plt.style.use('dark_background')
fig = plt.figure(1, figsize=(5, 3), dpi=200)
ax = fig.add_subplot(111)
ax2 = ax.twinx()
for i in np.arange(1950,2020):
    ax.plot(dmi_aux.sel(time=dmi_aux.time.dt.year.isin(i)), zorder=0, alpha=0.6)

ax.plot(mean, color='k', zorder=20, linewidth=2)
ax.fill_between(np.arange(0,12),sd/2, -sd/2, alpha=0.9, zorder=15, color='red')
# ax.fill_between(np.arange(0,12),
#                  [a for a in sd_eventos],
#                  [-a for a in sd_eventos],
#                  color='green', alpha=0.9, zorder=10)
ax.set_ylabel('DMI')
ax.set_xlabel('Meses')
ax.set_xticks(np.arange(0,12), np.arange(1,13))
ax.set_ylim(-1.5, 1.5)

ax2.bar(np.arange(0,12),cantidad, alpha=0.35, color='white')
ax2.set_ylim(0,40)

ax.grid()
plt.tight_layout()
#plt.savefig()



plt.show()






#plt.plot(mean, color='k')
plt.fill_between(np.arange(0,12),
                 [a for a in sd_eventos],
                 [-a for a in sd_eventos],
                 color='red', alpha=0.5)
plt.fill_between(np.arange(0,12),sd, -sd, alpha=0.5, zorder=1)
#plt.bar(np.arange(0,12),cantidad)
plt.show()










