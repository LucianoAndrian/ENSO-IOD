# P0
#----------------------------------------------------------------------------------------------------------------------#
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#----------------------------------------------------------------------------------------------------------------------#
data_dir = '/home/luciano.andrian/doc/climaext/p0/data/'
out_dir = '/home/luciano.andrian/doc/climaext/p0/'
save = False
#----------------------------------------------------------------------------------------------------------------------#
# 2
temp = pd.read_csv(data_dir + 'edc3deuttemp2007.txt', sep='\s+', skiprows=91) # \s+ separador de más de 1 espacio
temp_age = temp.loc[:,'Age']
temp_data = temp.loc[:,'Temperature']
co2 =  pd.read_csv(data_dir + 'edc3-composite-co2-2008-noaa.txt', sep='\s+', skiprows=275)
co2_age = co2.loc[:,'gas_ageBP']
co2_data = co2.loc[:,'CO2']

# Plot
fig = plt.figure(1, figsize=(8, 4), dpi=300)
ax = fig.add_subplot(111)
ax2 = ax.twinx()
co2_tl = ax.plot(co2_age, co2_data, label='CO2 [ppm]', color='dodgerblue')
temp_tl = ax2.plot(temp_age, temp_data, label='Temp. [ºC]', color='red')

ax.set_xlim(0,800000)
xticks = np.arange(0,9e5,1e5)
xlab=[]
for c,x in enumerate(xticks):
    xstr=format(x, '.0e').split('e')
    xlab.append(xstr[0] + '$x10^{5}$')
ax.set_xticks(xticks, xlab)
ax.set_ylim(0, 350)
ax2.set_ylim(-12,8)

lns = co2_tl + temp_tl
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)

if save:
    plt.savefig(out_dir + 'p0_c.jpg')
    print('Save')
    plt.close('all')
else:
    plt.show()


# time step check
def CheckTimeStep(serie):
    diff = []
    for i in range(0, len(serie)):
        t0 = np.abs(serie[i])
        try:
            t1 = np.abs(serie[i + 1])
        except:
            return diff
        if t0 > t1:
            t_dif = t0 - t1
        else:
            t_dif = t1 - t0

        diff.append(t_dif)

plt.plot(CheckTimeStep(co2_age));plt.show()
plt.plot(CheckTimeStep(temp_age));plt.show()