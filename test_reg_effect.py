import numpy as np
import matplotlib.pyplot as plt

from scipy.stats import pearsonr
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
# ---------------------------------------------------------------------------- #
def crear_serie_1(n_steps, b0, f0):

    t = np.linspace(0, 2 * np.pi, n_steps)
    serie = b0 * np.cos(f0 * t)
    return serie

def crear_serie_2(n_steps, b1, f1):

    t = np.linspace(0, 2 * np.pi, n_steps)
    serie = b1 * np.sin(f1 * t)
    return serie

def regre(series, intercept, coef=0):
    df = pd.DataFrame(series)
    if intercept:
        X = np.column_stack((np.ones_like(df[df.columns[1]]),
                             df[df.columns[1:]]))
    else:
        X = df[df.columns[1:]].values
    y = df[df.columns[0]]

    coefs = np.linalg.lstsq(X, y, rcond=None)[0]

    coefs_results = {}
    for ec, e in enumerate(series.keys()):
        if intercept and ec == 0:
            e = 'constant'
        if e != df.columns[0]:
            if intercept:
                coefs_results[e] = coefs[ec]
            else:
                coefs_results[e] = coefs[ec-1]

    if isinstance(coef, str):
        return coefs_results[coef]
    else:
        return coefs_results

def create_series(rs):

    np.random.seed(rs)

    n_steps = 100
    serie1 = crear_serie_1(n_steps, 1, 4) + 1 * np.random.normal(0, 1, n_steps)
    n34 = serie1 / serie1.std()

    serie2 = crear_serie_1(n_steps, 1, 4) + 1 * np.random.chisquare(1, n_steps)
    dmi = serie2 / serie2.std()

    serie = dmi + 3 * np.random.chisquare(1, n_steps)
    strato = serie / serie.std()

    return n34, dmi, strato
# ---------------------------------------------------------------------------- #
random_seed = []
random_seed_f = []
for rs_count, rs in enumerate(range(0,1000)):
    n34, dmi, strato = create_series(rs)
    seriesf = {'strato': strato, 'n34': n34, 'dmi': dmi}
    seriesn34 = {'strato': strato, 'n34': n34}
    seriesdmi = {'strato': strato, 'dmi': dmi}

    b_total = abs(regre(seriesn34, True, 'n34'))
    b_directo = abs(regre(seriesf, True, 'n34'))

    if (b_total*2 < b_directo):
        random_seed.append(rs)
        if (b_total*4 < b_directo):
            random_seed_f.append(rs)

print(f"se cumple un {round(100*len(random_seed)/rs_count,1)}% de las veces")


rs = random_seed_f[2]

n34, dmi, strato = create_series(rs)

seriesf = {'strato': strato, 'n34': n34, 'dmi': dmi}
seriesn34 = {'strato': strato, 'n34': n34}
seriesdmi = {'strato': strato, 'dmi': dmi}

plt.plot(strato, label='strato', color='k')
plt.plot(n34, label='n34', color='red')
plt.plot(dmi, label='dmi', color='green')
plt.legend()
plt.show()
# ---------------------------------------------------------------------------- #
print(f"Corr n34-dmi r={round(pearsonr(n34, dmi)[0],3)} "
      f"p-value={round(pearsonr(n34, dmi)[1],3)}")

print(f"Corr n34-strato r={round(pearsonr(n34, strato)[0],3)} "
      f"p-value={round(pearsonr(n34, strato)[1],3)}")

print(f"Corr strato-dmi r={round(pearsonr(strato, dmi)[0],3)} "
      f"p-value={round(pearsonr(strato, dmi)[1],3)}")
print('')
print(f"Efecto Total N34: {round(regre(seriesn34, True, 'n34'), 3)}")
print(f"Efecto Directo N34, ajustando DMI: "
      f"{round(regre(seriesf, True, 'n34'), 3)}")
print('')
print(f"Efecto Total DMI: {round(regre(seriesdmi, True, 'dmi'), 3)}")
print(f"Efecto Directo DMI, ajustando N34: "
      f"{round(regre(seriesf, True, 'dmi'), 3)}")
# ---------------------------------------------------------------------------- #

print('Variance Inflation Factor (VIF)')
df = pd.DataFrame({
    #'strato': strato,
    'n34': n34,
    'dmi': dmi
})
X = sm.add_constant(df)
vif = pd.DataFrame()
vif['Variable'] = X.columns
vif['VIF'] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
print(vif)
if (vif['VIF'].values>5).mean() == 0:
    print('No existe multicolinealidad')
# ---------------------------------------------------------------------------- #