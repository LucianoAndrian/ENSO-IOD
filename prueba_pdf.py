import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.stats
size = len(case)
hist, bins = np.histogram(case)
y = case
#h = plt.bar(bins[:-1], hist.astype(np.float32) / hist.sum(), width=(bins[1] - bins[0]), color='grey')
#h = plt.hist(y, density=True)

dist_names = ['gamma', 'beta', 'norm']

for dist_name in dist_names:
    dist = getattr(scipy.stats, dist_name)
    params = dist.fit(y)
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)

    if arg:
        pdf_fitted = dist.pdf(x, *arg, loc=loc, scale=scale)*size
    else:
        pdf_fitted = dist.pdf(x, loc=loc, scale=scale)*size

    pdf = pd.Series(pdf_fitted, x)
    plt.plot(pdf, label=dist_name)
    #plt.xlim(0,47)
plt.legend(loc='upper right')
plt.show()