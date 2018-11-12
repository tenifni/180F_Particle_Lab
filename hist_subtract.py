import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

plt.ion()

data = np.loadtxt('27FebTargetAl_3.txt')
bkgd = np.loadtxt('8MarNoTarg.txt')

# transpose, for indexing
data = data[:1000].T # only use the first 1000; it gets weird after this
bkgd = bkgd.T

# Convert to ns
data[:-1] *= 4
bkgd[:-1] *= 4

# column i + 1
i=2

# make cuts on tdc values
hdata = data[i][data[i]>150][:1000]
bdata = bkgd[i][bkgd[i]>150]

# bin data
xmax = 2e4
nbins = 80
dhist, dbins = np.histogram(hdata,bins=nbins,range=(0,xmax))
bhist, bbins = np.histogram(bdata,bins=nbins,range=(0,xmax))
binsz = np.median(np.diff(dbins)) # this is the binsize

# calculate run lengths
dlen = (data[-1][-1] - data[-1][0])
blen = (bkgd[-1][-1] - bkgd[-1][0])

# subtract background
dif = dhist/dlen - bhist/blen

# fit
def expo(x,p):
	return p[0] * np.exp(-x/p[1]) + p[2]
p0 = [1e-2,2.2e3,1e-4] # initial guess
# inverse variance
ivar = 1./dif
ivar[dif <= 0] = 0
func = lambda p: np.sum((dif-expo(bbins[:-1],p))**2*ivar)
res = minimize(func,x0=p0,method='Nelder-Mead')

# plot
plt.step(bbins[:-1],dif,where='post',color='k',label='bkgd')
plt.step(bbins[:-1],expo(bbins[:-1],res.x),'r-',label='fit')
plt.yscale('log')
plt.ylim(1e-4,1e-1)
plt.xlim(0,xmax)
plt.legend()
