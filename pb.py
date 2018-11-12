import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit,minimize

binsize = 20      # Size of histogram bins, in TDC time (ns/4)
cut = 0           # Don't fit times below this value
channel = 6       # TDC 4,5,6, or 7
cl = 3            # Number of neutron counter coincidences we require
spike = [200,300] # Range where spike occurs

bins = np.arange(0,4000,binsize)
taumu = 2197/4. # muon decay lifetime in tdc units

# get background histogram
bg_tdcs = np.loadtxt('data/12MarNoTarget.txt').T
bg = bg_tdcs[channel]
timestamps = bg_tdcs[-1]
mask = (bg > 0) * (np.sum((bg_tdcs[4:8] > 0),axis=0) >= cl)
bg = bg[mask]
timestamps = timestamps[mask]
bg_runlen = float(timestamps[-1]-timestamps[0])
# bin
bin_vals,bin_edges = np.histogram(bg,bins=bins)
bg_y = bin_vals.astype(float)
bg_x = bin_edges[:-1]
# guess error assuming Poisson statistics
bg_ysig = np.sqrt(bg_y)
bg_ysig[bg_y == 0] = np.inf

# get data histogram - similar to above
files = ['data/11MarTargetPb.txt',
		 'data/14MarTargetPb.txt']
pbs = []
pb_runlen = 0
for i in range(len(files)):
	pb_tdcs = np.loadtxt(files[i]).T
	pb = pb_tdcs[channel]
	timestamps = pb_tdcs[-1]
	mask = (pb > 0) * (np.sum((pb_tdcs[4:8] > 0),axis=0) >= cl)
	mask *= (pb_tdcs[min([channel,5])] > pb_tdcs[4])
	pb = pb[mask]
	timestamps = timestamps[mask]
	pb_runlen += float(timestamps[-1]-timestamps[0])
	pbs = np.append(pbs,pb)
# bin
bin_vals,bin_edges = np.histogram(pbs,bins=bins)
pb_y = bin_vals.astype(float)
pb_x = bin_edges[:-1]
# guess error assuming Poisson statistics
pb_ysig = np.sqrt(pb_y)
pb_ysig[pb_y == 0] = np.inf

# scale bg to pb using their run lengths
bg_y *= pb_runlen / bg_runlen
bg_ysig *= pb_runlen / bg_runlen

# subtract
y = pb_y - bg_y
ysig = np.sqrt(pb_ysig**2+bg_ysig**2)
x = bg_x+binsize/2.

# fit
expo = lambda x, a, b, tau, offset: a*np.exp(-x/tau)+b*np.exp(-x/taumu)+offset
p0 = [np.max(y), np.max(y)/4., 20, 0]
mask = (x < spike[0]) * (x > cut) + (x > spike[1]) * (x < 2000)
xfit,yfit,efit = x[mask], y[mask], ysig[mask]
ivar = 1./efit**2

minfunc = lambda p: np.sum((yfit - expo(xfit,*p))**2*ivar)
res = minimize(minfunc,x0=p0,method='Nelder-Mead')
p1 = res.x

pfit,pcov = curve_fit(expo, xfit, yfit, p1, efit)
perr = np.sqrt(np.diag(pcov))
print p0,p1,pfit

# plot
plt.subplot(211)
plt.errorbar(pb_x+binsize/2., pb_y, pb_ysig, label='Pb[{0}]'.format(channel),
	xerr=binsize/2,ls='none',capsize=0,markersize=0)
plt.errorbar(bg_x+binsize/2., bg_y, bg_ysig, label='BG[{0}]'.format(channel),
	xerr=binsize/2,ls='none',capsize=0,markersize=0)
plt.axvspan(spike[0],spike[1], alpha=0.1, color='red')
plt.ylim(0.5, 2e2)
plt.xlim(0,4000)
plt.legend()

plt.subplot(212)
plt.errorbar(bg_x+binsize/2., y, ysig, xerr=binsize/2., label='Pb-BG',
	ls='none',capsize=0,markersize=0,color='k')
plotx = np.linspace(min(xfit),max(xfit),1000)
plt.plot(plotx,expo(plotx,*pfit),ls='-',label='Total Fit')
plt.plot(plotx,pfit[0]*np.exp(-plotx/pfit[2]),ls='-',label='mu-')
plt.plot(plotx,pfit[1]*np.exp(-plotx/taumu),ls='-',label='mu+')
plt.axhline(pfit[3],label='tail',color='m')
plt.axvspan(spike[0],spike[1], alpha=0.2, color='red')
plt.text(1e2,1e2,'tau = {0:.1f} +- {1:.1f}'.format(4*pfit[2],4*perr[2]))
plt.ylim(-20, 250)
plt.xlim(0,500)
plt.xlabel('TDC time (ns/4)')
plt.legend()
plt.tight_layout()






