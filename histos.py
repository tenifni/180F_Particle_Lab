import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Generate 2000 random numbers, drawn from a normal distribution with
# mean 0 and standard deviation 1
N  = int(2e3)
g1 = np.random.normal(0, 1, N)

# Save g1 in a txt file
np.savetxt('histos.txt', g1, delimiter=' ')

# Load data
data = np.loadtxt('histos.txt')

# Make 1D histogram
hist, bin_edges = np.histogram(data, bins=50, range=(-5,5))
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.

# Estimate errors; ignore bins with 0 counts
errors = np.sqrt(hist)
errors[hist == 0] = np.inf

# Fit histogram and print parameters
def gauss(x,  amplitude, mean, sigma):
    return amplitude*np.exp(-0.5*((x-mean)/sigma)**2)
p0 = [max(hist), np.mean(data), np.std(data)]
optimal_pars, par_covariance = curve_fit(gauss, bin_centers, hist, p0, errors)
par_error = np.sqrt(np.diag(par_covariance))

# Print best-fit parameters
print 'Best-fit parameters:'
par_names = ['  A', ' mu', 'sig']
for i in range(3):
	print par_names[i],
	print '{0:.3f}'.format(optimal_pars[i]).rjust(8), '+-', end='')
	print '{0:.3f}'.format(par_error[i]).ljust(8)

# Plot + save output
plt.hist(g1,bin_centers,histtype='step',color='k')
plt.plot(bin_centers,gauss(bin_centers,*optimal_pars),'r-')
plt.xlim(-5,5)
plt.xlabel('X-values')
plt.ylabel('N')
plt.savefig('histos.pdf')
plt.show()

