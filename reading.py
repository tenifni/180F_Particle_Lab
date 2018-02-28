import numpy as np
import matplotlib.pyplot as plt

'''
mu+ example
'''
muplus_data = np.loadtxt('data/muplus.txt')
muplus_data = muplus_data.T # transpose, for indexing
# Description of data indices:
# 0              1       2      3         4
# CAMAC_SLOT_ID  NUMBER  NEVENT TDC_VALUE UNIX_TIME

# Convert to ns
muplus_data[3] *= 40
# Plot a histogram of TDC_VALUE
fig=plt.figure(figsize=(8,3))
plt.subplot(121)
plt.hist(muplus_data[3],bins=50,range=(0,1e4),histtype='step',color='k')
plt.xlabel('Time (ns)')
plt.ylabel('N')
# Plot TDC_VALUE vs NUMBER
plt.subplot(122)
plt.plot(muplus_data[1],muplus_data[3],'k.',alpha=0.3)
plt.xlim(0,max(muplus_data[1]))
plt.ylim(0,1e4)
plt.ylabel('Time (ns)')
plt.tight_layout()
plt.savefig('muplus.pdf')
plt.close()

'''
mu- example
'''
muminus_data = np.loadtxt('data/muminus.txt')
muminus_data = muminus_data.T # transpose, for indexing
# Description of data indices:
# 0    1    2    3    4    5    6    7    8
# TDC1 TDC2 TDC3 TDC4 TDC5 TDC6 TDC7 TDC8 UNIX_TIME

# Convert to ns
muplus_data[:-1] *= 4
# Plot histograms of all TDCs
fig=plt.figure(figsize=(8,6))
for i in range(8):
	plt.subplot(331+i)
	plt.hist(muminus_data[i],bins=50,range=(0,4e3),histtype='step',color='k')
	plt.ylim(0, max( [3, plt.ylim()[1]] ))
	plt.title('TDC {0}'.format(i+1))
plt.tight_layout()
plt.savefig('muminus.pdf')
plt.close()

'''
MCS example
'''
from os.path import isfile
if isfile('mcs.npy'):             # Loading from the txt file takes a while,
	mcs_data = np.load('mcs.npy') # so we'll use a npy file in later calls
else:
	print 'Loading MCS data! Zzzz...'
	mcs_data = np.loadtxt('data/mcs.txt').astype(int)
	mcs_data = mcs_data.T # transpose, for indexing
	np.save('mcs.npy',mcs_data)
# Description of data indices:
# 0  1  2
# x  y  data

# Convert data to array of images
x,y,z = mcs_data
data_len = np.count_nonzero((x+y)==0)
imgs = np.zeros((data_len,640,480))
img_num = 0
for i in xrange(1,len(x)):
	if x[i] == 0 and y[i] == 0:
		img_num += 1
		continue
	imgs[img_num][x[i]][y[i]] = z[i]
# Plot first picture, with cool color map
# cmap ref: https://matplotlib.org/examples/color/colormaps_reference.html
fig=plt.figure(figsize=(6,6))
plt.imshow(imgs[0],cmap='copper')
plt.colorbar()
plt.xlim(0,480)
plt.ylim(640,0)
plt.title('#nofilter')
plt.savefig('mcs.pdf')
plt.close()


