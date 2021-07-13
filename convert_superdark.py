import sys
import asdf
from math import *
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})



af = asdf.open("superdark_FUVA_167_300days.asdf")


print (af.tree)

#sys.exit(0)
#print (af['pha4-5'])

ii = []
val = []

#superdark = np.zeros ((71, 865, 3))
superdark = np.zeros ((12, 289))


for i in range (3, 28):

	for j in range (0,142):

		for k in range (0, 1730):

			l = 'pha' + str(i) + '-' + str(i+1)

			#ii = int(i / 2.0)
			jj = int(j / 12.0)
			kk = int(k / 6.0)

			superdark [jj, kk] = superdark [jj, kk] + af[l][j,k]

			#print (i, j, k, ii, jj, kk, af[l][j,k], superdark [jj, kk, ii] )


			#if k > 10:
			#	sys.exit(0)


print ('Double check - total number of photons, superdark ', np.sum(superdark))

np.save('sp_FUVA_167_300days', superdark)

neg = plt.imshow(superdark, aspect='auto')
plt.colorbar(neg)
plt.show()


print ('Statistics of the selection')
print ('min ', np.min(superdark))
print ('mean', np.mean(superdark))
print ('std',  np.std(superdark))

#	print (np.mean (af[l]))
#	print (np.median(af[l]))

#	ii.append (i)
#	val.append (np.mean (af[l]))

#plt.plot (ii, val)
#plt.show()

#plt.imshow(af['pha10-11'])
#plt.show()



