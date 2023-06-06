
# Visulise the distribution of probability for 1D integrator from t0 to t10
# Author: Jinyan Liu
# Inria Saclay and CMAP Ecole Polytechnique
# 2017

import numpy as np
import matplotlib.pyplot as plt

for i in range(1,141):
	filename = 'processLaw/stateDistribution/distributionProba_mode0.t' + str(141-i)
	data = np.genfromtxt(filename, dtype="float64", delimiter="\t")
	x = (data[:,0]-10)/10
	plt.figure('figure/proba_t' + str(141-i))
	plt.scatter(x, data[:,1])
	plt.plot(x, data[:,1])
	#plt.plot((0,0), (0,0.15), 'r--')
	#plt.xlim(-1,1)
	#plt.ylim(0,0.15)
	plt.xlabel('state')
	plt.ylabel('probability')

plt.show()



