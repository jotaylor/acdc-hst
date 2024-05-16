from math import * 
import numpy as np
from scipy.optimize import minimize
from scipy.special import factorial
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'legend.labelspacing':0.25, 'legend.fontsize': 12})
mpl.rcParams.update({'errorbar.capsize': 4})


## Exponential distribution
def exp_fun (E, k, Et):

	val = k * np.exp (-E / Et) #/ (Et *(1.0 - exp(-31.0/Et)) )

	return val

## Exponential + normal distribution
def exp_norm_fun (E, k, Et, w, mu, sigma):

	val = k *np.exp (-E / Et) + w * np.exp (-np.power(E - mu, 2.0) / 2./np.power(sigma,2.0)) #/ sqrt(2.*pi) / sigma

	return val

## Likelihood in form of C-statistics for exponential distribution
def Cstat (param, data):

	k, Et = param

	E = np.linspace (0,31, 32)


	model = exp_fun (E, k, Et)

	C = 2 * (model - data * np.log(model) + np.log(factorial(data)))

	C = np.sum(C)

	return C

## Likelihood in form of C-statistics for exponential + normal distribution
def Cstat_alt (param, data):

	k, Et, w, mu, sigma = param


	E = np.linspace (0,31, 32)

	model = exp_norm_fun (E, k, Et, w, mu, sigma)

	C = 2 * (model - data * np.log(model) + np.log(factorial(data)))

	C = np.sum(C)

	if isnan(C):

		C = 5e12
	


	return C


## Function which determines which of two models describe the spectra the best:
## (1) - exponential distribution, in this case the function returns two parameters: k, Et
## (2) - exponential + normal distribution, in this case the function returns five parameters: k, Et for exponential
## and w, mu and sigma for normal distribution

def study_spectra (data):

	x0 = (2, 10)

	res = minimize(Cstat, x0, method='L-BFGS-B', tol=1e-2, args=(data), bounds=((0, np.max(data)), (0.01, 60))   )

	Cstat_A = res.fun
	A_par   = res.x
	
	
	x0 = (2, 10, 20, 10, 3)

	res = minimize(Cstat_alt, x0, method='L-BFGS-B', tol=1e-2, args=(data), bounds=((0, np.max(data)), (0.01, 60), (0, np.max(data)), (3,31), (0.01,60)   ))

	Cstat_B = res.fun
	B_par   = res.x


	if Cstat_B < Cstat_A - 2*3: ##  AIC - model B has three more parameters

		return B_par

	else:	
	
		return A_par


if __name__ == "__main__":

	print ('Test: ',  factorial(0, exact=False))

	E = np.linspace (0,31, 32)

	k = 30.0

	Et = 7.0

	#data_pseudo = exp_fun (E, k, Et)

	data_pseudo = exp_norm_fun (E, k, Et, 15, 12, 5)

	for i in range (0, len(E)):

		data_pseudo[i] = int(data_pseudo[i])

		print (i, data_pseudo[i])


	
	

	val = study_spectra (data_pseudo)

	if len(val) < 3:

		model = exp_fun (E, val[0], val[1])

	else:

		model = exp_norm_fun (E, val[0], val[1], val[2], val[3], val[4])


	plt.plot (E, data_pseudo, 'k-', drawstyle='steps-mid', label='Pseudo data')
	plt.plot (E, model,       'b--', drawstyle='steps-mid',label='Model')
	plt.xlabel ('E')
	plt.ylabel ('Counts')
	plt.legend()
	plt.show()




