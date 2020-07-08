import sqlite3
from sqlite3 import Error
import sys
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

## Connecting to the database
def create_connection(db_file):
	""" create a database connection to a SQLite database """
	conn = None
	try:
		conn = sqlite3.connect(db_file)
		print(sqlite3.version)
	except Error as e:
		print(e)
	return conn

def select_all_tasks(conn, hv, segment):
	cur = conn.cursor()
	cur.execute("SELECT id, latitude, longitude, expstart, region, solar_flux, dark_pha0, dark_pha1, dark_pha2, dark_pha3, \
                 dark_pha4,  dark_pha5,  dark_pha6,  dark_pha7,  dark_pha8,  dark_pha9,  dark_pha10, dark_pha11, \
                 dark_pha12, dark_pha13, dark_pha14, dark_pha15, dark_pha16, dark_pha17, dark_pha18, dark_pha19, \
                 dark_pha20, dark_pha21, dark_pha22, dark_pha23, dark_pha24, dark_pha25, dark_pha26, dark_pha27, \
                 dark_pha28, dark_pha29, dark_pha30, dark_pha31 \
                 FROM Darks  \
                 WHERE region='inner' and segment='"+str(segment)+"'  and hv="+str(hv)+' and id > 0 and id < 50')


	rows = cur.fetchall()

	data = np.zeros(32)

	for row in rows:

		print (row)

		for i in range (0, 32):
	
			data[i] = data[i] + float(row[i+6])

		#break

	return data



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

#	for i in range (0, len(C)):
#		print (i, C[i])

#	print ('.. ', np.sum(C), np.sum(C[2:29]))

#	sys.exit(0)	

	C = np.sum(C[2:29])

	

	#print ('---', k, Et, C)

	return C

## Likelihood in form of C-statistics for exponential + normal distribution
def Cstat_alt (param, data):

	k, Et, w, mu, sigma = param


	E = np.linspace (0,31, 32)

	model = exp_norm_fun (E, k, Et, w, mu, sigma)

	C = 2 * (model - data * np.log(model) + np.log(factorial(data)))

	C = np.sum(C[2:29])

	print ('---', k, Et, w, mu, sigma, C)

	if isnan(C):

		C = 5e12
	


	return C


## Function which determines which of two models describe the spectra the best:
## (1) - exponential distribution, in this case the function returns two parameters: k, Et
## (2) - exponential + normal distribution, in this case the function returns five parameters: k, Et for exponential
## and w, mu and sigma for normal distribution

def study_spectra (data):

	x0 = (np.max(data), 10)

	res = minimize(Cstat, x0, method='L-BFGS-B', tol=1e-6, args=(data), bounds=((0, 2.0*np.max(data)), (0.01, 60))   )

	Cstat_A = res.fun
	A_par   = res.x
	
	
	x0 = (A_par[0], A_par[1], A_par[0]/3.0, 10, 3)

	res = minimize(Cstat_alt, x0, method='L-BFGS-B', tol=1e-6, args=(data), bounds=((0, 2.0*np.max(data)), (0.01, 60), (0, np.max(data)), (3,31), (0.01,60)   ))

	Cstat_B = res.fun
	B_par   = res.x

	print ('Cstat for A: ', Cstat_A, ' Cstat for B: ', Cstat_B)

	if Cstat_B < Cstat_A - 2*3: ##  AIC - model B has three more parameters

		return B_par

	else:	
	
		return A_par


if __name__ == "__main__":

	hv = 167
	segment = 'FUVA'


	conn = create_connection('cos_dark.db')
	data = select_all_tasks(conn, hv, segment)
	
	print ('-------------------------------------')
	
	for i in range (0, len(data)):

		print  (i, data[i])

	print ('-------------------------------------')

	#sys.exit(0)

	#print ('Test: ',  factorial(0, exact=False))

	E = np.linspace (0,31, 32)

	#k = 30.0

	#Et = 7.0

	#data_pseudo = exp_fun (E, k, Et)

	#data_pseudo = exp_norm_fun (E, k, Et, 15, 12, 5)

	#for i in range (0, len(E)):

	#	data_pseudo[i] = int(data_pseudo[i])

	#	print (i, data_pseudo[i])


	
	

	val = study_spectra (data)

	print ('Result of optimisation: ', val)

	if len(val) < 3:

		model = exp_fun (E, val[0], val[1])

	else:

		model = exp_norm_fun (E, val[0], val[1], val[2], val[3], val[4])


	plt.plot (E, data, 'k-', drawstyle='steps-mid', label='Data')
	plt.plot (E, model,       'b--', drawstyle='steps-mid',label='Model')
	plt.xlabel ('E')
	plt.ylabel ('Counts')
	plt.legend()
	plt.show()




